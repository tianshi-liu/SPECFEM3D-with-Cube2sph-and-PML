#ifndef SEM_SYNC_IO

#include "config.h"

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int size,myrank;
    long total_size,offset;
    MPI_Comm comm;
    MPI_File fh;
    MPI_Request req;
    float *buf;
} WFAsyncIO;

// build a per-node filename "<filename>.<nodeid>", where nodeid is the
// MPI_COMM_WORLD rank of the node-local leader. This makes the file unique
// per node so it also works on a shared parallel filesystem.
static void
build_node_filename(MPI_Comm comm,const char *filename,
                    char *out,size_t out_size)
{
    int world_size,comm_size;
    MPI_Comm_size(MPI_COMM_WORLD,&world_size);
    MPI_Comm_size(comm,&comm_size);

    // single node: the node-local communicator spans all ranks, so there is
    // no risk of collision -> keep the plain filename
    if(comm_size == world_size) {
        snprintf(out,out_size,"%s",filename);
        return;
    }

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&world_rank);

    // node leader (local rank 0) broadcasts its world rank as the node id
    int nodeid = world_rank;
    MPI_Bcast(&nodeid,1,MPI_INT,0,comm);

    snprintf(out,out_size,"%s.node.%d",filename,nodeid);
}

void 
FC_FUNC_(open_subsample_write,OPEN_SUBSAMPLE_WRITE) (
    long *io_ptr,int *f_size, const char *filename
)
{
    long n = *f_size;
    WFAsyncIO *fwdio = (WFAsyncIO*) malloc(sizeof(WFAsyncIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;
    fwdio->buf = (float*) malloc(sizeof(float)*fwdio->size);

    // node-local communicator: ranks sharing the same node write to the
    // same node-local file, so offsets are computed within this group only
    MPI_Comm_split_type(
        MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,
        MPI_INFO_NULL,&fwdio->comm
    );

    // rank within the node-local communicator
    int myrank;
    MPI_Comm_rank(fwdio->comm,&myrank);
    fwdio->myrank = myrank;
    fwdio->req = MPI_REQUEST_NULL;

    // gather total size
    MPI_Exscan(
        &n,&fwdio->offset,1,
        MPI_LONG,MPI_SUM,fwdio->comm
    );
    if(myrank == 0) fwdio->offset = 0;

    // total size
    MPI_Allreduce(
        &n,&fwdio->total_size,1,
        MPI_LONG,MPI_SUM,fwdio->comm
    );

    // per-node filename so it also works on a shared parallel filesystem
    char filename_new[512];
    build_node_filename(fwdio->comm,filename,filename_new,sizeof(filename_new));

    // open file
    MPI_File_open(
        fwdio->comm,filename_new,
        MPI_MODE_WRONLY + MPI_MODE_CREATE,
        MPI_INFO_NULL,&fwdio->fh
    );
}

void 
FC_FUNC_(open_subsample_read,OPEN_SUBSAMPLE_READ)(
    long *io_ptr,int *f_size,const char *filename,
    int *nstep, int *nstep_dump
)
{
    long n = *f_size;
    WFAsyncIO *fwdio = (WFAsyncIO*) malloc(sizeof(WFAsyncIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;
    fwdio->buf = (float*) malloc(sizeof(float)*fwdio->size);

    // node-local communicator: ranks sharing the same node read from the
    // same node-local file, so offsets are computed within this group only
    MPI_Comm_split_type(
        MPI_COMM_WORLD,MPI_COMM_TYPE_SHARED,0,
        MPI_INFO_NULL,&fwdio->comm
    );

    // rank within the node-local communicator
    int myrank = 0;
    MPI_Comm_rank(fwdio->comm,&myrank);
    fwdio->myrank = myrank;
    fwdio->req = MPI_REQUEST_NULL;

    // gather total size
    MPI_Exscan(
        &n,&fwdio->offset,1,
        MPI_LONG,MPI_SUM,fwdio->comm
    );
    if(myrank == 0) fwdio->offset = 0;

    // total size
    MPI_Allreduce(
        &n,&fwdio->total_size,1,
        MPI_LONG,MPI_SUM,fwdio->comm
    );

    // per-node filename so it also works on a shared parallel filesystem
    char filename_new[512];
    build_node_filename(fwdio->comm,filename,filename_new,sizeof(filename_new));

    // open file
    MPI_File_open(
        fwdio->comm,filename_new,
        MPI_MODE_RDONLY,
        MPI_INFO_NULL,&fwdio->fh
    );

    // pre read the first slice into buffer
    int i_save = *nstep / *nstep_dump;
    MPI_Offset file_offset = (i_save - 1) * 
                             fwdio->total_size + 
                             fwdio->offset;
    file_offset = file_offset * sizeof(float);
    MPI_File_iread_at(
        fwdio->fh,
        file_offset,
        fwdio->buf,fwdio->size,
        MPI_FLOAT,&fwdio->req
    );
}

void 
FC_FUNC_(close_subsample_file,CLOSE_SUBSAMPLE_FILE) (long *io_ptr)
{
    // convert to fwdio
    WFAsyncIO *fwdio = (WFAsyncIO *)(*io_ptr);

    // WAIT job finish
    MPI_Wait(&fwdio->req,MPI_STATUS_IGNORE);

    // closefile
    MPI_File_close(&fwdio->fh);

    // free node-local communicator
    MPI_Comm_free(&fwdio->comm);

    if(fwdio->buf != NULL) {
        free(fwdio->buf);
        fwdio->buf = NULL;
    }

    free(fwdio);
}

void 
FC_FUNC_(write_subsample_file,WRITE_SUBSAMPLE_FILE)(
    long *io_ptr,int *f_it_save,const float *displ
)
{
    WFAsyncIO *fwdio = (WFAsyncIO *) (*io_ptr);
    int i_save = *f_it_save;
    MPI_Offset file_offset = (i_save - 1) * 
                             fwdio->total_size + 
                             fwdio->offset;
    file_offset = file_offset * sizeof(float);

    // block request
    MPI_Wait(&fwdio->req,MPI_STATUS_IGNORE);

    // copy data to buffer
    memcpy(fwdio->buf,displ,sizeof(float)*fwdio->size);
    
    MPI_File_iwrite_at(
        fwdio->fh,
        file_offset,
        fwdio->buf,fwdio->size,
        MPI_FLOAT,&fwdio->req
    );
}

void 
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    long *io_ptr,int *f_it_save, float* displ
)
{
    WFAsyncIO *fwdio = (WFAsyncIO*) (*io_ptr);

    // wait job to finish
    MPI_Wait(&fwdio->req,MPI_STATUS_IGNORE);
    memcpy(displ,fwdio->buf,sizeof(float)*fwdio->size);

    // read previous buffer if required
    int i_save = *f_it_save - 1;
    if(i_save >=1) {
        MPI_Offset file_offset = (i_save - 1) * 
                                fwdio->total_size + 
                                fwdio->offset;
        file_offset = file_offset * sizeof(float);
        
        MPI_File_iread_at(
            fwdio->fh,
            file_offset,
            fwdio->buf,fwdio->size,
            MPI_FLOAT,&fwdio->req
        );
    }
}

#endif