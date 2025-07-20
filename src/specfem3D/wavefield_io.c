#include "config.h"

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int size,myrank;
    long total_size,offset;
    MPI_File fh;
    MPI_Request req;
    float *buf;
} SubSampleIO ;

SubSampleIO *subio;

void 
FC_FUNC_(open_subsample_write,OPEN_SUBSAMPLE_WRITE) (
    int *f_size, const char *filename
)
{
    long n = *f_size;
    subio = (SubSampleIO*) malloc(sizeof(SubSampleIO));
    subio->size = *f_size;
    subio->buf = (float*) malloc(sizeof(float)*subio->size);

    // mpisize
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    subio->myrank = myrank;
    subio->req = MPI_REQUEST_NULL;

    // gather total size
    MPI_Exscan(
        &n,&subio->offset,1,
        MPI_LONG,MPI_SUM,MPI_COMM_WORLD
    );
    if(myrank == 0) subio->offset = 0;

    // total size
    MPI_Allreduce(
        &n,&subio->total_size,1,
        MPI_LONG,MPI_SUM,MPI_COMM_WORLD
    );

    // open file
    MPI_File_open(
        MPI_COMM_WORLD,filename,
        MPI_MODE_WRONLY + MPI_MODE_CREATE,
        MPI_INFO_NULL,&subio->fh
    );
}

void 
FC_FUNC_(open_subsample_read,OPEN_SUBSAMPLE_READ)(
    int *f_size,const char *filename,
    int *nstep, int *nstep_dump
)
{
    long n = *f_size;
    subio = (SubSampleIO*) malloc(sizeof(SubSampleIO));
    subio->size = *f_size;
    subio->buf = (float*) malloc(sizeof(float)*subio->size);

    // mpisize
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    subio->myrank = myrank;
    subio->req = MPI_REQUEST_NULL;

    // gather total size
    MPI_Exscan(
        &n,&subio->offset,1,
        MPI_LONG,MPI_SUM,MPI_COMM_WORLD
    );
    if(myrank == 0) subio->offset = 0;

    // total size
    MPI_Allreduce(
        &n,&subio->total_size,1,
        MPI_LONG,MPI_SUM,MPI_COMM_WORLD
    );

    // open file
    MPI_File_open(
        MPI_COMM_WORLD,filename,
        MPI_MODE_RDONLY,
        MPI_INFO_NULL,&subio->fh
    );

    // pre read the first slice into buffer
    int i_save = *nstep / *nstep_dump;
    MPI_Offset file_offset = (i_save - 1) * 
                             subio->total_size + 
                             subio->offset;
    file_offset = file_offset * sizeof(float);
    MPI_File_iread_at(
        subio->fh,
        file_offset,
        subio->buf,subio->size,
        MPI_FLOAT,&subio->req
    );
}

void 
FC_FUNC_(close_subsample_file,CLOSE_SUBSAMPLE_FILE) ()
{
    // WAIT job finish
    MPI_Wait(&subio->req,MPI_STATUS_IGNORE);

    // closefile
    MPI_File_close(&subio->fh);

    if(subio->buf != NULL) {
        free(subio->buf);
        subio->buf = NULL;
    }

    free(subio);
}


void 
FC_FUNC_(write_subsample_file,WRITE_SUBSAMPLE_FILE)(
    int *f_it_save,const float *displ
)
{
    int i_save = *f_it_save;
    MPI_Offset file_offset = (i_save - 1) * 
                             subio->total_size + 
                             subio->offset;
    file_offset = file_offset * sizeof(float);

    // block request
    MPI_Wait(&subio->req,MPI_STATUS_IGNORE);

    // copy data to buffer
    memcpy(subio->buf,displ,sizeof(float)*subio->size);
    
    MPI_File_iwrite_at(
        subio->fh,
        file_offset,
        subio->buf,subio->size,
        MPI_FLOAT,&subio->req
    );
}

void 
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    int *f_it_save, float* displ
)
{
    // wait job to finish
    MPI_Wait(&subio->req,MPI_STATUS_IGNORE);
    memcpy(displ,subio->buf,sizeof(float)*subio->size);

    // read previous buffer if required
    int i_save = *f_it_save - 1;
    if(i_save >=1) {
        MPI_Offset file_offset = (i_save - 1) * 
                                subio->total_size + 
                                subio->offset;
        file_offset = file_offset * sizeof(float);
        
        MPI_File_iread_at(
            subio->fh,
            file_offset,
            subio->buf,subio->size,
            MPI_FLOAT,&subio->req
        );
    }
}
