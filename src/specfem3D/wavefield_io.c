#ifndef SEM_ASYNC_IO

#include "config.h"

#include <mpi.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int size,myrank;
    FILE *fp;
} WFsyncIO;

void 
FC_FUNC_(open_subsample_write,OPEN_SUBSAMPLE_WRITE) (
    long *io_ptr,int *f_size, const char *filename
)
{
    long n = *f_size;
    WFsyncIO *fwdio = (WFsyncIO*) malloc(sizeof(WFsyncIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;

    // mpisize
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    fwdio->myrank = myrank;

    // create new filename
    size_t len_f = strlen(filename);
    char filename_new[512];
    snprintf(filename_new,sizeof(filename_new),"%s.%06d",filename,myrank); 

    // open file
    fwdio -> fp = fopen(filename_new,"wb");
}

void 
FC_FUNC_(open_subsample_read,OPEN_SUBSAMPLE_READ)(
    long *io_ptr,int *f_size,const char *filename,
    int *nstep, int *nstep_dump
)
{
    long n = *f_size;
    WFsyncIO *fwdio = (WFsyncIO*) malloc(sizeof(WFsyncIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;

    // mpisize
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    fwdio->myrank = myrank;

    // create new filename
    size_t len_f = strlen(filename);
    char filename_new[512];
    snprintf(filename_new,sizeof(filename_new),"%s.%06d",filename,myrank); 

    // open file
    fwdio -> fp = fopen(filename_new,"rb");
}

void 
FC_FUNC_(close_subsample_file,CLOSE_SUBSAMPLE_FILE) (long *io_ptr)
{
    // convert to fwdio
    WFsyncIO *fwdio = (WFsyncIO *)(*io_ptr);

    // closefile
    fclose(fwdio->fp);

    free(fwdio);
}

void 
FC_FUNC_(write_subsample_file,WRITE_SUBSAMPLE_FILE)(
    long *io_ptr,int *f_it_save,const float *displ
)
{
    WFsyncIO *fwdio = (WFsyncIO *) (*io_ptr);

    fwrite(displ,sizeof(float),fwdio->size,fwdio->fp);
}

void 
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    long *io_ptr,int *f_it_save, float* displ
)
{
    WFsyncIO *fwdio = (WFsyncIO*) (*io_ptr);

    // get current loation
    int i_save = *f_it_save - 1;
    long offset = i_save * fwdio->size * sizeof(float);

    fseek(fwdio->fp,offset,SEEK_SET);
    size_t nread = fread(displ,sizeof(float),fwdio->size,fwdio->fp);
    if (nread != fwdio->size) {
        if (feof(fwdio->fp)) {
            fprintf(stderr, "Unexpected end of file: got %zu of %d elements\n",
                    nread, fwdio->size);
        } else if (ferror(fwdio->fp)) {
            perror("fread");
        }
    }
}

#endif