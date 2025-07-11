#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"

typedef struct {
    int size,myrank;
    long total_size,offset;
    MPI_File fh;
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

    // mpisize
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    subio -> myrank = myrank;

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

    // printf("myrank = %d %d %ld %ld\n",myrank,*f_size,subio->offset,subio->total_size); 
    

    // open file
    MPI_File_open(
        MPI_COMM_WORLD,filename,
        MPI_MODE_WRONLY + MPI_MODE_CREATE,
        MPI_INFO_NULL,&subio->fh
    );
}

void 
FC_FUNC_(open_subsample_read,OPEN_SUBSAMPLE_READ)(
    int *f_size,const char *filename
)
{
    long n = *f_size;
    subio = (SubSampleIO*) malloc(sizeof(SubSampleIO));
    subio->size = *f_size;

    // mpisize
    int myrank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    subio -> myrank = myrank;

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
}

void 
FC_FUNC_(close_subsample_file,CLOSE_SUBSAMPLE_FILE) ()
{
    MPI_File_close(&subio->fh);

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
    
    int ierr = 
    MPI_File_write_at(
        subio->fh,
        file_offset,
        displ,subio->size,
        MPI_FLOAT,MPI_STATUS_IGNORE
    );

    if(ierr != MPI_SUCCESS) {
        char err_str[MPI_MAX_ERROR_STRING];
        int len;
        MPI_Error_string(ierr,err_str,&len);
        printf("MPI error: %s\n", err_str);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
}

void 
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    int *f_it_save, float* displ
)
{
    int i_save = *f_it_save;
    MPI_Offset file_offset = (i_save - 1) * 
                             subio->total_size + 
                             subio->offset;
    file_offset = file_offset * sizeof(float);

    MPI_File_read_at(
        subio->fh,
        file_offset,
        displ,subio->size,
        MPI_FLOAT,MPI_STATUS_IGNORE
    );
}
