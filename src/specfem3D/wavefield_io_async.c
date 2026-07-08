#ifndef SEM_SYNC_IO

#include "config.h"

#include <mpi.h>

#include <pthread.h>
#include <fcntl.h>
#include <unistd.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

// kind of the pending background request
typedef enum { IO_NONE = 0, IO_WRITE, IO_READ } io_op_t;

// Thread-based asynchronous IO. Instead of using MPI-IO non-blocking calls to
// overlap file access with computation, a dedicated worker thread performs the
// pread()/pwrite() while the main (compute) thread continues. Each MPI rank
// writes to its own private file "<filename>.<worldrank>", so there is no
// cross-rank offset bookkeeping and no collective communication at all.
typedef struct {
    int  size;          // number of floats per record
    int  myrank;        // MPI_COMM_WORLD rank -> used for the filename
    int  fd;            // POSIX file descriptor
    float *buf;         // internal staging buffer (mirrors the MPI version)

    // background worker plumbing
    pthread_t       worker;
    pthread_mutex_t lock;
    pthread_cond_t  cv_submit;   // main -> worker : a request is ready
    pthread_cond_t  cv_done;     // worker -> main : request finished
    io_op_t         op;          // pending op (IO_NONE when idle)
    off_t           offset;      // byte offset of the pending op
    int             busy;        // 1 while a request is in flight
    int             shutdown;    // 1 asks the worker to exit
} WFThreadsIO;

// ---------------------------------------------------------------------------
// full pread/pwrite helpers (handle short transfers)
// ---------------------------------------------------------------------------
static void
pwrite_full(int fd,const void *buf,size_t nbytes,off_t off)
{
    const char *p = (const char*)buf;
    while(nbytes > 0) {
        ssize_t n = pwrite(fd,p,nbytes,off);
        if(n < 0) {
            if(errno == EINTR) continue;
            fprintf(stderr,"pwrite failed: %s\n",strerror(errno));
            return;
        }
        p      += n;
        off    += n;
        nbytes -= (size_t)n;
    }
}

static void
pread_full(int fd,void *buf,size_t nbytes,off_t off)
{
    char *p = (char*)buf;
    while(nbytes > 0) {
        ssize_t n = pread(fd,p,nbytes,off);
        if(n < 0) {
            if(errno == EINTR) continue;
            fprintf(stderr,"pread failed: %s\n",strerror(errno));
            return;
        }
        if(n == 0) break;            // EOF (reading past end-of-file)
        p      += n;
        off    += n;
        nbytes -= (size_t)n;
    }
}

// ---------------------------------------------------------------------------
// worker thread: sleeps until a request is submitted, performs it, repeats
// ---------------------------------------------------------------------------
static void *
io_worker(void *arg)
{
    WFThreadsIO *io = (WFThreadsIO*)arg;
    size_t nbytes = (size_t)io->size * sizeof(float);

    pthread_mutex_lock(&io->lock);
    for(;;) {
        while(io->op == IO_NONE && !io->shutdown)
            pthread_cond_wait(&io->cv_submit,&io->lock);

        if(io->op == IO_NONE && io->shutdown) break;

        io_op_t op  = io->op;
        off_t   off = io->offset;

        // do the (potentially slow) IO without holding the lock so the main
        // thread can keep computing
        pthread_mutex_unlock(&io->lock);
        if(op == IO_WRITE) pwrite_full(io->fd,io->buf,nbytes,off);
        else               pread_full (io->fd,io->buf,nbytes,off);
        pthread_mutex_lock(&io->lock);

        io->op   = IO_NONE;
        io->busy = 0;
        pthread_cond_signal(&io->cv_done);
    }
    pthread_mutex_unlock(&io->lock);
    return NULL;
}

// block until the in-flight request (if any) has finished
static void
io_wait(WFThreadsIO *io)
{
    pthread_mutex_lock(&io->lock);
    while(io->busy) pthread_cond_wait(&io->cv_done,&io->lock);
    pthread_mutex_unlock(&io->lock);
}

// hand a new request to the worker (caller must have io_wait'ed first)
static void
io_submit(WFThreadsIO *io,io_op_t op,off_t off)
{
    pthread_mutex_lock(&io->lock);
    io->op     = op;
    io->offset = off;
    io->busy   = 1;
    pthread_cond_signal(&io->cv_submit);
    pthread_mutex_unlock(&io->lock);
}

// spin up the worker + init synchronisation primitives
static void
io_start(WFThreadsIO *io)
{
    pthread_mutex_init(&io->lock,NULL);
    pthread_cond_init(&io->cv_submit,NULL);
    pthread_cond_init(&io->cv_done,NULL);
    io->op       = IO_NONE;
    io->busy     = 0;
    io->shutdown = 0;
    pthread_create(&io->worker,NULL,io_worker,io);
}

// build the per-rank filename "<filename>.<worldrank>"
static void
build_rank_filename(int rank,const char *filename,char *out,size_t out_size)
{
    snprintf(out,out_size,"%s.%06d",filename,rank);
}

// ---------------------------------------------------------------------------
// public API (mirrors asyncio.c)
// ---------------------------------------------------------------------------
void
FC_FUNC_(open_subsample_write,OPEN_SUBSAMPLE_WRITE) (
    long *io_ptr,int *f_size,const char *filename
)
{
    WFThreadsIO *fwdio = (WFThreadsIO*) malloc(sizeof(WFThreadsIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;
    fwdio->buf  = (float*) malloc(sizeof(float)*fwdio->size);

    MPI_Comm_rank(MPI_COMM_WORLD,&fwdio->myrank);

    // per-rank private file
    char filename_new[512];
    build_rank_filename(fwdio->myrank,filename,filename_new,sizeof(filename_new));

    fwdio->fd = open(filename_new,O_WRONLY | O_CREAT,0644);
    if(fwdio->fd < 0)
        fprintf(stderr,"cannot open %s for write: %s\n",
                filename_new,strerror(errno));

    io_start(fwdio);
}

void
FC_FUNC_(open_subsample_read,OPEN_SUBSAMPLE_READ)(
    long *io_ptr,int *f_size,const char *filename,
    int *nstep,int *nstep_dump
)
{
    WFThreadsIO *fwdio = (WFThreadsIO*) malloc(sizeof(WFThreadsIO));
    *io_ptr = (long)fwdio;
    fwdio->size = *f_size;
    fwdio->buf  = (float*) malloc(sizeof(float)*fwdio->size);

    MPI_Comm_rank(MPI_COMM_WORLD,&fwdio->myrank);

    // per-rank private file
    char filename_new[512];
    build_rank_filename(fwdio->myrank,filename,filename_new,sizeof(filename_new));

    fwdio->fd = open(filename_new,O_RDONLY);
    if(fwdio->fd < 0)
        fprintf(stderr,"cannot open %s for read: %s\n",
                filename_new,strerror(errno));

    io_start(fwdio);

    // pre-read the first slice into the buffer. With a private per-rank file
    // the stride between records is just the record size.
    int i_save = *nstep / *nstep_dump;
    off_t file_offset = (off_t)(i_save - 1) * fwdio->size * sizeof(float);
    io_submit(fwdio,IO_READ,file_offset);
}

void
FC_FUNC_(close_subsample_file,CLOSE_SUBSAMPLE_FILE) (long *io_ptr)
{
    WFThreadsIO *fwdio = (WFThreadsIO *)(*io_ptr);

    // wait for the outstanding request to finish
    io_wait(fwdio);

    // ask the worker to exit and join it
    pthread_mutex_lock(&fwdio->lock);
    fwdio->shutdown = 1;
    pthread_cond_signal(&fwdio->cv_submit);
    pthread_mutex_unlock(&fwdio->lock);
    pthread_join(fwdio->worker,NULL);

    pthread_mutex_destroy(&fwdio->lock);
    pthread_cond_destroy(&fwdio->cv_submit);
    pthread_cond_destroy(&fwdio->cv_done);

    close(fwdio->fd);

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
    WFThreadsIO *fwdio = (WFThreadsIO *) (*io_ptr);
    int i_save = *f_it_save;
    off_t file_offset = (off_t)(i_save - 1) * fwdio->size * sizeof(float);

    // block until the previous write drained out of the shared buffer
    io_wait(fwdio);

    // stage the new data, then let the worker flush it in the background
    memcpy(fwdio->buf,displ,sizeof(float)*fwdio->size);
    io_submit(fwdio,IO_WRITE,file_offset);
}

void
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    long *io_ptr,int *f_it_save,float *displ
)
{
    WFThreadsIO *fwdio = (WFThreadsIO*) (*io_ptr);

    // wait for the prefetch to complete, then consume it
    io_wait(fwdio);
    memcpy(displ,fwdio->buf,sizeof(float)*fwdio->size);

    // prefetch the previous slice (data is consumed in reverse, as in asyncio.c)
    int i_save = *f_it_save - 1;
    if(i_save >= 1) {
        off_t file_offset = (off_t)(i_save - 1) * fwdio->size * sizeof(float);
        io_submit(fwdio,IO_READ,file_offset);
    }
}

#endif
