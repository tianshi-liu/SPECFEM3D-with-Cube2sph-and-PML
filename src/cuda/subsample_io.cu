#include "mesh_constants_cuda.h"
#include "config.h"

extern "C" {

void 
FC_FUNC_(write_subsample_file,WRITE_SUBSAMPLE_FILE)(
    int *f_it_save,const float *displ
);

void
FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE) (
    int *f_it_save, float* displ
);

}


extern "C" void 
FC_FUNC_(write_subsample_file_cuda,WRITE_SUBSAMPLE_FILE_cuda)(
    int *f_it_save, long *Mesh_pointer
)
{
    Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper
    cudaMemcpy(
        mp->h_sub_buffer,
        mp->d_displ,mp->NGLOB_AB*3*sizeof(realw),
        cudaMemcpyDeviceToHost
    );

    FC_FUNC_(write_subsample_file,WRITE_SUBSAMPLE_FILE)(
        f_it_save,mp->h_sub_buffer
    );
}

extern "C" void 
FC_FUNC_(read_subsample_file_cuda,READ_SUBSAMPLE_FILE_cuda)(
    int *f_it_save, long *Mesh_pointer
)
{
    Mesh* mp = (Mesh*)(*Mesh_pointer); // get Mesh from fortran integer wrapper

    FC_FUNC_(read_subsample_file,READ_SUBSAMPLE_FILE)(
        f_it_save,mp->h_sub_buffer
    );

    cudaMemcpy(
        mp->d_b_displ,
        mp->h_sub_buffer,mp->NGLOB_AB*3*sizeof(realw),
        cudaMemcpyHostToDevice
    );


}