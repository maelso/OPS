//
// auto-generated by ops.py//

//header
#define OPS_3D
#define OPS_SOA
#define OPS_ACC_MD_MACROS
#include "ops_lib_cpp.h"
#ifdef OPS_MPI
#include "ops_mpi_core.h"
#endif

//set max number of OMP threads for reductions
#ifndef MAX_REDUCT_THREADS
#define MAX_REDUCT_THREADS 64
#endif
//global constants


void ops_init_backend() {}

//user kernel files
#include "multidim_kernel_omp_kernel.cpp"
#include "multidim_copy_kernel_omp_kernel.cpp"
#include "multidim_reduce_kernel_omp_kernel.cpp"
