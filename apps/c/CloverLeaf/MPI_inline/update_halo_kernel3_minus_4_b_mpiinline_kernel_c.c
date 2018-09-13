//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_update_halo_kernel3_minus_4_b;
int xdim1_update_halo_kernel3_minus_4_b;


#define OPS_ACC0(x,y) (n_x*1+n_y*xdim0_update_halo_kernel3_minus_4_b*1+x+xdim0_update_halo_kernel3_minus_4_b*(y))
#define OPS_ACC1(x,y) (n_x*1+n_y*xdim1_update_halo_kernel3_minus_4_b*1+x+xdim1_update_halo_kernel3_minus_4_b*(y))

//user function



void update_halo_kernel3_minus_4_b_c_wrapper(
  double * restrict vol_flux_x,
  double * restrict mass_flux_x,
  const int * restrict fields,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      
  if(fields[FIELD_VOL_FLUX_X] == 1)  vol_flux_x[OPS_ACC0(0,0)]  = -(vol_flux_x[OPS_ACC0(-4,0)]);
  if(fields[FIELD_MASS_FLUX_X] == 1) mass_flux_x[OPS_ACC1(0,0)] = -(mass_flux_x[OPS_ACC1(-4,0)]);

    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1

