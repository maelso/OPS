//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_update_halo_kernel2_xvel_plus_4_top;
int ydim0_update_halo_kernel2_xvel_plus_4_top;
int xdim1_update_halo_kernel2_xvel_plus_4_top;
int ydim1_update_halo_kernel2_xvel_plus_4_top;


#define OPS_ACC0(x,y,z) (n_x*1+n_y*xdim0_update_halo_kernel2_xvel_plus_4_top*1+n_z*xdim0_update_halo_kernel2_xvel_plus_4_top*ydim0_update_halo_kernel2_xvel_plus_4_top*1+x+xdim0_update_halo_kernel2_xvel_plus_4_top*(y)+xdim0_update_halo_kernel2_xvel_plus_4_top*ydim0_update_halo_kernel2_xvel_plus_4_top*(z))
#define OPS_ACC1(x,y,z) (n_x*1+n_y*xdim1_update_halo_kernel2_xvel_plus_4_top*1+n_z*xdim1_update_halo_kernel2_xvel_plus_4_top*ydim1_update_halo_kernel2_xvel_plus_4_top*1+x+xdim1_update_halo_kernel2_xvel_plus_4_top*(y)+xdim1_update_halo_kernel2_xvel_plus_4_top*ydim1_update_halo_kernel2_xvel_plus_4_top*(z))

//user function



void update_halo_kernel2_xvel_plus_4_top_c_wrapper(
  double * restrict xvel0,
  double * restrict xvel1,
  const int * restrict fields,
  int x_size, int y_size, int z_size) {
  #pragma omp parallel for
  for ( int n_z=0; n_z<z_size; n_z++ ){
    for ( int n_y=0; n_y<y_size; n_y++ ){
      for ( int n_x=0; n_x<x_size; n_x++ ){
        
  if(fields[FIELD_XVEL0] == 1) xvel0[OPS_ACC0(0,0,0)] = xvel0[OPS_ACC0(0,-4,0)];
  if(fields[FIELD_XVEL1] == 1) xvel1[OPS_ACC1(0,0,0)] = xvel1[OPS_ACC1(0,-4,0)];

      }
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1

