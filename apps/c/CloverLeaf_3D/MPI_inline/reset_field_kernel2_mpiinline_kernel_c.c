//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_reset_field_kernel2;
int ydim0_reset_field_kernel2;
int xdim1_reset_field_kernel2;
int ydim1_reset_field_kernel2;
int xdim2_reset_field_kernel2;
int ydim2_reset_field_kernel2;
int xdim3_reset_field_kernel2;
int ydim3_reset_field_kernel2;
int xdim4_reset_field_kernel2;
int ydim4_reset_field_kernel2;
int xdim5_reset_field_kernel2;
int ydim5_reset_field_kernel2;


#define OPS_ACC0(x,y,z) (n_x*1+n_y*xdim0_reset_field_kernel2*1+n_z*xdim0_reset_field_kernel2*ydim0_reset_field_kernel2*1+x+xdim0_reset_field_kernel2*(y)+xdim0_reset_field_kernel2*ydim0_reset_field_kernel2*(z))
#define OPS_ACC1(x,y,z) (n_x*1+n_y*xdim1_reset_field_kernel2*1+n_z*xdim1_reset_field_kernel2*ydim1_reset_field_kernel2*1+x+xdim1_reset_field_kernel2*(y)+xdim1_reset_field_kernel2*ydim1_reset_field_kernel2*(z))
#define OPS_ACC2(x,y,z) (n_x*1+n_y*xdim2_reset_field_kernel2*1+n_z*xdim2_reset_field_kernel2*ydim2_reset_field_kernel2*1+x+xdim2_reset_field_kernel2*(y)+xdim2_reset_field_kernel2*ydim2_reset_field_kernel2*(z))
#define OPS_ACC3(x,y,z) (n_x*1+n_y*xdim3_reset_field_kernel2*1+n_z*xdim3_reset_field_kernel2*ydim3_reset_field_kernel2*1+x+xdim3_reset_field_kernel2*(y)+xdim3_reset_field_kernel2*ydim3_reset_field_kernel2*(z))
#define OPS_ACC4(x,y,z) (n_x*1+n_y*xdim4_reset_field_kernel2*1+n_z*xdim4_reset_field_kernel2*ydim4_reset_field_kernel2*1+x+xdim4_reset_field_kernel2*(y)+xdim4_reset_field_kernel2*ydim4_reset_field_kernel2*(z))
#define OPS_ACC5(x,y,z) (n_x*1+n_y*xdim5_reset_field_kernel2*1+n_z*xdim5_reset_field_kernel2*ydim5_reset_field_kernel2*1+x+xdim5_reset_field_kernel2*(y)+xdim5_reset_field_kernel2*ydim5_reset_field_kernel2*(z))

//user function



void reset_field_kernel2_c_wrapper(
  double * restrict xvel0,
  const double * restrict xvel1,
  double * restrict yvel0,
  const double * restrict yvel1,
  double * restrict zvel0,
  const double * restrict zvel1,
  int x_size, int y_size, int z_size) {
  #pragma omp parallel for
  for ( int n_z=0; n_z<z_size; n_z++ ){
    for ( int n_y=0; n_y<y_size; n_y++ ){
      for ( int n_x=0; n_x<x_size; n_x++ ){
        

  xvel0[OPS_ACC0(0,0,0)]  = xvel1[OPS_ACC1(0,0,0)] ;
  yvel0[OPS_ACC2(0,0,0)]  = yvel1[OPS_ACC3(0,0,0)] ;
  zvel0[OPS_ACC4(0,0,0)]  = zvel1[OPS_ACC5(0,0,0)] ;

      }
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5

