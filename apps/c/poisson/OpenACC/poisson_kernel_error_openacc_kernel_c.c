//
// auto-generated by ops.py
//
#include "./OpenACC/poisson_common.h"

#define OPS_GPU

int xdim0_poisson_kernel_error;
int xdim1_poisson_kernel_error;


#undef OPS_ACC0
#undef OPS_ACC1


#define OPS_ACC0(x,y) (x+xdim0_poisson_kernel_error*(y))
#define OPS_ACC1(x,y) (x+xdim1_poisson_kernel_error*(y))

//user function
inline
void poisson_kernel_error(const double *u, const double *ref, double *err) {
  *err = *err + (u[OPS_ACC0(0,0)]-ref[OPS_ACC1(0,0)])*(u[OPS_ACC0(0,0)]-ref[OPS_ACC1(0,0)]);
}


#undef OPS_ACC0
#undef OPS_ACC1



void poisson_kernel_error_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  int x_size, int y_size) {
  double p_a2_0 = p_a2[0];
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1) reduction(+:p_a2_0)
  #pragma acc loop reduction(+:p_a2_0)
  #endif
  for ( int n_y=0; n_y<y_size; n_y++ ){
    #ifdef OPS_GPU
    #pragma acc loop reduction(+:p_a2_0)
    #endif
    for ( int n_x=0; n_x<x_size; n_x++ ){
      poisson_kernel_error(  p_a0 + n_x*1*1 + n_y*xdim0_poisson_kernel_error*1*1,
           p_a1 + n_x*1*1 + n_y*xdim1_poisson_kernel_error*1*1, &p_a2_0 );

    }
  }
  p_a2[0] = p_a2_0;
}
