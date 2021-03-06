//
// auto-generated by ops.py
//

int xdim0_poisson_kernel_initialguess;


//user function



void poisson_kernel_initialguess_c_wrapper(
  double * restrict u_p,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      ptr_double u = { u_p + n_x*1 + n_y * xdim0_poisson_kernel_initialguess*1, xdim0_poisson_kernel_initialguess};

      OPS_ACC(u, 0, 0) = 0.0;
    }
  }
}
