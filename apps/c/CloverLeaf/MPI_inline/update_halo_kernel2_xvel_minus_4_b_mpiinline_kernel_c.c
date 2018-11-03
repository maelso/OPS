//
// auto-generated by ops.py
//

int xdim0_update_halo_kernel2_xvel_minus_4_b;
int xdim1_update_halo_kernel2_xvel_minus_4_b;


//user function



void update_halo_kernel2_xvel_minus_4_b_c_wrapper(
  double * restrict xvel0_p,
  double * restrict xvel1_p,
  const int * restrict fields,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      ptr_double xvel0 = { xvel0_p + n_x*1 + n_y * xdim0_update_halo_kernel2_xvel_minus_4_b*1, xdim0_update_halo_kernel2_xvel_minus_4_b};
      ptr_double xvel1 = { xvel1_p + n_x*1 + n_y * xdim1_update_halo_kernel2_xvel_minus_4_b*1, xdim1_update_halo_kernel2_xvel_minus_4_b};
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1

