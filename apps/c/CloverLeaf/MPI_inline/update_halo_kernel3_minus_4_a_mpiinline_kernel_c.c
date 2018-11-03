//
// auto-generated by ops.py
//

int xdim0_update_halo_kernel3_minus_4_a;
int xdim1_update_halo_kernel3_minus_4_a;


//user function



void update_halo_kernel3_minus_4_a_c_wrapper(
  double * restrict vol_flux_x_p,
  double * restrict mass_flux_x_p,
  const int * restrict fields,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      ptr_double vol_flux_x = { vol_flux_x_p + n_x*1 + n_y * xdim0_update_halo_kernel3_minus_4_a*1, xdim0_update_halo_kernel3_minus_4_a};
      ptr_double mass_flux_x = { mass_flux_x_p + n_x*1 + n_y * xdim1_update_halo_kernel3_minus_4_a*1, xdim1_update_halo_kernel3_minus_4_a};
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1

