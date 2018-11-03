//
// auto-generated by ops.py
//

int xdim0_advec_mom_kernel1_y_nonvector;
int xdim1_advec_mom_kernel1_y_nonvector;
int xdim2_advec_mom_kernel1_y_nonvector;
int xdim3_advec_mom_kernel1_y_nonvector;
int xdim4_advec_mom_kernel1_y_nonvector;


//user function



void advec_mom_kernel1_y_nonvector_c_wrapper(
  double * restrict node_flux_p,
  double * restrict node_mass_pre_p,
  double * restrict mom_flux_p,
  double * restrict celldy_p,
  double * restrict vel1_p,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      const ptr_double node_flux = { node_flux_p + n_x*1 + n_y * xdim0_advec_mom_kernel1_y_nonvector*1, xdim0_advec_mom_kernel1_y_nonvector};
      const ptr_double node_mass_pre = { node_mass_pre_p + n_x*1 + n_y * xdim1_advec_mom_kernel1_y_nonvector*1, xdim1_advec_mom_kernel1_y_nonvector};
      ptr_double mom_flux = { mom_flux_p + n_x*1 + n_y * xdim2_advec_mom_kernel1_y_nonvector*1, xdim2_advec_mom_kernel1_y_nonvector};
      const ptr_double celldy = { celldy_p + n_x*0 + n_y * xdim3_advec_mom_kernel1_y_nonvector*1, xdim3_advec_mom_kernel1_y_nonvector};
      const ptr_double vel1 = { vel1_p + n_x*1 + n_y * xdim4_advec_mom_kernel1_y_nonvector*1, xdim4_advec_mom_kernel1_y_nonvector};
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4

