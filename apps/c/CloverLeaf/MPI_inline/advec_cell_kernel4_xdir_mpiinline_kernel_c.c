//
// auto-generated by ops.py
//

int xdim0_advec_cell_kernel4_xdir;
int xdim1_advec_cell_kernel4_xdir;
int xdim2_advec_cell_kernel4_xdir;
int xdim3_advec_cell_kernel4_xdir;
int xdim4_advec_cell_kernel4_xdir;
int xdim5_advec_cell_kernel4_xdir;
int xdim6_advec_cell_kernel4_xdir;
int xdim7_advec_cell_kernel4_xdir;
int xdim8_advec_cell_kernel4_xdir;
int xdim9_advec_cell_kernel4_xdir;
int xdim10_advec_cell_kernel4_xdir;


//user function



void advec_cell_kernel4_xdir_c_wrapper(
  double * restrict density1_p,
  double * restrict energy1_p,
  double * restrict mass_flux_x_p,
  double * restrict vol_flux_x_p,
  double * restrict pre_vol_p,
  double * restrict post_vol_p,
  double * restrict pre_mass_p,
  double * restrict post_mass_p,
  double * restrict advec_vol_p,
  double * restrict post_ener_p,
  double * restrict ener_flux_p,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      ptr_double density1 = { density1_p + n_x*1 + n_y * xdim0_advec_cell_kernel4_xdir*1, xdim0_advec_cell_kernel4_xdir};
      ptr_double energy1 = { energy1_p + n_x*1 + n_y * xdim1_advec_cell_kernel4_xdir*1, xdim1_advec_cell_kernel4_xdir};
      const ptr_double mass_flux_x = { mass_flux_x_p + n_x*1 + n_y * xdim2_advec_cell_kernel4_xdir*1, xdim2_advec_cell_kernel4_xdir};
      const ptr_double vol_flux_x = { vol_flux_x_p + n_x*1 + n_y * xdim3_advec_cell_kernel4_xdir*1, xdim3_advec_cell_kernel4_xdir};
      const ptr_double pre_vol = { pre_vol_p + n_x*1 + n_y * xdim4_advec_cell_kernel4_xdir*1, xdim4_advec_cell_kernel4_xdir};
      const ptr_double post_vol = { post_vol_p + n_x*1 + n_y * xdim5_advec_cell_kernel4_xdir*1, xdim5_advec_cell_kernel4_xdir};
      ptr_double pre_mass = { pre_mass_p + n_x*1 + n_y * xdim6_advec_cell_kernel4_xdir*1, xdim6_advec_cell_kernel4_xdir};
      ptr_double post_mass = { post_mass_p + n_x*1 + n_y * xdim7_advec_cell_kernel4_xdir*1, xdim7_advec_cell_kernel4_xdir};
      ptr_double advec_vol = { advec_vol_p + n_x*1 + n_y * xdim8_advec_cell_kernel4_xdir*1, xdim8_advec_cell_kernel4_xdir};
      ptr_double post_ener = { post_ener_p + n_x*1 + n_y * xdim9_advec_cell_kernel4_xdir*1, xdim9_advec_cell_kernel4_xdir};
      const ptr_double ener_flux = { ener_flux_p + n_x*1 + n_y * xdim10_advec_cell_kernel4_xdir*1, xdim10_advec_cell_kernel4_xdir};
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6
#undef OPS_ACC7
#undef OPS_ACC8
#undef OPS_ACC9
#undef OPS_ACC10

