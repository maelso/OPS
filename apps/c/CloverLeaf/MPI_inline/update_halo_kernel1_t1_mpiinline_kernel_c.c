//
// auto-generated by ops.py
//

int xdim0_update_halo_kernel1_t1;
int xdim1_update_halo_kernel1_t1;
int xdim2_update_halo_kernel1_t1;
int xdim3_update_halo_kernel1_t1;
int xdim4_update_halo_kernel1_t1;
int xdim5_update_halo_kernel1_t1;
int xdim6_update_halo_kernel1_t1;


//user function



void update_halo_kernel1_t1_c_wrapper(
  double * restrict density0_p,
  double * restrict density1_p,
  double * restrict energy0_p,
  double * restrict energy1_p,
  double * restrict pressure_p,
  double * restrict viscosity_p,
  double * restrict soundspeed_p,
  const int * restrict fields,
  int x_size, int y_size) {
  #pragma omp parallel for
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      ptr_double density0 = { density0_p + n_x*1 + n_y * xdim0_update_halo_kernel1_t1*1, xdim0_update_halo_kernel1_t1};
      ptr_double density1 = { density1_p + n_x*1 + n_y * xdim1_update_halo_kernel1_t1*1, xdim1_update_halo_kernel1_t1};
      ptr_double energy0 = { energy0_p + n_x*1 + n_y * xdim2_update_halo_kernel1_t1*1, xdim2_update_halo_kernel1_t1};
      ptr_double energy1 = { energy1_p + n_x*1 + n_y * xdim3_update_halo_kernel1_t1*1, xdim3_update_halo_kernel1_t1};
      ptr_double pressure = { pressure_p + n_x*1 + n_y * xdim4_update_halo_kernel1_t1*1, xdim4_update_halo_kernel1_t1};
      ptr_double viscosity = { viscosity_p + n_x*1 + n_y * xdim5_update_halo_kernel1_t1*1, xdim5_update_halo_kernel1_t1};
      ptr_double soundspeed = { soundspeed_p + n_x*1 + n_y * xdim6_update_halo_kernel1_t1*1, xdim6_update_halo_kernel1_t1};
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

