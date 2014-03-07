//
// auto-generated by op2.py on 2014-03-07 14:18
//

//header
#include "ops_lib_cpp.h"
#include "ops_cuda_rt_support.h"
#include "ops_cuda_reduction.h"

// global constants
__constant__ double g_small;
__constant__ double g_big;
__constant__ double dtc_safe;
__constant__ double dtu_safe;
__constant__ double dtv_safe;
__constant__ double dtdiv_safe;
__constant__ int x_max;
__constant__ int y_max;
__constant__ double dt;

#define FIELD_DENSITY0 0
#define FIELD_DENSITY1 1
#define FIELD_ENERGY0 2
#define FIELD_ENERGY1 3
#define FIELD_PRESSURE 4
#define FIELD_VISCOSITY 5
#define FIELD_SOUNDSPEED 6
#define FIELD_XVEL0 7
#define FIELD_XVEL1 8
#define FIELD_YVEL0 9
#define FIELD_YVEL1 10
#define FIELD_VOL_FLUX_X 11
#define FIELD_VOL_FLUX_Y 12
#define FIELD_MASS_FLUX_X 13
#define FIELD_MASS_FLUX_Y 14
#define NUM_FIELDS 15

void ops_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!strcmp(name,"g_small")) {
    cutilSafeCall(cudaMemcpyToSymbol(g_small, dat, dim*size));
  }
  else
  if (!strcmp(name,"g_big")) {
    cutilSafeCall(cudaMemcpyToSymbol(g_big, dat, dim*size));
  }
  else
  if (!strcmp(name,"dtc_safe")) {
    cutilSafeCall(cudaMemcpyToSymbol(dtc_safe, dat, dim*size));
  }
  else
  if (!strcmp(name,"dtu_safe")) {
    cutilSafeCall(cudaMemcpyToSymbol(dtu_safe, dat, dim*size));
  }
  else
  if (!strcmp(name,"dtv_safe")) {
    cutilSafeCall(cudaMemcpyToSymbol(dtv_safe, dat, dim*size));
  }
  else
  if (!strcmp(name,"dtdiv_safe")) {
    cutilSafeCall(cudaMemcpyToSymbol(dtdiv_safe, dat, dim*size));
  }
  else
  if (!strcmp(name,"x_max")) {
    cutilSafeCall(cudaMemcpyToSymbol(x_max, dat, dim*size));
  }
  else
  if (!strcmp(name,"y_max")) {
    cutilSafeCall(cudaMemcpyToSymbol(y_max, dat, dim*size));
  }
  else
  if (!strcmp(name,"dt")) {
    cutilSafeCall(cudaMemcpyToSymbol(dt, dat, dim*size));
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}


//user kernel files
#include "revert_kernel_cuda_kernel.cu"
#include "reset_field_kernel1_cuda_kernel.cu"
#include "reset_field_kernel2_cuda_kernel.cu"
#include "ideal_gas_kernel_cuda_kernel.cu"
#include "PdV_kernel_predict_cuda_kernel.cu"
#include "PdV_kernel_nopredict_cuda_kernel.cu"
#include "accelerate_kernel_stepbymass_cuda_kernel.cu"
#include "accelerate_kernelx1_cuda_kernel.cu"
#include "accelerate_kernely1_cuda_kernel.cu"
#include "accelerate_kernelx2_cuda_kernel.cu"
#include "accelerate_kernely2_cuda_kernel.cu"
#include "advec_cell_kernel1_xdir_cuda_kernel.cu"
#include "advec_cell_kernel2_xdir_cuda_kernel.cu"
#include "advec_cell_kernel3_xdir_cuda_kernel.cu"
#include "advec_cell_kernel4_xdir_cuda_kernel.cu"
#include "advec_cell_kernel1_ydir_cuda_kernel.cu"
#include "advec_cell_kernel2_ydir_cuda_kernel.cu"
#include "advec_cell_kernel3_ydir_cuda_kernel.cu"
#include "advec_cell_kernel4_ydir_cuda_kernel.cu"
#include "advec_mom_kernel_x1_cuda_kernel.cu"
#include "advec_mom_kernel_y1_cuda_kernel.cu"
#include "advec_mom_kernel_x2_cuda_kernel.cu"
#include "advec_mom_kernel_y2_cuda_kernel.cu"
#include "advec_mom_kernel_mass_flux_x_cuda_kernel.cu"
#include "advec_mom_kernel_post_advec_x_cuda_kernel.cu"
#include "advec_mom_kernel_pre_advec_x_cuda_kernel.cu"
#include "advec_mom_kernel1_x_nonvector_cuda_kernel.cu"
#include "advec_mom_kernel2_x_cuda_kernel.cu"
#include "advec_mom_kernel_mass_flux_y_cuda_kernel.cu"
#include "advec_mom_kernel_post_advec_y_cuda_kernel.cu"
#include "advec_mom_kernel_pre_advec_y_cuda_kernel.cu"
#include "advec_mom_kernel1_y_nonvector_cuda_kernel.cu"
#include "advec_mom_kernel2_y_cuda_kernel.cu"
#include "calc_dt_kernel_cuda_kernel.cu"
#include "calc_dt_kernel_min_cuda_kernel.cu"
#include "calc_dt_kernel_get_cuda_kernel.cu"
#include "calc_dt_kernel_print_cuda_kernel.cu"
#include "field_summary_kernel_cuda_kernel.cu"
#include "flux_calc_kernelx_cuda_kernel.cu"
#include "flux_calc_kernely_cuda_kernel.cu"
#include "viscosity_kernel_cuda_kernel.cu"
#include "update_halo_kernel1_b2_cuda_kernel.cu"
#include "update_halo_kernel1_b1_cuda_kernel.cu"
#include "update_halo_kernel1_t2_cuda_kernel.cu"
#include "update_halo_kernel1_t1_cuda_kernel.cu"
#include "update_halo_kernel1_l2_cuda_kernel.cu"
#include "update_halo_kernel1_l1_cuda_kernel.cu"
#include "update_halo_kernel1_r2_cuda_kernel.cu"
#include "update_halo_kernel1_r1_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_plus_4_a_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_plus_2_a_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_plus_4_b_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_plus_2_b_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_minus_4_a_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_minus_2_a_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_minus_4_b_cuda_kernel.cu"
#include "update_halo_kernel2_xvel_minus_2_b_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_minus_4_a_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_minus_2_a_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_minus_4_b_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_minus_2_b_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_plus_4_a_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_plus_2_a_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_plus_4_b_cuda_kernel.cu"
#include "update_halo_kernel2_yvel_plus_2_b_cuda_kernel.cu"
#include "update_halo_kernel3_plus_4_a_cuda_kernel.cu"
#include "update_halo_kernel3_plus_2_a_cuda_kernel.cu"
#include "update_halo_kernel3_plus_4_b_cuda_kernel.cu"
#include "update_halo_kernel3_plus_2_b_cuda_kernel.cu"
#include "update_halo_kernel3_minus_4_a_cuda_kernel.cu"
#include "update_halo_kernel3_minus_2_a_cuda_kernel.cu"
#include "update_halo_kernel3_minus_4_b_cuda_kernel.cu"
#include "update_halo_kernel3_minus_2_b_cuda_kernel.cu"
#include "update_halo_kernel4_minus_4_a_cuda_kernel.cu"
#include "update_halo_kernel4_minus_2_a_cuda_kernel.cu"
#include "update_halo_kernel4_minus_4_b_cuda_kernel.cu"
#include "update_halo_kernel4_minus_2_b_cuda_kernel.cu"
#include "update_halo_kernel4_plus_4_a_cuda_kernel.cu"
#include "update_halo_kernel4_plus_2_a_cuda_kernel.cu"
#include "update_halo_kernel4_plus_4_b_cuda_kernel.cu"
#include "update_halo_kernel4_plus_2_b_cuda_kernel.cu"
