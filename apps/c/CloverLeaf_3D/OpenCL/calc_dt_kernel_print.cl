//
// auto-generated by ops.py
//

#ifdef OCL_FMA
#pragma OPENCL FP_CONTRACT ON
#else
#pragma OPENCL FP_CONTRACT OFF
#endif
#pragma OPENCL EXTENSION cl_khr_fp64:enable

#include "user_types.h"
#define OPS_3D
#define OPS_API 2
#define OPS_NO_GLOBALS
#include "ops_macros.h"
#include "ops_opencl_reduction.h"

#ifndef MIN
#define MIN(a,b) ((a<b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a>b) ? (a) : (b))
#endif
#ifndef SIGN
#define SIGN(a,b) ((b<0.0) ? (a*(-1)) : (a))
#endif
#define OPS_READ 0
#define OPS_WRITE 1
#define OPS_RW 2
#define OPS_INC 3
#define OPS_MIN 4
#define OPS_MAX 5
#define ZERO_double 0.0;
#define INFINITY_double INFINITY;
#define ZERO_float 0.0f;
#define INFINITY_float INFINITY;
#define ZERO_int 0;
#define INFINITY_int INFINITY;
#define ZERO_uint 0;
#define INFINITY_uint INFINITY;
#define ZERO_ll 0;
#define INFINITY_ll INFINITY;
#define ZERO_ull 0;
#define INFINITY_ull INFINITY;
#define ZERO_bool 0;

//user function

void calc_dt_kernel_print(const ptr_double xvel0,
  const ptr_double yvel0,
  const ptr_double zvel0,
  const ptr_double density0,
  const ptr_double energy0,
  const ptr_double pressure,
  const ptr_double soundspeed,
  double *output) {
  output[0] = OPS_ACCS(xvel0, 0,0,0);
  output[1] = OPS_ACCS(yvel0, 0,0,0);
  output[2] = OPS_ACCS(zvel0, 0,0,0);
  output[3] = OPS_ACCS(xvel0, 1,0,0);
  output[4] = OPS_ACCS(yvel0, 1,0,0);
  output[5] = OPS_ACCS(zvel0, 0,0,0);
  output[6] = OPS_ACCS(xvel0, 1,1,0);
  output[7] = OPS_ACCS(yvel0, 1,1,0);
  output[8] = OPS_ACCS(zvel0, 0,0,0);
  output[9] = OPS_ACCS(xvel0, 0,1,0);
  output[10] = OPS_ACCS(yvel0, 0,1,0);
  output[11] = OPS_ACCS(zvel0, 0,0,0);
  output[12] = OPS_ACCS(xvel0, 0,0,1);
  output[13] = OPS_ACCS(yvel0, 0,0,1);
  output[14] = OPS_ACCS(zvel0, 0,0,1);
  output[15] = OPS_ACCS(xvel0, 1,0,1);
  output[16] = OPS_ACCS(yvel0, 1,0,1);
  output[17] = OPS_ACCS(zvel0, 0,0,1);
  output[18] = OPS_ACCS(xvel0, 1,1,1);
  output[19] = OPS_ACCS(yvel0, 1,1,1);
  output[20] = OPS_ACCS(zvel0, 0,0,1);
  output[21] = OPS_ACCS(xvel0, 0,1,1);
  output[22] = OPS_ACCS(yvel0, 0,1,1);
  output[23] = OPS_ACCS(zvel0, 0,0,1);
  output[24] = OPS_ACCS(density0, 0,0,0);
  output[25] = OPS_ACCS(energy0, 0,0,0);
  output[26] = OPS_ACCS(pressure, 0,0,0);
  output[27] = OPS_ACCS(soundspeed, 0,0,0);

}


__kernel void ops_calc_dt_kernel_print(
__global const double* restrict arg0,
__global const double* restrict arg1,
__global const double* restrict arg2,
__global const double* restrict arg3,
__global const double* restrict arg4,
__global const double* restrict arg5,
__global const double* restrict arg6,
__global double* restrict arg7,
__local double* scratch7,
int r_bytes7,
const int base0,
const int base1,
const int base2,
const int base3,
const int base4,
const int base5,
const int base6,
const int size0,
const int size1,
const int size2 ){

  arg7 += r_bytes7;
  double arg7_l[28];
  for (int d=0; d<28; d++) arg7_l[d] = ZERO_double;

  int idx_y = get_global_id(1);
  int idx_z = get_global_id(2);
  int idx_x = get_global_id(0);

  if (idx_x < size0 && idx_y < size1 && idx_z < size2) {
    const ptr_double ptr0 = { &arg0[base0 + idx_x * 1*1 + idx_y * 1*1 * xdim0_calc_dt_kernel_print + idx_z * 1*1 * xdim0_calc_dt_kernel_print * ydim0_calc_dt_kernel_print], xdim0_calc_dt_kernel_print, ydim0_calc_dt_kernel_print};
    const ptr_double ptr1 = { &arg1[base1 + idx_x * 1*1 + idx_y * 1*1 * xdim1_calc_dt_kernel_print + idx_z * 1*1 * xdim1_calc_dt_kernel_print * ydim1_calc_dt_kernel_print], xdim1_calc_dt_kernel_print, ydim1_calc_dt_kernel_print};
    const ptr_double ptr2 = { &arg2[base2 + idx_x * 1*1 + idx_y * 1*1 * xdim2_calc_dt_kernel_print + idx_z * 1*1 * xdim2_calc_dt_kernel_print * ydim2_calc_dt_kernel_print], xdim2_calc_dt_kernel_print, ydim2_calc_dt_kernel_print};
    const ptr_double ptr3 = { &arg3[base3 + idx_x * 1*1 + idx_y * 1*1 * xdim3_calc_dt_kernel_print + idx_z * 1*1 * xdim3_calc_dt_kernel_print * ydim3_calc_dt_kernel_print], xdim3_calc_dt_kernel_print, ydim3_calc_dt_kernel_print};
    const ptr_double ptr4 = { &arg4[base4 + idx_x * 1*1 + idx_y * 1*1 * xdim4_calc_dt_kernel_print + idx_z * 1*1 * xdim4_calc_dt_kernel_print * ydim4_calc_dt_kernel_print], xdim4_calc_dt_kernel_print, ydim4_calc_dt_kernel_print};
    const ptr_double ptr5 = { &arg5[base5 + idx_x * 1*1 + idx_y * 1*1 * xdim5_calc_dt_kernel_print + idx_z * 1*1 * xdim5_calc_dt_kernel_print * ydim5_calc_dt_kernel_print], xdim5_calc_dt_kernel_print, ydim5_calc_dt_kernel_print};
    const ptr_double ptr6 = { &arg6[base6 + idx_x * 1*1 + idx_y * 1*1 * xdim6_calc_dt_kernel_print + idx_z * 1*1 * xdim6_calc_dt_kernel_print * ydim6_calc_dt_kernel_print], xdim6_calc_dt_kernel_print, ydim6_calc_dt_kernel_print};
    calc_dt_kernel_print(ptr0,
                   ptr1,
                   ptr2,
                   ptr3,
                   ptr4,
                   ptr5,
                   ptr6,
                   arg7_l);
  }
  int group_index = get_group_id(0) + get_group_id(1)*get_num_groups(0)+ get_group_id(2)*get_num_groups(0)*get_num_groups(1);
  for (int d=0; d<28; d++)
    reduce_double(arg7_l[d], scratch7, &arg7[group_index*28+d], OPS_INC);

}
