//
// auto-generated by ops.py
//

#ifdef OCL_FMA
#pragma OPENCL FP_CONTRACT ON
#else
#pragma OPENCL FP_CONTRACT OFF
#endif
#pragma OPENCL EXTENSION cl_khr_fp64:enable

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
#define OPS_ACC0(x) (x)
#define OPS_ACC1(x) (x)
#define OPS_ACC2(x) (x)
#define OPS_ACC3(x) (x)
#define OPS_ACC4(x) (x)


//user function
void calvar_kernel(const __global double * restrict rho_new,const __global double * restrict rhou_new,const __global double * restrict rhoE_new,
__global double * restrict workarray2,__global double * restrict workarray3,
  const double gam1)

 {
  double p, rhoi, u;
  rhoi = 1/rho_new[OPS_ACC0(0)];
  u = rhou_new[OPS_ACC1(0)] * rhoi;
  p = gam1 * (rhoE_new[OPS_ACC2(0)] - 0.5 * rho_new[OPS_ACC0(0)]* u * u);

  workarray2[OPS_ACC3(0)] = p + rhou_new[OPS_ACC1(0)] * u ;
  workarray3[OPS_ACC4(0)] = (p + rhoE_new[OPS_ACC2(0)]) * u ;
  }



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4



__kernel void ops_calvar_kernel(
__global const double* restrict arg0,
__global const double* restrict arg1,
__global const double* restrict arg2,
__global double* restrict arg3,
__global double* restrict arg4,
const double gam1,
const int base0,
const int base1,
const int base2,
const int base3,
const int base4,
const int size0 ){


  int idx_x = get_global_id(0);

  if (idx_x < size0) {
    calvar_kernel(&arg0[base0 + idx_x * 1*1],
                  &arg1[base1 + idx_x * 1*1],
                  &arg2[base2 + idx_x * 1*1],
                  &arg3[base3 + idx_x * 1*1],
                  &arg4[base4 + idx_x * 1*1],
                  gam1);
  }

}
