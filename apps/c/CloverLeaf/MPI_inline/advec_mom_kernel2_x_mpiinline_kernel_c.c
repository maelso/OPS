//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_advec_mom_kernel2_x;
int xdim1_advec_mom_kernel2_x;
int xdim2_advec_mom_kernel2_x;
int xdim3_advec_mom_kernel2_x;

#define OPS_ACC0(x, y)                                                         \
  (n_x * 1 + n_y * xdim0_advec_mom_kernel2_x * 1 + x +                         \
   xdim0_advec_mom_kernel2_x * (y))
#define OPS_ACC1(x, y)                                                         \
  (n_x * 1 + n_y * xdim1_advec_mom_kernel2_x * 1 + x +                         \
   xdim1_advec_mom_kernel2_x * (y))
#define OPS_ACC2(x, y)                                                         \
  (n_x * 1 + n_y * xdim2_advec_mom_kernel2_x * 1 + x +                         \
   xdim2_advec_mom_kernel2_x * (y))
#define OPS_ACC3(x, y)                                                         \
  (n_x * 1 + n_y * xdim3_advec_mom_kernel2_x * 1 + x +                         \
   xdim3_advec_mom_kernel2_x * (y))

// user function

void advec_mom_kernel2_x_c_wrapper(double *restrict vel1,
                                   const double *restrict node_mass_post,
                                   const double *restrict node_mass_pre,
                                   const double *restrict mom_flux, int x_size,
                                   int y_size) {
#pragma omp parallel for
  for (int n_y = 0; n_y < y_size; n_y++) {
    for (int n_x = 0; n_x < x_size; n_x++) {

      vel1[OPS_ACC0(0, 0)] =
          (vel1[OPS_ACC0(0, 0)] * node_mass_pre[OPS_ACC2(0, 0)] +
           mom_flux[OPS_ACC3(-1, 0)] - mom_flux[OPS_ACC3(0, 0)]) /
          node_mass_post[OPS_ACC1(0, 0)];
    }
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
