//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_field_summary_kernel;
int xdim1_field_summary_kernel;
int xdim2_field_summary_kernel;
int xdim3_field_summary_kernel;
int xdim4_field_summary_kernel;
int xdim5_field_summary_kernel;

#define OPS_ACC0(x, y)                                                         \
  (n_x * 1 + n_y * xdim0_field_summary_kernel * 1 + x +                        \
   xdim0_field_summary_kernel * (y))
#define OPS_ACC1(x, y)                                                         \
  (n_x * 1 + n_y * xdim1_field_summary_kernel * 1 + x +                        \
   xdim1_field_summary_kernel * (y))
#define OPS_ACC2(x, y)                                                         \
  (n_x * 1 + n_y * xdim2_field_summary_kernel * 1 + x +                        \
   xdim2_field_summary_kernel * (y))
#define OPS_ACC3(x, y)                                                         \
  (n_x * 1 + n_y * xdim3_field_summary_kernel * 1 + x +                        \
   xdim3_field_summary_kernel * (y))
#define OPS_ACC4(x, y)                                                         \
  (n_x * 1 + n_y * xdim4_field_summary_kernel * 1 + x +                        \
   xdim4_field_summary_kernel * (y))
#define OPS_ACC5(x, y)                                                         \
  (n_x * 1 + n_y * xdim5_field_summary_kernel * 1 + x +                        \
   xdim5_field_summary_kernel * (y))

// user function

void field_summary_kernel_c_wrapper(
    const double *restrict volume, const double *restrict density0,
    const double *restrict energy0, const double *restrict pressure,
    const double *restrict xvel0, const double *restrict yvel0,
    double *restrict vol_g, double *restrict mass_g, double *restrict ie_g,
    double *restrict ke_g, double *restrict press_g, int x_size, int y_size) {
  double vol_0 = vol_g[0];
  double mass_0 = mass_g[0];
  double ie_0 = ie_g[0];
  double ke_0 = ke_g[0];
  double press_0 = press_g[0];
#pragma omp parallel for reduction(+ : vol_0) reduction(+ : mass_0) reduction( \
    + : ie_0) reduction(+ : ke_0) reduction(+ : press_0)
  for (int n_y = 0; n_y < y_size; n_y++) {
    for (int n_x = 0; n_x < x_size; n_x++) {
      double vol[1];
      vol[0] = ZERO_double;
      double mass[1];
      mass[0] = ZERO_double;
      double ie[1];
      ie[0] = ZERO_double;
      double ke[1];
      ke[0] = ZERO_double;
      double press[1];
      press[0] = ZERO_double;

      double vsqrd, cell_vol, cell_mass;

      vsqrd = 0.0;
      vsqrd = vsqrd +
              0.25 * (xvel0[OPS_ACC4(0, 0)] * xvel0[OPS_ACC4(0, 0)] +
                      yvel0[OPS_ACC5(0, 0)] * yvel0[OPS_ACC5(0, 0)]);
      vsqrd = vsqrd +
              0.25 * (xvel0[OPS_ACC4(1, 0)] * xvel0[OPS_ACC4(1, 0)] +
                      yvel0[OPS_ACC5(1, 0)] * yvel0[OPS_ACC5(1, 0)]);
      vsqrd = vsqrd +
              0.25 * (xvel0[OPS_ACC4(0, 1)] * xvel0[OPS_ACC4(0, 1)] +
                      yvel0[OPS_ACC5(0, 1)] * yvel0[OPS_ACC5(0, 1)]);
      vsqrd = vsqrd +
              0.25 * (xvel0[OPS_ACC4(1, 1)] * xvel0[OPS_ACC4(1, 1)] +
                      yvel0[OPS_ACC5(1, 1)] * yvel0[OPS_ACC5(1, 1)]);

      cell_vol = volume[OPS_ACC0(0, 0)];
      cell_mass = cell_vol * density0[OPS_ACC1(0, 0)];
      *vol = *vol + cell_vol;
      *mass = *mass + cell_mass;
      *ie = *ie + cell_mass * energy0[OPS_ACC2(0, 0)];
      *ke = *ke + cell_mass * 0.5 * vsqrd;
      *press = *press + cell_vol * pressure[OPS_ACC3(0, 0)];

      vol_0 += vol[0];
      mass_0 += mass[0];
      ie_0 += ie[0];
      ke_0 += ke[0];
      press_0 += press[0];
    }
  }
  vol_g[0] = vol_0;
  mass_g[0] = mass_0;
  ie_g[0] = ie_0;
  ke_g[0] = ke_0;
  press_g[0] = press_0;
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
