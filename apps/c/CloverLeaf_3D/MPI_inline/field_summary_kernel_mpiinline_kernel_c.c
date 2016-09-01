//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_field_summary_kernel;
int ydim0_field_summary_kernel;
int xdim1_field_summary_kernel;
int ydim1_field_summary_kernel;
int xdim2_field_summary_kernel;
int ydim2_field_summary_kernel;
int xdim3_field_summary_kernel;
int ydim3_field_summary_kernel;
int xdim4_field_summary_kernel;
int ydim4_field_summary_kernel;
int xdim5_field_summary_kernel;
int ydim5_field_summary_kernel;
int xdim6_field_summary_kernel;
int ydim6_field_summary_kernel;

#define OPS_ACC0(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim0_field_summary_kernel * 1 +                            \
   n_z * xdim0_field_summary_kernel * ydim0_field_summary_kernel * 1 + x +     \
   xdim0_field_summary_kernel * (y) +                                          \
   xdim0_field_summary_kernel * ydim0_field_summary_kernel * (z))
#define OPS_ACC1(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim1_field_summary_kernel * 1 +                            \
   n_z * xdim1_field_summary_kernel * ydim1_field_summary_kernel * 1 + x +     \
   xdim1_field_summary_kernel * (y) +                                          \
   xdim1_field_summary_kernel * ydim1_field_summary_kernel * (z))
#define OPS_ACC2(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim2_field_summary_kernel * 1 +                            \
   n_z * xdim2_field_summary_kernel * ydim2_field_summary_kernel * 1 + x +     \
   xdim2_field_summary_kernel * (y) +                                          \
   xdim2_field_summary_kernel * ydim2_field_summary_kernel * (z))
#define OPS_ACC3(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim3_field_summary_kernel * 1 +                            \
   n_z * xdim3_field_summary_kernel * ydim3_field_summary_kernel * 1 + x +     \
   xdim3_field_summary_kernel * (y) +                                          \
   xdim3_field_summary_kernel * ydim3_field_summary_kernel * (z))
#define OPS_ACC4(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim4_field_summary_kernel * 1 +                            \
   n_z * xdim4_field_summary_kernel * ydim4_field_summary_kernel * 1 + x +     \
   xdim4_field_summary_kernel * (y) +                                          \
   xdim4_field_summary_kernel * ydim4_field_summary_kernel * (z))
#define OPS_ACC5(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim5_field_summary_kernel * 1 +                            \
   n_z * xdim5_field_summary_kernel * ydim5_field_summary_kernel * 1 + x +     \
   xdim5_field_summary_kernel * (y) +                                          \
   xdim5_field_summary_kernel * ydim5_field_summary_kernel * (z))
#define OPS_ACC6(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim6_field_summary_kernel * 1 +                            \
   n_z * xdim6_field_summary_kernel * ydim6_field_summary_kernel * 1 + x +     \
   xdim6_field_summary_kernel * (y) +                                          \
   xdim6_field_summary_kernel * ydim6_field_summary_kernel * (z))

// user function

void field_summary_kernel_c_wrapper(
    const double *restrict volume, const double *restrict density0,
    const double *restrict energy0, const double *restrict pressure,
    const double *restrict xvel0, const double *restrict yvel0,
    const double *restrict zvel0, double *restrict vol_g,
    double *restrict mass_g, double *restrict ie_g, double *restrict ke_g,
    double *restrict press_g, int x_size, int y_size, int z_size) {
  double vol_v = *vol_g;
  double mass_v = *mass_g;
  double ie_v = *ie_g;
  double ke_v = *ke_g;
  double press_v = *press_g;
#pragma omp parallel for reduction(+ : vol_v) reduction(+ : mass_v) reduction( \
    + : ie_v) reduction(+ : ke_v) reduction(+ : press_v)
  for (int n_z = 0; n_z < z_size; n_z++) {
    for (int n_y = 0; n_y < y_size; n_y++) {
      for (int n_x = 0; n_x < x_size; n_x++) {
        double *restrict vol = &vol_v;
        double *restrict mass = &mass_v;
        double *restrict ie = &ie_v;
        double *restrict ke = &ke_v;
        double *restrict press = &press_v;

        double vsqrd, cell_vol, cell_mass;

        vsqrd = 0.0;
        vsqrd += 0.125 * (xvel0[OPS_ACC4(0, 0, 0)] * xvel0[OPS_ACC4(0, 0, 0)] +
                          yvel0[OPS_ACC5(0, 0, 0)] * yvel0[OPS_ACC5(0, 0, 0)] +
                          zvel0[OPS_ACC6(0, 0, 0)] * zvel0[OPS_ACC6(0, 0, 0)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(1, 0, 0)] * xvel0[OPS_ACC4(1, 0, 0)] +
                          yvel0[OPS_ACC5(1, 0, 0)] * yvel0[OPS_ACC5(1, 0, 0)] +
                          zvel0[OPS_ACC6(1, 0, 0)] * zvel0[OPS_ACC6(1, 0, 0)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(0, 1, 0)] * xvel0[OPS_ACC4(0, 1, 0)] +
                          yvel0[OPS_ACC5(0, 1, 0)] * yvel0[OPS_ACC5(0, 1, 0)] +
                          zvel0[OPS_ACC6(0, 1, 0)] * zvel0[OPS_ACC6(0, 1, 0)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(1, 1, 0)] * xvel0[OPS_ACC4(1, 1, 0)] +
                          yvel0[OPS_ACC5(1, 1, 0)] * yvel0[OPS_ACC5(1, 1, 0)] +
                          zvel0[OPS_ACC6(1, 1, 0)] * zvel0[OPS_ACC6(1, 1, 0)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(0, 0, 1)] * xvel0[OPS_ACC4(0, 0, 1)] +
                          yvel0[OPS_ACC5(0, 0, 1)] * yvel0[OPS_ACC5(0, 0, 1)] +
                          zvel0[OPS_ACC6(0, 0, 1)] * zvel0[OPS_ACC6(0, 0, 1)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(1, 0, 1)] * xvel0[OPS_ACC4(1, 0, 1)] +
                          yvel0[OPS_ACC5(1, 0, 1)] * yvel0[OPS_ACC5(1, 0, 1)] +
                          zvel0[OPS_ACC6(1, 0, 1)] * zvel0[OPS_ACC6(1, 0, 1)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(0, 1, 1)] * xvel0[OPS_ACC4(0, 1, 1)] +
                          yvel0[OPS_ACC5(0, 1, 1)] * yvel0[OPS_ACC5(0, 1, 1)] +
                          zvel0[OPS_ACC6(0, 1, 1)] * zvel0[OPS_ACC6(0, 1, 1)]);
        vsqrd += 0.125 * (xvel0[OPS_ACC4(1, 1, 1)] * xvel0[OPS_ACC4(1, 1, 1)] +
                          yvel0[OPS_ACC5(1, 1, 1)] * yvel0[OPS_ACC5(1, 1, 1)] +
                          zvel0[OPS_ACC6(1, 1, 1)] * zvel0[OPS_ACC6(1, 1, 1)]);

        cell_vol = volume[OPS_ACC0(0, 0, 0)];
        cell_mass = cell_vol * density0[OPS_ACC1(0, 0, 0)];
        *vol = *vol + cell_vol;
        *mass = *mass + cell_mass;
        *ie = *ie + cell_mass * energy0[OPS_ACC2(0, 0, 0)];
        *ke = *ke + cell_mass * 0.5 * vsqrd;
        *press = *press + cell_vol * pressure[OPS_ACC3(0, 0, 0)];
      }
    }
  }
  *vol_g = vol_v;
  *mass_g = mass_v;
  *ie_g = ie_v;
  *ke_g = ke_v;
  *press_g = press_v;
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6