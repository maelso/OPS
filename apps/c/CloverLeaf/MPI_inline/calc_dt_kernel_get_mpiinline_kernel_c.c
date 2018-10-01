//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

int xdim0_calc_dt_kernel_get;
int xdim1_calc_dt_kernel_get;

#define OPS_ACC0(x, y)                                                         \
  (n_x * 1 + n_y * xdim0_calc_dt_kernel_get * 0 + x +                          \
   xdim0_calc_dt_kernel_get * (y))
#define OPS_ACC1(x, y)                                                         \
  (n_x * 0 + n_y * xdim1_calc_dt_kernel_get * 1 + x +                          \
   xdim1_calc_dt_kernel_get * (y))

// user function

void calc_dt_kernel_get_c_wrapper(const double *restrict cellx,
                                  const double *restrict celly,
                                  double *restrict xl_pos_g,
                                  double *restrict yl_pos_g, int x_size,
                                  int y_size) {
  double xl_pos_0 = xl_pos_g[0];
  double yl_pos_0 = yl_pos_g[0];
#pragma omp parallel for reduction(+ : xl_pos_0) reduction(+ : yl_pos_0)
  for (int n_y = 0; n_y < y_size; n_y++) {
    for (int n_x = 0; n_x < x_size; n_x++) {
      double xl_pos[1];
      xl_pos[0] = ZERO_double;
      double yl_pos[1];
      yl_pos[0] = ZERO_double;

      *xl_pos = cellx[OPS_ACC0(0, 0)];
      *yl_pos = celly[OPS_ACC1(0, 0)];

      xl_pos_0 += xl_pos[0];
      yl_pos_0 += yl_pos[0];
    }
  }
  xl_pos_g[0] = xl_pos_0;
  yl_pos_g[0] = yl_pos_0;
}
#undef OPS_ACC0
#undef OPS_ACC1
