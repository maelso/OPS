//
// auto-generated by ops.py
//

#define OPS_GPU

int xdim0_flux_calc_kernelx;
int ydim0_flux_calc_kernelx;
int xdim1_flux_calc_kernelx;
int ydim1_flux_calc_kernelx;
int xdim2_flux_calc_kernelx;
int ydim2_flux_calc_kernelx;
int xdim3_flux_calc_kernelx;
int ydim3_flux_calc_kernelx;

#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3

#define OPS_ACC0(x, y, z)                                                      \
  (x + xdim0_flux_calc_kernelx * (y) +                                         \
   xdim0_flux_calc_kernelx * ydim0_flux_calc_kernelx * (z))
#define OPS_ACC1(x, y, z)                                                      \
  (x + xdim1_flux_calc_kernelx * (y) +                                         \
   xdim1_flux_calc_kernelx * ydim1_flux_calc_kernelx * (z))
#define OPS_ACC2(x, y, z)                                                      \
  (x + xdim2_flux_calc_kernelx * (y) +                                         \
   xdim2_flux_calc_kernelx * ydim2_flux_calc_kernelx * (z))
#define OPS_ACC3(x, y, z)                                                      \
  (x + xdim3_flux_calc_kernelx * (y) +                                         \
   xdim3_flux_calc_kernelx * ydim3_flux_calc_kernelx * (z))

// user function
inline void flux_calc_kernelx(double *vol_flux_x, const double *xarea,
                              const double *xvel0, const double *xvel1) {

  vol_flux_x[OPS_ACC0(0, 0, 0)] =
      0.125 * dt * (xarea[OPS_ACC1(0, 0, 0)]) *
      (xvel0[OPS_ACC2(0, 0, 0)] + xvel0[OPS_ACC2(0, 1, 0)] +
       xvel0[OPS_ACC2(0, 0, 1)] + xvel0[OPS_ACC2(0, 1, 1)] +
       xvel1[OPS_ACC3(0, 0, 0)] + xvel1[OPS_ACC3(0, 1, 0)] +
       xvel1[OPS_ACC3(0, 0, 1)] + xvel1[OPS_ACC3(0, 1, 1)]);
}

#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3

void flux_calc_kernelx_c_wrapper(double *p_a0, double *p_a1, double *p_a2,
                                 double *p_a3, int x_size, int y_size,
                                 int z_size) {
#ifdef OPS_GPU
#pragma acc parallel deviceptr(p_a0, p_a1, p_a2, p_a3)
#pragma acc loop
#endif
  for (int n_z = 0; n_z < z_size; n_z++) {
#ifdef OPS_GPU
#pragma acc loop
#endif
    for (int n_y = 0; n_y < y_size; n_y++) {
#ifdef OPS_GPU
#pragma acc loop
#endif
      for (int n_x = 0; n_x < x_size; n_x++) {
        flux_calc_kernelx(
            p_a0 + n_x * 1 * 1 + n_y * xdim0_flux_calc_kernelx * 1 * 1 +
                n_z * xdim0_flux_calc_kernelx * ydim0_flux_calc_kernelx * 1 * 1,
            p_a1 + n_x * 1 * 1 + n_y * xdim1_flux_calc_kernelx * 1 * 1 +
                n_z * xdim1_flux_calc_kernelx * ydim1_flux_calc_kernelx * 1 * 1,
            p_a2 + n_x * 1 * 1 + n_y * xdim2_flux_calc_kernelx * 1 * 1 +
                n_z * xdim2_flux_calc_kernelx * ydim2_flux_calc_kernelx * 1 * 1,
            p_a3 + n_x * 1 * 1 + n_y * xdim3_flux_calc_kernelx * 1 * 1 +
                n_z * xdim3_flux_calc_kernelx * ydim3_flux_calc_kernelx * 1 *
                    1);
      }
    }
  }
}
