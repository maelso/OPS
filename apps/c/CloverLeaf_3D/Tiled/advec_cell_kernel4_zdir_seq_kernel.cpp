//
// auto-generated by ops.py
//
#define OPS_ACC0(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim0_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim0_advec_cell_kernel4_zdir * ydim0_advec_cell_kernel4_zdir * 1 +   \
   x + xdim0_advec_cell_kernel4_zdir * (y) +                                   \
   xdim0_advec_cell_kernel4_zdir * ydim0_advec_cell_kernel4_zdir * (z))
#define OPS_ACC1(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim1_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim1_advec_cell_kernel4_zdir * ydim1_advec_cell_kernel4_zdir * 1 +   \
   x + xdim1_advec_cell_kernel4_zdir * (y) +                                   \
   xdim1_advec_cell_kernel4_zdir * ydim1_advec_cell_kernel4_zdir * (z))
#define OPS_ACC2(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim2_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim2_advec_cell_kernel4_zdir * ydim2_advec_cell_kernel4_zdir * 1 +   \
   x + xdim2_advec_cell_kernel4_zdir * (y) +                                   \
   xdim2_advec_cell_kernel4_zdir * ydim2_advec_cell_kernel4_zdir * (z))
#define OPS_ACC3(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim3_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim3_advec_cell_kernel4_zdir * ydim3_advec_cell_kernel4_zdir * 1 +   \
   x + xdim3_advec_cell_kernel4_zdir * (y) +                                   \
   xdim3_advec_cell_kernel4_zdir * ydim3_advec_cell_kernel4_zdir * (z))
#define OPS_ACC4(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim4_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim4_advec_cell_kernel4_zdir * ydim4_advec_cell_kernel4_zdir * 1 +   \
   x + xdim4_advec_cell_kernel4_zdir * (y) +                                   \
   xdim4_advec_cell_kernel4_zdir * ydim4_advec_cell_kernel4_zdir * (z))
#define OPS_ACC5(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim5_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim5_advec_cell_kernel4_zdir * ydim5_advec_cell_kernel4_zdir * 1 +   \
   x + xdim5_advec_cell_kernel4_zdir * (y) +                                   \
   xdim5_advec_cell_kernel4_zdir * ydim5_advec_cell_kernel4_zdir * (z))
#define OPS_ACC6(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim6_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim6_advec_cell_kernel4_zdir * ydim6_advec_cell_kernel4_zdir * 1 +   \
   x + xdim6_advec_cell_kernel4_zdir * (y) +                                   \
   xdim6_advec_cell_kernel4_zdir * ydim6_advec_cell_kernel4_zdir * (z))
#define OPS_ACC7(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim7_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim7_advec_cell_kernel4_zdir * ydim7_advec_cell_kernel4_zdir * 1 +   \
   x + xdim7_advec_cell_kernel4_zdir * (y) +                                   \
   xdim7_advec_cell_kernel4_zdir * ydim7_advec_cell_kernel4_zdir * (z))
#define OPS_ACC8(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim8_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim8_advec_cell_kernel4_zdir * ydim8_advec_cell_kernel4_zdir * 1 +   \
   x + xdim8_advec_cell_kernel4_zdir * (y) +                                   \
   xdim8_advec_cell_kernel4_zdir * ydim8_advec_cell_kernel4_zdir * (z))
#define OPS_ACC9(x, y, z)                                                      \
  (n_x * 1 + n_y * xdim9_advec_cell_kernel4_zdir * 1 +                         \
   n_z * xdim9_advec_cell_kernel4_zdir * ydim9_advec_cell_kernel4_zdir * 1 +   \
   x + xdim9_advec_cell_kernel4_zdir * (y) +                                   \
   xdim9_advec_cell_kernel4_zdir * ydim9_advec_cell_kernel4_zdir * (z))
#define OPS_ACC10(x, y, z)                                                     \
  (n_x * 1 + n_y * xdim10_advec_cell_kernel4_zdir * 1 +                        \
   n_z * xdim10_advec_cell_kernel4_zdir * ydim10_advec_cell_kernel4_zdir * 1 + \
   x + xdim10_advec_cell_kernel4_zdir * (y) +                                  \
   xdim10_advec_cell_kernel4_zdir * ydim10_advec_cell_kernel4_zdir * (z))

// user function

// host stub function
void ops_par_loop_advec_cell_kernel4_zdir_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];
  ops_arg arg3 = desc->args[3];
  ops_arg arg4 = desc->args[4];
  ops_arg arg5 = desc->args[5];
  ops_arg arg6 = desc->args[6];
  ops_arg arg7 = desc->args[7];
  ops_arg arg8 = desc->args[8];
  ops_arg arg9 = desc->args[9];
  ops_arg arg10 = desc->args[10];

  // Timing
  double t1, t2, c1, c2;

  ops_arg args[11] = {arg0, arg1, arg2, arg3, arg4, arg5,
                      arg6, arg7, arg8, arg9, arg10};

#ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args, 11, range, 18))
    return;
#endif

  if (OPS_diags > 1) {
    OPS_kernels[18].count++;
    ops_timers_core(&c2, &t2);
  }

  // compute locally allocated range for the sub-block
  int start[3];
  int end[3];

  for (int n = 0; n < 3; n++) {
    start[n] = range[2 * n];
    end[n] = range[2 * n + 1];
  }

#ifdef OPS_DEBUG
  ops_register_args(args, "advec_cell_kernel4_zdir");
#endif

  // set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  double *__restrict__ density1 = (double *)(args[0].data + base0);

  int base1 = args[1].dat->base_offset;
  double *__restrict__ energy1 = (double *)(args[1].data + base1);

  int base2 = args[2].dat->base_offset;
  const double *__restrict__ mass_flux_z = (double *)(args[2].data + base2);

  int base3 = args[3].dat->base_offset;
  const double *__restrict__ vol_flux_z = (double *)(args[3].data + base3);

  int base4 = args[4].dat->base_offset;
  const double *__restrict__ pre_vol = (double *)(args[4].data + base4);

  int base5 = args[5].dat->base_offset;
  const double *__restrict__ post_vol = (double *)(args[5].data + base5);

  int base6 = args[6].dat->base_offset;
  double *__restrict__ pre_mass = (double *)(args[6].data + base6);

  int base7 = args[7].dat->base_offset;
  double *__restrict__ post_mass = (double *)(args[7].data + base7);

  int base8 = args[8].dat->base_offset;
  double *__restrict__ advec_vol = (double *)(args[8].data + base8);

  int base9 = args[9].dat->base_offset;
  double *__restrict__ post_ener = (double *)(args[9].data + base9);

  int base10 = args[10].dat->base_offset;
  const double *__restrict__ ener_flux = (double *)(args[10].data + base10);

  // initialize global variable with the dimension of dats
  int xdim0_advec_cell_kernel4_zdir = args[0].dat->size[0];
  int ydim0_advec_cell_kernel4_zdir = args[0].dat->size[1];
  int xdim1_advec_cell_kernel4_zdir = args[1].dat->size[0];
  int ydim1_advec_cell_kernel4_zdir = args[1].dat->size[1];
  int xdim2_advec_cell_kernel4_zdir = args[2].dat->size[0];
  int ydim2_advec_cell_kernel4_zdir = args[2].dat->size[1];
  int xdim3_advec_cell_kernel4_zdir = args[3].dat->size[0];
  int ydim3_advec_cell_kernel4_zdir = args[3].dat->size[1];
  int xdim4_advec_cell_kernel4_zdir = args[4].dat->size[0];
  int ydim4_advec_cell_kernel4_zdir = args[4].dat->size[1];
  int xdim5_advec_cell_kernel4_zdir = args[5].dat->size[0];
  int ydim5_advec_cell_kernel4_zdir = args[5].dat->size[1];
  int xdim6_advec_cell_kernel4_zdir = args[6].dat->size[0];
  int ydim6_advec_cell_kernel4_zdir = args[6].dat->size[1];
  int xdim7_advec_cell_kernel4_zdir = args[7].dat->size[0];
  int ydim7_advec_cell_kernel4_zdir = args[7].dat->size[1];
  int xdim8_advec_cell_kernel4_zdir = args[8].dat->size[0];
  int ydim8_advec_cell_kernel4_zdir = args[8].dat->size[1];
  int xdim9_advec_cell_kernel4_zdir = args[9].dat->size[0];
  int ydim9_advec_cell_kernel4_zdir = args[9].dat->size[1];
  int xdim10_advec_cell_kernel4_zdir = args[10].dat->size[0];
  int ydim10_advec_cell_kernel4_zdir = args[10].dat->size[1];

  if (OPS_diags > 1) {
    ops_timers_core(&c1, &t1);
    OPS_kernels[18].mpi_time += t1 - t2;
  }

#pragma omp parallel for collapse(2)
  for (int n_z = start[2]; n_z < end[2]; n_z++) {
    for (int n_y = start[1]; n_y < end[1]; n_y++) {
#ifdef intel
#pragma loop_count(10000)
#pragma omp simd aligned(density1, energy1, mass_flux_z, vol_flux_z, pre_vol,  \
                         post_vol, pre_mass, post_mass, advec_vol, post_ener,  \
                         ener_flux)
#else
#pragma simd
#endif
      for (int n_x = start[0]; n_x < end[0]; n_x++) {

        pre_mass[OPS_ACC6(0, 0, 0)] =
            density1[OPS_ACC0(0, 0, 0)] * pre_vol[OPS_ACC4(0, 0, 0)];
        post_mass[OPS_ACC7(0, 0, 0)] = pre_mass[OPS_ACC6(0, 0, 0)] +
                                       mass_flux_z[OPS_ACC2(0, 0, 0)] -
                                       mass_flux_z[OPS_ACC2(0, 0, 1)];
        post_ener[OPS_ACC9(0, 0, 0)] =
            (energy1[OPS_ACC1(0, 0, 0)] * pre_mass[OPS_ACC6(0, 0, 0)] +
             ener_flux[OPS_ACC10(0, 0, 0)] - ener_flux[OPS_ACC10(0, 0, 1)]) /
            post_mass[OPS_ACC7(0, 0, 0)];
        advec_vol[OPS_ACC8(0, 0, 0)] = pre_vol[OPS_ACC4(0, 0, 0)] +
                                       vol_flux_z[OPS_ACC3(0, 0, 0)] -
                                       vol_flux_z[OPS_ACC3(0, 0, 1)];
        density1[OPS_ACC0(0, 0, 0)] =
            post_mass[OPS_ACC7(0, 0, 0)] / advec_vol[OPS_ACC8(0, 0, 0)];
        energy1[OPS_ACC1(0, 0, 0)] = post_ener[OPS_ACC9(0, 0, 0)];
      }
    }
  }
  if (OPS_diags > 1) {
    ops_timers_core(&c2, &t2);
    OPS_kernels[18].time += t2 - t1;
  }

  if (OPS_diags > 1) {
    // Update kernel record
    ops_timers_core(&c1, &t1);
    OPS_kernels[18].mpi_time += t1 - t2;
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg1);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg2);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg3);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg4);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg5);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg6);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg7);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg8);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg9);
    OPS_kernels[18].transfer += ops_compute_transfer(dim, start, end, &arg10);
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

void ops_par_loop_advec_cell_kernel4_zdir(
    char const *name, ops_block block, int dim, int *range, ops_arg arg0,
    ops_arg arg1, ops_arg arg2, ops_arg arg3, ops_arg arg4, ops_arg arg5,
    ops_arg arg6, ops_arg arg7, ops_arg arg8, ops_arg arg9, ops_arg arg10) {
  ops_kernel_descriptor *desc =
      (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->index = 18;
  desc->hash = 5381;
  desc->hash = ((desc->hash << 5) + desc->hash) + 18;
  for (int i = 0; i < 6; i++) {
    desc->range[i] = range[i];
    desc->orig_range[i] = range[i];
    desc->hash = ((desc->hash << 5) + desc->hash) + range[i];
  }
  desc->nargs = 11;
  desc->args = (ops_arg *)malloc(11 * sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg0.dat->index;
  desc->args[1] = arg1;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg1.dat->index;
  desc->args[2] = arg2;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg2.dat->index;
  desc->args[3] = arg3;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg3.dat->index;
  desc->args[4] = arg4;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg4.dat->index;
  desc->args[5] = arg5;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg5.dat->index;
  desc->args[6] = arg6;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg6.dat->index;
  desc->args[7] = arg7;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg7.dat->index;
  desc->args[8] = arg8;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg8.dat->index;
  desc->args[9] = arg9;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg9.dat->index;
  desc->args[10] = arg10;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg10.dat->index;
  desc->function = ops_par_loop_advec_cell_kernel4_zdir_execute;
  if (OPS_diags > 1) {
    ops_timing_realloc(18, "advec_cell_kernel4_zdir");
  }
  ops_enqueue_kernel(desc);
}
