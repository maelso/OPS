//
// auto-generated by ops.py
//
#include "./MPI_inline/clover_leaf_common.h"

extern int xdim0_PdV_kernel_nopredict;
int xdim0_PdV_kernel_nopredict_h = -1;
extern int xdim1_PdV_kernel_nopredict;
int xdim1_PdV_kernel_nopredict_h = -1;
extern int xdim2_PdV_kernel_nopredict;
int xdim2_PdV_kernel_nopredict_h = -1;
extern int xdim3_PdV_kernel_nopredict;
int xdim3_PdV_kernel_nopredict_h = -1;
extern int xdim4_PdV_kernel_nopredict;
int xdim4_PdV_kernel_nopredict_h = -1;
extern int xdim5_PdV_kernel_nopredict;
int xdim5_PdV_kernel_nopredict_h = -1;
extern int xdim6_PdV_kernel_nopredict;
int xdim6_PdV_kernel_nopredict_h = -1;
extern int xdim7_PdV_kernel_nopredict;
int xdim7_PdV_kernel_nopredict_h = -1;
extern int xdim8_PdV_kernel_nopredict;
int xdim8_PdV_kernel_nopredict_h = -1;
extern int xdim9_PdV_kernel_nopredict;
int xdim9_PdV_kernel_nopredict_h = -1;
extern int xdim10_PdV_kernel_nopredict;
int xdim10_PdV_kernel_nopredict_h = -1;
extern int xdim11_PdV_kernel_nopredict;
int xdim11_PdV_kernel_nopredict_h = -1;
extern int xdim12_PdV_kernel_nopredict;
int xdim12_PdV_kernel_nopredict_h = -1;
extern int xdim13_PdV_kernel_nopredict;
int xdim13_PdV_kernel_nopredict_h = -1;

#ifdef __cplusplus
extern "C" {
#endif
void PdV_kernel_nopredict_c_wrapper(double *p_a0, double *p_a1, double *p_a2,
                                    double *p_a3, double *p_a4, double *p_a5,
                                    double *p_a6, double *p_a7, double *p_a8,
                                    double *p_a9, double *p_a10, double *p_a11,
                                    double *p_a12, double *p_a13, int x_size,
                                    int y_size);

#ifdef __cplusplus
}
#endif

// host stub function
void ops_par_loop_PdV_kernel_nopredict(char const *name, ops_block block,
                                       int dim, int *range, ops_arg arg0,
                                       ops_arg arg1, ops_arg arg2, ops_arg arg3,
                                       ops_arg arg4, ops_arg arg5, ops_arg arg6,
                                       ops_arg arg7, ops_arg arg8, ops_arg arg9,
                                       ops_arg arg10, ops_arg arg11,
                                       ops_arg arg12, ops_arg arg13) {

  ops_arg args[14] = {arg0, arg1, arg2, arg3,  arg4,  arg5,  arg6,
                      arg7, arg8, arg9, arg10, arg11, arg12, arg13};

#ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args, 14, range, 56))
    return;
#endif

  if (OPS_diags > 1) {
    ops_timing_realloc(56, "PdV_kernel_nopredict");
    OPS_kernels[56].count++;
  }

  // compute localy allocated range for the sub-block
  int start[2];
  int end[2];
  int arg_idx[2];

#ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (compute_ranges(args, 14, block, range, start, end, arg_idx) < 0)
    return;
#else
  for (int n = 0; n < 2; n++) {
    start[n] = range[2 * n];
    end[n] = range[2 * n + 1];
    arg_idx[n] = start[n];
  }
#endif

  int x_size = MAX(0, end[0] - start[0]);
  int y_size = MAX(0, end[1] - start[1]);

  xdim0 = args[0].dat->size[0];
  xdim1 = args[1].dat->size[0];
  xdim2 = args[2].dat->size[0];
  xdim3 = args[3].dat->size[0];
  xdim4 = args[4].dat->size[0];
  xdim5 = args[5].dat->size[0];
  xdim6 = args[6].dat->size[0];
  xdim7 = args[7].dat->size[0];
  xdim8 = args[8].dat->size[0];
  xdim9 = args[9].dat->size[0];
  xdim10 = args[10].dat->size[0];
  xdim11 = args[11].dat->size[0];
  xdim12 = args[12].dat->size[0];
  xdim13 = args[13].dat->size[0];

  // Timing
  double t1, t2, c1, c2;
  if (OPS_diags > 1) {
    ops_timers_core(&c2, &t2);
  }

  if (xdim0 != xdim0_PdV_kernel_nopredict_h ||
      xdim1 != xdim1_PdV_kernel_nopredict_h ||
      xdim2 != xdim2_PdV_kernel_nopredict_h ||
      xdim3 != xdim3_PdV_kernel_nopredict_h ||
      xdim4 != xdim4_PdV_kernel_nopredict_h ||
      xdim5 != xdim5_PdV_kernel_nopredict_h ||
      xdim6 != xdim6_PdV_kernel_nopredict_h ||
      xdim7 != xdim7_PdV_kernel_nopredict_h ||
      xdim8 != xdim8_PdV_kernel_nopredict_h ||
      xdim9 != xdim9_PdV_kernel_nopredict_h ||
      xdim10 != xdim10_PdV_kernel_nopredict_h ||
      xdim11 != xdim11_PdV_kernel_nopredict_h ||
      xdim12 != xdim12_PdV_kernel_nopredict_h ||
      xdim13 != xdim13_PdV_kernel_nopredict_h) {
    xdim0_PdV_kernel_nopredict = xdim0;
    xdim0_PdV_kernel_nopredict_h = xdim0;
    xdim1_PdV_kernel_nopredict = xdim1;
    xdim1_PdV_kernel_nopredict_h = xdim1;
    xdim2_PdV_kernel_nopredict = xdim2;
    xdim2_PdV_kernel_nopredict_h = xdim2;
    xdim3_PdV_kernel_nopredict = xdim3;
    xdim3_PdV_kernel_nopredict_h = xdim3;
    xdim4_PdV_kernel_nopredict = xdim4;
    xdim4_PdV_kernel_nopredict_h = xdim4;
    xdim5_PdV_kernel_nopredict = xdim5;
    xdim5_PdV_kernel_nopredict_h = xdim5;
    xdim6_PdV_kernel_nopredict = xdim6;
    xdim6_PdV_kernel_nopredict_h = xdim6;
    xdim7_PdV_kernel_nopredict = xdim7;
    xdim7_PdV_kernel_nopredict_h = xdim7;
    xdim8_PdV_kernel_nopredict = xdim8;
    xdim8_PdV_kernel_nopredict_h = xdim8;
    xdim9_PdV_kernel_nopredict = xdim9;
    xdim9_PdV_kernel_nopredict_h = xdim9;
    xdim10_PdV_kernel_nopredict = xdim10;
    xdim10_PdV_kernel_nopredict_h = xdim10;
    xdim11_PdV_kernel_nopredict = xdim11;
    xdim11_PdV_kernel_nopredict_h = xdim11;
    xdim12_PdV_kernel_nopredict = xdim12;
    xdim12_PdV_kernel_nopredict_h = xdim12;
    xdim13_PdV_kernel_nopredict = xdim13;
    xdim13_PdV_kernel_nopredict_h = xdim13;
  }

  int dat0 = (OPS_soa ? args[0].dat->type_size : args[0].dat->elem_size);
  int dat1 = (OPS_soa ? args[1].dat->type_size : args[1].dat->elem_size);
  int dat2 = (OPS_soa ? args[2].dat->type_size : args[2].dat->elem_size);
  int dat3 = (OPS_soa ? args[3].dat->type_size : args[3].dat->elem_size);
  int dat4 = (OPS_soa ? args[4].dat->type_size : args[4].dat->elem_size);
  int dat5 = (OPS_soa ? args[5].dat->type_size : args[5].dat->elem_size);
  int dat6 = (OPS_soa ? args[6].dat->type_size : args[6].dat->elem_size);
  int dat7 = (OPS_soa ? args[7].dat->type_size : args[7].dat->elem_size);
  int dat8 = (OPS_soa ? args[8].dat->type_size : args[8].dat->elem_size);
  int dat9 = (OPS_soa ? args[9].dat->type_size : args[9].dat->elem_size);
  int dat10 = (OPS_soa ? args[10].dat->type_size : args[10].dat->elem_size);
  int dat11 = (OPS_soa ? args[11].dat->type_size : args[11].dat->elem_size);
  int dat12 = (OPS_soa ? args[12].dat->type_size : args[12].dat->elem_size);
  int dat13 = (OPS_soa ? args[13].dat->type_size : args[13].dat->elem_size);

  // set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset +
              (OPS_soa ? args[0].dat->type_size : args[0].dat->elem_size) *
                  start[0] * args[0].stencil->stride[0];
  base0 = base0 +
          (OPS_soa ? args[0].dat->type_size : args[0].dat->elem_size) *
              args[0].dat->size[0] * start[1] * args[0].stencil->stride[1];
  double *p_a0 = (double *)(args[0].data + base0);

  int base1 = args[1].dat->base_offset +
              (OPS_soa ? args[1].dat->type_size : args[1].dat->elem_size) *
                  start[0] * args[1].stencil->stride[0];
  base1 = base1 +
          (OPS_soa ? args[1].dat->type_size : args[1].dat->elem_size) *
              args[1].dat->size[0] * start[1] * args[1].stencil->stride[1];
  double *p_a1 = (double *)(args[1].data + base1);

  int base2 = args[2].dat->base_offset +
              (OPS_soa ? args[2].dat->type_size : args[2].dat->elem_size) *
                  start[0] * args[2].stencil->stride[0];
  base2 = base2 +
          (OPS_soa ? args[2].dat->type_size : args[2].dat->elem_size) *
              args[2].dat->size[0] * start[1] * args[2].stencil->stride[1];
  double *p_a2 = (double *)(args[2].data + base2);

  int base3 = args[3].dat->base_offset +
              (OPS_soa ? args[3].dat->type_size : args[3].dat->elem_size) *
                  start[0] * args[3].stencil->stride[0];
  base3 = base3 +
          (OPS_soa ? args[3].dat->type_size : args[3].dat->elem_size) *
              args[3].dat->size[0] * start[1] * args[3].stencil->stride[1];
  double *p_a3 = (double *)(args[3].data + base3);

  int base4 = args[4].dat->base_offset +
              (OPS_soa ? args[4].dat->type_size : args[4].dat->elem_size) *
                  start[0] * args[4].stencil->stride[0];
  base4 = base4 +
          (OPS_soa ? args[4].dat->type_size : args[4].dat->elem_size) *
              args[4].dat->size[0] * start[1] * args[4].stencil->stride[1];
  double *p_a4 = (double *)(args[4].data + base4);

  int base5 = args[5].dat->base_offset +
              (OPS_soa ? args[5].dat->type_size : args[5].dat->elem_size) *
                  start[0] * args[5].stencil->stride[0];
  base5 = base5 +
          (OPS_soa ? args[5].dat->type_size : args[5].dat->elem_size) *
              args[5].dat->size[0] * start[1] * args[5].stencil->stride[1];
  double *p_a5 = (double *)(args[5].data + base5);

  int base6 = args[6].dat->base_offset +
              (OPS_soa ? args[6].dat->type_size : args[6].dat->elem_size) *
                  start[0] * args[6].stencil->stride[0];
  base6 = base6 +
          (OPS_soa ? args[6].dat->type_size : args[6].dat->elem_size) *
              args[6].dat->size[0] * start[1] * args[6].stencil->stride[1];
  double *p_a6 = (double *)(args[6].data + base6);

  int base7 = args[7].dat->base_offset +
              (OPS_soa ? args[7].dat->type_size : args[7].dat->elem_size) *
                  start[0] * args[7].stencil->stride[0];
  base7 = base7 +
          (OPS_soa ? args[7].dat->type_size : args[7].dat->elem_size) *
              args[7].dat->size[0] * start[1] * args[7].stencil->stride[1];
  double *p_a7 = (double *)(args[7].data + base7);

  int base8 = args[8].dat->base_offset +
              (OPS_soa ? args[8].dat->type_size : args[8].dat->elem_size) *
                  start[0] * args[8].stencil->stride[0];
  base8 = base8 +
          (OPS_soa ? args[8].dat->type_size : args[8].dat->elem_size) *
              args[8].dat->size[0] * start[1] * args[8].stencil->stride[1];
  double *p_a8 = (double *)(args[8].data + base8);

  int base9 = args[9].dat->base_offset +
              (OPS_soa ? args[9].dat->type_size : args[9].dat->elem_size) *
                  start[0] * args[9].stencil->stride[0];
  base9 = base9 +
          (OPS_soa ? args[9].dat->type_size : args[9].dat->elem_size) *
              args[9].dat->size[0] * start[1] * args[9].stencil->stride[1];
  double *p_a9 = (double *)(args[9].data + base9);

  int base10 = args[10].dat->base_offset +
               (OPS_soa ? args[10].dat->type_size : args[10].dat->elem_size) *
                   start[0] * args[10].stencil->stride[0];
  base10 = base10 +
           (OPS_soa ? args[10].dat->type_size : args[10].dat->elem_size) *
               args[10].dat->size[0] * start[1] * args[10].stencil->stride[1];
  double *p_a10 = (double *)(args[10].data + base10);

  int base11 = args[11].dat->base_offset +
               (OPS_soa ? args[11].dat->type_size : args[11].dat->elem_size) *
                   start[0] * args[11].stencil->stride[0];
  base11 = base11 +
           (OPS_soa ? args[11].dat->type_size : args[11].dat->elem_size) *
               args[11].dat->size[0] * start[1] * args[11].stencil->stride[1];
  double *p_a11 = (double *)(args[11].data + base11);

  int base12 = args[12].dat->base_offset +
               (OPS_soa ? args[12].dat->type_size : args[12].dat->elem_size) *
                   start[0] * args[12].stencil->stride[0];
  base12 = base12 +
           (OPS_soa ? args[12].dat->type_size : args[12].dat->elem_size) *
               args[12].dat->size[0] * start[1] * args[12].stencil->stride[1];
  double *p_a12 = (double *)(args[12].data + base12);

  int base13 = args[13].dat->base_offset +
               (OPS_soa ? args[13].dat->type_size : args[13].dat->elem_size) *
                   start[0] * args[13].stencil->stride[0];
  base13 = base13 +
           (OPS_soa ? args[13].dat->type_size : args[13].dat->elem_size) *
               args[13].dat->size[0] * start[1] * args[13].stencil->stride[1];
  double *p_a13 = (double *)(args[13].data + base13);

  ops_H_D_exchanges_host(args, 14);
  ops_halo_exchanges(args, 14, range);

  if (OPS_diags > 1) {
    ops_timers_core(&c1, &t1);
    OPS_kernels[56].mpi_time += t1 - t2;
  }

  PdV_kernel_nopredict_c_wrapper(p_a0, p_a1, p_a2, p_a3, p_a4, p_a5, p_a6, p_a7,
                                 p_a8, p_a9, p_a10, p_a11, p_a12, p_a13, x_size,
                                 y_size);

  if (OPS_diags > 1) {
    ops_timers_core(&c2, &t2);
    OPS_kernels[56].time += t2 - t1;
  }
  ops_set_dirtybit_host(args, 14);
  ops_set_halo_dirtybit3(&args[6], range);
  ops_set_halo_dirtybit3(&args[10], range);
  ops_set_halo_dirtybit3(&args[13], range);

  // Update kernel record
  if (OPS_diags > 1) {
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg1);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg2);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg3);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg4);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg5);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg6);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg7);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg8);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg9);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg10);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg11);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg12);
    OPS_kernels[56].transfer += ops_compute_transfer(dim, start, end, &arg13);
  }
}
