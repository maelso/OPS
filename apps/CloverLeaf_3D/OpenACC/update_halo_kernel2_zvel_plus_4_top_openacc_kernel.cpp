//
// auto-generated by ops.py
//
#include "./OpenACC/clover_leaf_common.h"

#define OPS_GPU

extern int xdim0_update_halo_kernel2_zvel_plus_4_top;
int xdim0_update_halo_kernel2_zvel_plus_4_top_h = -1;
extern int ydim0_update_halo_kernel2_zvel_plus_4_top;
int ydim0_update_halo_kernel2_zvel_plus_4_top_h = -1;
extern int xdim1_update_halo_kernel2_zvel_plus_4_top;
int xdim1_update_halo_kernel2_zvel_plus_4_top_h = -1;
extern int ydim1_update_halo_kernel2_zvel_plus_4_top;
int ydim1_update_halo_kernel2_zvel_plus_4_top_h = -1;

#ifdef __cplusplus
extern "C" {
#endif
void update_halo_kernel2_zvel_plus_4_top_c_wrapper(
  double *p_a0,
  double *p_a1,
  int *p_a2,
  int x_size, int y_size, int z_size);

#ifdef __cplusplus
}
#endif

// host stub function
void ops_par_loop_update_halo_kernel2_zvel_plus_4_top(char const *name, ops_block Block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {

  ops_arg args[3] = { arg0, arg1, arg2};


  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,3,range,95)) return;
  #endif

  ops_timing_realloc(95,"update_halo_kernel2_zvel_plus_4_top");
  OPS_kernels[95].count++;

  //compute localy allocated range for the sub-block
  int start[3];
  int end[3];
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<3; n++ ){
    start[n] = sb->decomp_disp[n];end[n] = sb->decomp_disp[n]+sb->decomp_size[n];
    if (start[n] >= range[2*n]) {
      start[n] = 0;
    }
    else {
      start[n] = range[2*n] - start[n];
    }
    if (sb->id_m[n]==MPI_PROC_NULL && range[2*n] < 0) start[n] = range[2*n];
    if (end[n] >= range[2*n+1]) {
      end[n] = range[2*n+1] - sb->decomp_disp[n];
    }
    else {
      end[n] = sb->decomp_size[n];
    }
    if (sb->id_p[n]==MPI_PROC_NULL && (range[2*n+1] > sb->decomp_disp[n]+sb->decomp_size[n]))
      end[n] += (range[2*n+1]-sb->decomp_disp[n]-sb->decomp_size[n]);
  }
  #else //OPS_MPI
  for ( int n=0; n<3; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif //OPS_MPI

  int x_size = MAX(0,end[0]-start[0]);
  int y_size = MAX(0,end[1]-start[1]);
  int z_size = MAX(0,end[2]-start[2]);


  xdim0 = args[0].dat->size[0]*args[0].dat->dim;
  ydim0 = args[0].dat->size[1];
  xdim1 = args[1].dat->size[0]*args[1].dat->dim;
  ydim1 = args[1].dat->size[1];

  //Timing
  double t1,t2,c1,c2;
  ops_timers_core(&c2,&t2);

  if (xdim0 != xdim0_update_halo_kernel2_zvel_plus_4_top_h || ydim0 != ydim0_update_halo_kernel2_zvel_plus_4_top_h || xdim1 != xdim1_update_halo_kernel2_zvel_plus_4_top_h || ydim1 != ydim1_update_halo_kernel2_zvel_plus_4_top_h) {
    xdim0_update_halo_kernel2_zvel_plus_4_top = xdim0;
    xdim0_update_halo_kernel2_zvel_plus_4_top_h = xdim0;
    ydim0_update_halo_kernel2_zvel_plus_4_top = ydim0;
    ydim0_update_halo_kernel2_zvel_plus_4_top_h = ydim0;
    xdim1_update_halo_kernel2_zvel_plus_4_top = xdim1;
    xdim1_update_halo_kernel2_zvel_plus_4_top_h = xdim1;
    ydim1_update_halo_kernel2_zvel_plus_4_top = ydim1;
    ydim1_update_halo_kernel2_zvel_plus_4_top_h = ydim1;
  }

  int dat0 = args[0].dat->elem_size;
  int dat1 = args[1].dat->elem_size;

  int *arg2h = (int *)arg2.data;
  //Upload large globals
  int consts_bytes = 0;
  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));
  reallocConstArrays(consts_bytes);
  consts_bytes = 0;
  args[2].data = OPS_consts_h + consts_bytes;
  args[2].data_d = OPS_consts_d + consts_bytes;
  for (int d=0; d<NUM_FIELDS; d++) ((int *)args[2].data)[d] = arg2h[d];
  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));
  mvConstArraysToDevice(consts_bytes);

  //set up initial pointers
  int d_m[OPS_MAX_DIM];
  #ifdef OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[0].dat->d_m[d] + OPS_sub_dat_list[args[0].dat->index]->d_im[d];
  #else //OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[0].dat->d_m[d];
  #endif //OPS_MPI
  int base0 = dat0 * 1 * 
    (start[0] * args[0].stencil->stride[0] - args[0].dat->base[0] - d_m[0]);
  base0 = base0+ dat0 *
    args[0].dat->size[0] *
    (start[1] * args[0].stencil->stride[1] - args[0].dat->base[1] - d_m[1]);
  base0 = base0+ dat0 *
    args[0].dat->size[0] *
    args[0].dat->size[1] *
    (start[2] * args[0].stencil->stride[2] - args[0].dat->base[2] - d_m[2]);
  #ifdef OPS_GPU
  double *p_a0 = (double *)((char *)args[0].data_d + base0);
  #else
  double *p_a0 = (double *)((char *)args[0].data + base0);
  #endif

  #ifdef OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[1].dat->d_m[d] + OPS_sub_dat_list[args[1].dat->index]->d_im[d];
  #else //OPS_MPI
  for (int d = 0; d < dim; d++) d_m[d] = args[1].dat->d_m[d];
  #endif //OPS_MPI
  int base1 = dat1 * 1 * 
    (start[0] * args[1].stencil->stride[0] - args[1].dat->base[0] - d_m[0]);
  base1 = base1+ dat1 *
    args[1].dat->size[0] *
    (start[1] * args[1].stencil->stride[1] - args[1].dat->base[1] - d_m[1]);
  base1 = base1+ dat1 *
    args[1].dat->size[0] *
    args[1].dat->size[1] *
    (start[2] * args[1].stencil->stride[2] - args[1].dat->base[2] - d_m[2]);
  #ifdef OPS_GPU
  double *p_a1 = (double *)((char *)args[1].data_d + base1);
  #else
  double *p_a1 = (double *)((char *)args[1].data + base1);
  #endif

  #ifdef OPS_GPU
  int *p_a2 = (int *)args[2].data_d;
  #else
  int *p_a2 = arg2h;
  #endif

  #ifdef OPS_GPU
  ops_H_D_exchanges_device(args, 3);
  #else
  ops_H_D_exchanges_host(args, 3);
  #endif
  ops_halo_exchanges(args,3,range);

  ops_timers_core(&c1,&t1);
  OPS_kernels[95].mpi_time += t1-t2;

  update_halo_kernel2_zvel_plus_4_top_c_wrapper(
    p_a0,
    p_a1,
    p_a2,
    x_size, y_size, z_size);

  ops_timers_core(&c2,&t2);
  OPS_kernels[95].time += t2-t1;
  #ifdef OPS_GPU
  ops_set_dirtybit_device(args, 3);
  #else
  ops_set_dirtybit_host(args, 3);
  #endif
  ops_set_halo_dirtybit3(&args[0],range);
  ops_set_halo_dirtybit3(&args[1],range);

  //Update kernel record
  OPS_kernels[95].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[95].transfer += ops_compute_transfer(dim, range, &arg1);
}
