//
// auto-generated by ops.py on 2014-07-10 10:37
//

#include "./MPI_OpenMP_XeonPhi/clover_leaf_common.h"

extern int xdim0_advec_mom_kernel1_y_nonvector;
extern int xdim1_advec_mom_kernel1_y_nonvector;
extern int xdim2_advec_mom_kernel1_y_nonvector;
extern int xdim3_advec_mom_kernel1_y_nonvector;
extern int xdim4_advec_mom_kernel1_y_nonvector;

#ifdef __cplusplus
extern "C" {
#endif
void advec_mom_kernel1_y_nonvector_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  double *p_a4,
  int x_size, int y_size);

#ifdef __cplusplus
}
#endif

// host stub function
void ops_par_loop_advec_mom_kernel1_y_nonvector(char const *name, ops_block Block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3, ops_arg arg4) {

  ops_arg args[5] = { arg0, arg1, arg2, arg3, arg4};

  sub_block_list sb = OPS_sub_block_list[Block->index];
  //compute localy allocated range for the sub-block
  int start_add[2];
  int end_add[2];
  for ( int n=0; n<2; n++ ){
    start_add[n] = sb->istart[n];end_add[n] = sb->iend[n]+1;
    if (start_add[n] >= range[2*n]) {
      start_add[n] = 0;
    }
    else {
      start_add[n] = range[2*n] - start_add[n];
    }
    if (end_add[n] >= range[2*n+1]) {
      end_add[n] = range[2*n+1] - sb->istart[n];
    }
    else {
      end_add[n] = sb->sizes[n];
    }
  }


  int x_size = MAX(0,end_add[0]-start_add[0]);
  int y_size = MAX(0,end_add[1]-start_add[1]);


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(20,"advec_mom_kernel1_y_nonvector");
  ops_timers_core(&c2,&t2);

  if (OPS_kernels[20].count == 0) {
    xdim0_advec_mom_kernel1_y_nonvector = args[0].dat->block_size[0]*args[0].dat->dim;
    xdim1_advec_mom_kernel1_y_nonvector = args[1].dat->block_size[0]*args[1].dat->dim;
    xdim2_advec_mom_kernel1_y_nonvector = args[2].dat->block_size[0]*args[2].dat->dim;
    xdim3_advec_mom_kernel1_y_nonvector = args[3].dat->block_size[0]*args[3].dat->dim;
    xdim4_advec_mom_kernel1_y_nonvector = args[4].dat->block_size[0]*args[4].dat->dim;
  }

  int dat0 = args[0].dat->size;
  int dat1 = args[1].dat->size;
  int dat2 = args[2].dat->size;
  int dat3 = args[3].dat->size;
  int dat4 = args[4].dat->size;


  //set up initial pointers
  int base0 = dat0 * 1 * 
    (start_add[0] * args[0].stencil->stride[0] - args[0].dat->offset[0]);
  base0 = base0+ dat0 *
    args[0].dat->block_size[0] *
    (start_add[1] * args[0].stencil->stride[1] - args[0].dat->offset[1]);
  double *p_a0 = (double *)((char *)args[0].data + base0);

  int base1 = dat1 * 1 * 
    (start_add[0] * args[1].stencil->stride[0] - args[1].dat->offset[0]);
  base1 = base1+ dat1 *
    args[1].dat->block_size[0] *
    (start_add[1] * args[1].stencil->stride[1] - args[1].dat->offset[1]);
  double *p_a1 = (double *)((char *)args[1].data + base1);

  int base2 = dat2 * 1 * 
    (start_add[0] * args[2].stencil->stride[0] - args[2].dat->offset[0]);
  base2 = base2+ dat2 *
    args[2].dat->block_size[0] *
    (start_add[1] * args[2].stencil->stride[1] - args[2].dat->offset[1]);
  double *p_a2 = (double *)((char *)args[2].data + base2);

  int base3 = dat3 * 1 * 
    (start_add[0] * args[3].stencil->stride[0] - args[3].dat->offset[0]);
  base3 = base3+ dat3 *
    args[3].dat->block_size[0] *
    (start_add[1] * args[3].stencil->stride[1] - args[3].dat->offset[1]);
  double *p_a3 = (double *)((char *)args[3].data + base3);

  int base4 = dat4 * 1 * 
    (start_add[0] * args[4].stencil->stride[0] - args[4].dat->offset[0]);
  base4 = base4+ dat4 *
    args[4].dat->block_size[0] *
    (start_add[1] * args[4].stencil->stride[1] - args[4].dat->offset[1]);
  double *p_a4 = (double *)((char *)args[4].data + base4);


  ops_H_D_exchanges(args, 5);
  ops_halo_exchanges(args,5,range);

  ops_timers_core(&c1,&t1);
  OPS_kernels[20].mpi_time += t1-t2;

  advec_mom_kernel1_y_nonvector_c_wrapper(
    p_a0,
    p_a1,
    p_a2,
    p_a3,
    p_a4,
    x_size, y_size);

  ops_timers_core(&c2,&t2);
  OPS_kernels[20].time += t2-t1;
  ops_set_dirtybit_host(args, 5);
  ops_set_halo_dirtybit3(&args[2],range);

  //Update kernel record
  OPS_kernels[20].count++;
  OPS_kernels[20].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[20].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[20].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[20].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[20].transfer += ops_compute_transfer(dim, range, &arg4);
}
