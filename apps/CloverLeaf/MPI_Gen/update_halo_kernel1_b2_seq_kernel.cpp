//
// auto-generated by ops.py on 2014-02-24 14:30
//

//user function
#include "update_halo_kernel.h"

// host stub function
void ops_par_loop_update_halo_kernel1_b2(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6) {

  char *p_a[7];
  int  offs[7][2];
  ops_arg args[7] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6};


  sub_block_list sb = OPS_sub_block_list[block->index];
  //compute localy allocated range for the sub-block
  int ndim = sb->ndim;
  int start[ndim*7];
  int end[ndim*7];

  int s[ndim];
  int e[ndim];

  for ( int n=0; n<ndim; n++ ){
    s[n] = sb->istart[n];e[n] = sb->iend[n]+1;
    if (s[n] >= range[2*n]) {
      s[n] = 0;
    }
    else {
      s[n] = range[2*n] - s[n];
    }
    if (e[n] >= range[2*n+1]) {
      e[n] = range[2*n+1] - sb->istart[n];
    }
    else {
      e[n] = sb->sizes[n];
    }
  }

  for ( int i=0; i<7; i++ ){
    for ( int n=0; n<ndim; n++ ){
      start[i*ndim+n] = s[n];
      end[i*ndim+n]   = e[n];
    }
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "update_halo_kernel1_b2");
  #endif

  offs[0][0] = args[0].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[0][n] = off2(ndim, n, &start[0*ndim],
    &end[0*ndim],args[0].dat->block_size, args[0].stencil->stride);
  }
  offs[1][0] = args[1].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[1][n] = off2(ndim, n, &start[1*ndim],
    &end[1*ndim],args[1].dat->block_size, args[1].stencil->stride);
  }
  offs[2][0] = args[2].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[2][n] = off2(ndim, n, &start[2*ndim],
    &end[2*ndim],args[2].dat->block_size, args[2].stencil->stride);
  }
  offs[3][0] = args[3].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[3][n] = off2(ndim, n, &start[3*ndim],
    &end[3*ndim],args[3].dat->block_size, args[3].stencil->stride);
  }
  offs[4][0] = args[4].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[4][n] = off2(ndim, n, &start[4*ndim],
    &end[4*ndim],args[4].dat->block_size, args[4].stencil->stride);
  }
  offs[5][0] = args[5].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[5][n] = off2(ndim, n, &start[5*ndim],
    &end[5*ndim],args[5].dat->block_size, args[5].stencil->stride);
  }
  offs[6][0] = args[6].stencil->stride[0]*1;  //unit step in x dimension
  for ( int n=1; n<ndim; n++ ){
    offs[6][n] = off2(ndim, n, &start[6*ndim],
    &end[6*ndim],args[6].dat->block_size, args[6].stencil->stride);
  }


  //set up initial pointers
  p_a[0] = (char *)args[0].data
  + address2(ndim, args[0].dat->size, &start[0*ndim],
  args[0].dat->block_size, args[0].stencil->stride, args[0].dat->offset);
  ops_exchange_halo(&args[0],2);

  //set up initial pointers
  p_a[1] = (char *)args[1].data
  + address2(ndim, args[1].dat->size, &start[1*ndim],
  args[1].dat->block_size, args[1].stencil->stride, args[1].dat->offset);
  ops_exchange_halo(&args[1],2);

  //set up initial pointers
  p_a[2] = (char *)args[2].data
  + address2(ndim, args[2].dat->size, &start[2*ndim],
  args[2].dat->block_size, args[2].stencil->stride, args[2].dat->offset);
  ops_exchange_halo(&args[2],2);

  //set up initial pointers
  p_a[3] = (char *)args[3].data
  + address2(ndim, args[3].dat->size, &start[3*ndim],
  args[3].dat->block_size, args[3].stencil->stride, args[3].dat->offset);
  ops_exchange_halo(&args[3],2);

  //set up initial pointers
  p_a[4] = (char *)args[4].data
  + address2(ndim, args[4].dat->size, &start[4*ndim],
  args[4].dat->block_size, args[4].stencil->stride, args[4].dat->offset);
  ops_exchange_halo(&args[4],2);

  //set up initial pointers
  p_a[5] = (char *)args[5].data
  + address2(ndim, args[5].dat->size, &start[5*ndim],
  args[5].dat->block_size, args[5].stencil->stride, args[5].dat->offset);
  ops_exchange_halo(&args[5],2);

  //set up initial pointers
  p_a[6] = (char *)args[6].data
  + address2(ndim, args[6].dat->size, &start[6*ndim],
  args[6].dat->block_size, args[6].stencil->stride, args[6].dat->offset);
  ops_exchange_halo(&args[6],2);


  int off0_1 = offs[0][0];
  int off0_2 = offs[0][1];
  int dat0 = args[0].dat->size;
  int off1_1 = offs[1][0];
  int off1_2 = offs[1][1];
  int dat1 = args[1].dat->size;
  int off2_1 = offs[2][0];
  int off2_2 = offs[2][1];
  int dat2 = args[2].dat->size;
  int off3_1 = offs[3][0];
  int off3_2 = offs[3][1];
  int dat3 = args[3].dat->size;
  int off4_1 = offs[4][0];
  int off4_2 = offs[4][1];
  int dat4 = args[4].dat->size;
  int off5_1 = offs[5][0];
  int off5_2 = offs[5][1];
  int dat5 = args[5].dat->size;
  int off6_1 = offs[6][0];
  int off6_2 = offs[6][1];
  int dat6 = args[6].dat->size;

  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(0,"update_halo_kernel1_b2");
  ops_timers_core(&c1,&t1);

  xdim0 = args[0].dat->block_size[0];
  xdim1 = args[1].dat->block_size[0];
  xdim2 = args[2].dat->block_size[0];
  xdim3 = args[3].dat->block_size[0];
  xdim4 = args[4].dat->block_size[0];
  xdim5 = args[5].dat->block_size[0];
  xdim6 = args[6].dat->block_size[0];

  int n_x;
  for ( int n_y=s[1]; n_y<e[1]; n_y++ ){
    for ( int n_x=s[0]; n_x<s[0]+(e[0]-s[0])/SIMD_VEC; n_x++ ){
        //call kernel function, passing in pointers to data -vectorised
        #pragma simd
        for ( int i=0; i<SIMD_VEC; i++ ){
          update_halo_kernel1_b2(  (double *)p_a[0]+ i*1, (double *)p_a[1]+ i*1, (double *)p_a[2]+ i*1,
           (double *)p_a[3]+ i*1, (double *)p_a[4]+ i*1, (double *)p_a[5]+ i*1, (double *)p_a[6]+ i*1 );

        }

        //shift pointers to data x direction
        p_a[0]= p_a[0] + (dat0 * off0_1)*SIMD_VEC;
        p_a[1]= p_a[1] + (dat1 * off1_1)*SIMD_VEC;
        p_a[2]= p_a[2] + (dat2 * off2_1)*SIMD_VEC;
        p_a[3]= p_a[3] + (dat3 * off3_1)*SIMD_VEC;
        p_a[4]= p_a[4] + (dat4 * off4_1)*SIMD_VEC;
        p_a[5]= p_a[5] + (dat5 * off5_1)*SIMD_VEC;
        p_a[6]= p_a[6] + (dat6 * off6_1)*SIMD_VEC;
      }

      for ( int n_x=s[0]+((e[0]-s[0])/SIMD_VEC)*SIMD_VEC; n_x<e[0]; n_x++ ){
          //call kernel function, passing in pointers to data - remainder
          update_halo_kernel1_b2(  (double *)p_a[0], (double *)p_a[1], (double *)p_a[2],
           (double *)p_a[3], (double *)p_a[4], (double *)p_a[5], (double *)p_a[6] );


          //shift pointers to data x direction
          p_a[0]= p_a[0] + (dat0 * off0_1);
          p_a[1]= p_a[1] + (dat1 * off1_1);
          p_a[2]= p_a[2] + (dat2 * off2_1);
          p_a[3]= p_a[3] + (dat3 * off3_1);
          p_a[4]= p_a[4] + (dat4 * off4_1);
          p_a[5]= p_a[5] + (dat5 * off5_1);
          p_a[6]= p_a[6] + (dat6 * off6_1);
        }

        //shift pointers to data y direction
        p_a[0]= p_a[0] + (dat0 * off0_2);
        p_a[1]= p_a[1] + (dat1 * off1_2);
        p_a[2]= p_a[2] + (dat2 * off2_2);
        p_a[3]= p_a[3] + (dat3 * off3_2);
        p_a[4]= p_a[4] + (dat4 * off4_2);
        p_a[5]= p_a[5] + (dat5 * off5_2);
        p_a[6]= p_a[6] + (dat6 * off6_2);
      }
      ops_set_halo_dirtybit(&args[0]);
      ops_set_halo_dirtybit(&args[1]);
      ops_set_halo_dirtybit(&args[2]);
      ops_set_halo_dirtybit(&args[3]);
      ops_set_halo_dirtybit(&args[4]);
      ops_set_halo_dirtybit(&args[5]);
      ops_set_halo_dirtybit(&args[6]);

      //Update kernel record
      ops_timers_core(&c2,&t2);
      OPS_kernels[0].count++;
      OPS_kernels[0].time += t2-t1;
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg0);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg1);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg2);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg3);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg4);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg5);
      OPS_kernels[0].transfer += ops_compute_transfer(dim, range, &arg6);
    }
