//
// auto-generated by ops.py on 2014-07-11 14:02
//

__constant__ int xdim0_advec_cell_kernel4_xdir;
__constant__ int xdim1_advec_cell_kernel4_xdir;
__constant__ int xdim2_advec_cell_kernel4_xdir;
__constant__ int xdim3_advec_cell_kernel4_xdir;
__constant__ int xdim4_advec_cell_kernel4_xdir;
__constant__ int xdim5_advec_cell_kernel4_xdir;
__constant__ int xdim6_advec_cell_kernel4_xdir;
__constant__ int xdim7_advec_cell_kernel4_xdir;
__constant__ int xdim8_advec_cell_kernel4_xdir;
__constant__ int xdim9_advec_cell_kernel4_xdir;
__constant__ int xdim10_advec_cell_kernel4_xdir;

#define OPS_ACC0(x,y) (x+xdim0_advec_cell_kernel4_xdir*(y))
#define OPS_ACC1(x,y) (x+xdim1_advec_cell_kernel4_xdir*(y))
#define OPS_ACC2(x,y) (x+xdim2_advec_cell_kernel4_xdir*(y))
#define OPS_ACC3(x,y) (x+xdim3_advec_cell_kernel4_xdir*(y))
#define OPS_ACC4(x,y) (x+xdim4_advec_cell_kernel4_xdir*(y))
#define OPS_ACC5(x,y) (x+xdim5_advec_cell_kernel4_xdir*(y))
#define OPS_ACC6(x,y) (x+xdim6_advec_cell_kernel4_xdir*(y))
#define OPS_ACC7(x,y) (x+xdim7_advec_cell_kernel4_xdir*(y))
#define OPS_ACC8(x,y) (x+xdim8_advec_cell_kernel4_xdir*(y))
#define OPS_ACC9(x,y) (x+xdim9_advec_cell_kernel4_xdir*(y))
#define OPS_ACC10(x,y) (x+xdim10_advec_cell_kernel4_xdir*(y))

//user function
__device__

inline void advec_cell_kernel4_xdir( double *density1, double *energy1,
                         const double *mass_flux_x, const double *vol_flux_x,
                         const double *pre_vol, const double *post_vol,
                         double *pre_mass, double *post_mass,
                         double *advec_vol, double *post_ener,
                         const double *ener_flux) {

  pre_mass[OPS_ACC6(0,0)] = density1[OPS_ACC0(0,0)] * pre_vol[OPS_ACC4(0,0)];
  post_mass[OPS_ACC7(0,0)] = pre_mass[OPS_ACC6(0,0)] + mass_flux_x[OPS_ACC2(0,0)] - mass_flux_x[OPS_ACC2(1,0)];
  post_ener[OPS_ACC9(0,0)] = ( energy1[OPS_ACC1(0,0)] * pre_mass[OPS_ACC6(0,0)] + ener_flux[OPS_ACC10(0,0)] - ener_flux[OPS_ACC10(1,0)])/post_mass[OPS_ACC7(0,0)];
  advec_vol[OPS_ACC8(0,0)] = pre_vol[OPS_ACC4(0,0)] + vol_flux_x[OPS_ACC3(0,0)] - vol_flux_x[OPS_ACC3(1,0)];
  density1[OPS_ACC0(0,0)] = post_mass[OPS_ACC7(0,0)]/advec_vol[OPS_ACC8(0,0)];
  energy1[OPS_ACC1(0,0)] = post_ener[OPS_ACC9(0,0)];

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


__global__ void ops_advec_cell_kernel4_xdir(
double* __restrict arg0,
double* __restrict arg1,
const double* __restrict arg2,
const double* __restrict arg3,
const double* __restrict arg4,
const double* __restrict arg5,
double* __restrict arg6,
double* __restrict arg7,
double* __restrict arg8,
double* __restrict arg9,
const double* __restrict arg10,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 + idx_y * 1 * xdim0_advec_cell_kernel4_xdir;
  arg1 += idx_x * 1 + idx_y * 1 * xdim1_advec_cell_kernel4_xdir;
  arg2 += idx_x * 1 + idx_y * 1 * xdim2_advec_cell_kernel4_xdir;
  arg3 += idx_x * 1 + idx_y * 1 * xdim3_advec_cell_kernel4_xdir;
  arg4 += idx_x * 1 + idx_y * 1 * xdim4_advec_cell_kernel4_xdir;
  arg5 += idx_x * 1 + idx_y * 1 * xdim5_advec_cell_kernel4_xdir;
  arg6 += idx_x * 1 + idx_y * 1 * xdim6_advec_cell_kernel4_xdir;
  arg7 += idx_x * 1 + idx_y * 1 * xdim7_advec_cell_kernel4_xdir;
  arg8 += idx_x * 1 + idx_y * 1 * xdim8_advec_cell_kernel4_xdir;
  arg9 += idx_x * 1 + idx_y * 1 * xdim9_advec_cell_kernel4_xdir;
  arg10 += idx_x * 1 + idx_y * 1 * xdim10_advec_cell_kernel4_xdir;

  if (idx_x < size0 && idx_y < size1) {
    advec_cell_kernel4_xdir(arg0, arg1, arg2, arg3,
                   arg4, arg5, arg6, arg7, arg8,
                   arg9, arg10);
  }

}

// host stub function
void ops_par_loop_advec_cell_kernel4_xdir(char const *name, ops_block Block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7, ops_arg arg8,
 ops_arg arg9, ops_arg arg10) {

  ops_arg args[11] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10};

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  for ( int n=0; n<2; n++ ){
    start[n] = sb->decomp_disp[n];end[n] = sb->decomp_disp[n]+sb->decomp_size[n];
    if (start[n] >= range[2*n]) {
      start[n] = 0;
    }
    else {
      start[n] = range[2*n] - start[n];
    }
    if (end[n] >= range[2*n+1]) {
      end[n] = range[2*n+1] - sb->decomp_disp[n];
    }
    else {
      end[n] = sb->decomp_size[n];
    }
  }
  #else //OPS_MPI
  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif //OPS_MPI

  int x_size = MAX(0,end[0]-start[0]);
  int y_size = MAX(0,end[1]-start[1]);

  int xdim0 = args[0].dat->size[0]*args[0].dat->dim;
  int xdim1 = args[1].dat->size[0]*args[1].dat->dim;
  int xdim2 = args[2].dat->size[0]*args[2].dat->dim;
  int xdim3 = args[3].dat->size[0]*args[3].dat->dim;
  int xdim4 = args[4].dat->size[0]*args[4].dat->dim;
  int xdim5 = args[5].dat->size[0]*args[5].dat->dim;
  int xdim6 = args[6].dat->size[0]*args[6].dat->dim;
  int xdim7 = args[7].dat->size[0]*args[7].dat->dim;
  int xdim8 = args[8].dat->size[0]*args[8].dat->dim;
  int xdim9 = args[9].dat->size[0]*args[9].dat->dim;
  int xdim10 = args[10].dat->size[0]*args[10].dat->dim;


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(10,"advec_cell_kernel4_xdir");
  ops_timers_core(&c2,&t2);

  if (OPS_kernels[10].count == 0) {
    cudaMemcpyToSymbol( xdim0_advec_cell_kernel4_xdir, &xdim0, sizeof(int) );
    cudaMemcpyToSymbol( xdim1_advec_cell_kernel4_xdir, &xdim1, sizeof(int) );
    cudaMemcpyToSymbol( xdim2_advec_cell_kernel4_xdir, &xdim2, sizeof(int) );
    cudaMemcpyToSymbol( xdim3_advec_cell_kernel4_xdir, &xdim3, sizeof(int) );
    cudaMemcpyToSymbol( xdim4_advec_cell_kernel4_xdir, &xdim4, sizeof(int) );
    cudaMemcpyToSymbol( xdim5_advec_cell_kernel4_xdir, &xdim5, sizeof(int) );
    cudaMemcpyToSymbol( xdim6_advec_cell_kernel4_xdir, &xdim6, sizeof(int) );
    cudaMemcpyToSymbol( xdim7_advec_cell_kernel4_xdir, &xdim7, sizeof(int) );
    cudaMemcpyToSymbol( xdim8_advec_cell_kernel4_xdir, &xdim8, sizeof(int) );
    cudaMemcpyToSymbol( xdim9_advec_cell_kernel4_xdir, &xdim9, sizeof(int) );
    cudaMemcpyToSymbol( xdim10_advec_cell_kernel4_xdir, &xdim10, sizeof(int) );
  }



  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 block(OPS_block_size_x,OPS_block_size_y,1);



  int dat0 = args[0].dat->elem_size;
  int dat1 = args[1].dat->elem_size;
  int dat2 = args[2].dat->elem_size;
  int dat3 = args[3].dat->elem_size;
  int dat4 = args[4].dat->elem_size;
  int dat5 = args[5].dat->elem_size;
  int dat6 = args[6].dat->elem_size;
  int dat7 = args[7].dat->elem_size;
  int dat8 = args[8].dat->elem_size;
  int dat9 = args[9].dat->elem_size;
  int dat10 = args[10].dat->elem_size;

  char *p_a[11];

  //set up initial pointers
  int base0 = dat0 * 1 * 
  (start[0] * args[0].stencil->stride[0] - args[0].dat->base[0] - args[0].dat->d_m[0]);
  base0 = base0+ dat0 *
    args[0].dat->size[0] *
    (start[1] * args[0].stencil->stride[1] - args[0].dat->base[1] - args[0].dat->d_m[1]);
  p_a[0] = (char *)args[0].data_d + base0;

  int base1 = dat1 * 1 * 
  (start[0] * args[1].stencil->stride[0] - args[1].dat->base[0] - args[1].dat->d_m[0]);
  base1 = base1+ dat1 *
    args[1].dat->size[0] *
    (start[1] * args[1].stencil->stride[1] - args[1].dat->base[1] - args[1].dat->d_m[1]);
  p_a[1] = (char *)args[1].data_d + base1;

  int base2 = dat2 * 1 * 
  (start[0] * args[2].stencil->stride[0] - args[2].dat->base[0] - args[2].dat->d_m[0]);
  base2 = base2+ dat2 *
    args[2].dat->size[0] *
    (start[1] * args[2].stencil->stride[1] - args[2].dat->base[1] - args[2].dat->d_m[1]);
  p_a[2] = (char *)args[2].data_d + base2;

  int base3 = dat3 * 1 * 
  (start[0] * args[3].stencil->stride[0] - args[3].dat->base[0] - args[3].dat->d_m[0]);
  base3 = base3+ dat3 *
    args[3].dat->size[0] *
    (start[1] * args[3].stencil->stride[1] - args[3].dat->base[1] - args[3].dat->d_m[1]);
  p_a[3] = (char *)args[3].data_d + base3;

  int base4 = dat4 * 1 * 
  (start[0] * args[4].stencil->stride[0] - args[4].dat->base[0] - args[4].dat->d_m[0]);
  base4 = base4+ dat4 *
    args[4].dat->size[0] *
    (start[1] * args[4].stencil->stride[1] - args[4].dat->base[1] - args[4].dat->d_m[1]);
  p_a[4] = (char *)args[4].data_d + base4;

  int base5 = dat5 * 1 * 
  (start[0] * args[5].stencil->stride[0] - args[5].dat->base[0] - args[5].dat->d_m[0]);
  base5 = base5+ dat5 *
    args[5].dat->size[0] *
    (start[1] * args[5].stencil->stride[1] - args[5].dat->base[1] - args[5].dat->d_m[1]);
  p_a[5] = (char *)args[5].data_d + base5;

  int base6 = dat6 * 1 * 
  (start[0] * args[6].stencil->stride[0] - args[6].dat->base[0] - args[6].dat->d_m[0]);
  base6 = base6+ dat6 *
    args[6].dat->size[0] *
    (start[1] * args[6].stencil->stride[1] - args[6].dat->base[1] - args[6].dat->d_m[1]);
  p_a[6] = (char *)args[6].data_d + base6;

  int base7 = dat7 * 1 * 
  (start[0] * args[7].stencil->stride[0] - args[7].dat->base[0] - args[7].dat->d_m[0]);
  base7 = base7+ dat7 *
    args[7].dat->size[0] *
    (start[1] * args[7].stencil->stride[1] - args[7].dat->base[1] - args[7].dat->d_m[1]);
  p_a[7] = (char *)args[7].data_d + base7;

  int base8 = dat8 * 1 * 
  (start[0] * args[8].stencil->stride[0] - args[8].dat->base[0] - args[8].dat->d_m[0]);
  base8 = base8+ dat8 *
    args[8].dat->size[0] *
    (start[1] * args[8].stencil->stride[1] - args[8].dat->base[1] - args[8].dat->d_m[1]);
  p_a[8] = (char *)args[8].data_d + base8;

  int base9 = dat9 * 1 * 
  (start[0] * args[9].stencil->stride[0] - args[9].dat->base[0] - args[9].dat->d_m[0]);
  base9 = base9+ dat9 *
    args[9].dat->size[0] *
    (start[1] * args[9].stencil->stride[1] - args[9].dat->base[1] - args[9].dat->d_m[1]);
  p_a[9] = (char *)args[9].data_d + base9;

  int base10 = dat10 * 1 * 
  (start[0] * args[10].stencil->stride[0] - args[10].dat->base[0] - args[10].dat->d_m[0]);
  base10 = base10+ dat10 *
    args[10].dat->size[0] *
    (start[1] * args[10].stencil->stride[1] - args[10].dat->base[1] - args[10].dat->d_m[1]);
  p_a[10] = (char *)args[10].data_d + base10;


  ops_H_D_exchanges_device(args, 11);
  ops_halo_exchanges(args,11,range);

  ops_timers_core(&c1,&t1);
  OPS_kernels[10].mpi_time += t1-t2;


  //call kernel wrapper function, passing in pointers to data
  ops_advec_cell_kernel4_xdir<<<grid, block >>> (  (double *)p_a[0], (double *)p_a[1],
           (double *)p_a[2], (double *)p_a[3],
           (double *)p_a[4], (double *)p_a[5],
           (double *)p_a[6], (double *)p_a[7],
           (double *)p_a[8], (double *)p_a[9],
           (double *)p_a[10],x_size, y_size);

  if (OPS_diags>1) {
    cutilSafeCall(cudaDeviceSynchronize());
  }
  ops_timers_core(&c2,&t2);
  OPS_kernels[10].time += t2-t1;
  ops_set_dirtybit_device(args, 11);
  ops_set_halo_dirtybit3(&args[0],range);
  ops_set_halo_dirtybit3(&args[1],range);
  ops_set_halo_dirtybit3(&args[6],range);
  ops_set_halo_dirtybit3(&args[7],range);
  ops_set_halo_dirtybit3(&args[8],range);
  ops_set_halo_dirtybit3(&args[9],range);

  //Update kernel record
  OPS_kernels[10].count++;
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg4);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg5);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg6);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg7);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg8);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg9);
  OPS_kernels[10].transfer += ops_compute_transfer(dim, range, &arg10);
}
