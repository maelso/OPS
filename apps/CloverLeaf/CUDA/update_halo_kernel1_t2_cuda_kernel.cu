//
// auto-generated by ops.py on 2014-03-07 14:30
//

__constant__ int xdim0_update_halo_kernel1_t2;
__constant__ int xdim1_update_halo_kernel1_t2;
__constant__ int xdim2_update_halo_kernel1_t2;
__constant__ int xdim3_update_halo_kernel1_t2;
__constant__ int xdim4_update_halo_kernel1_t2;
__constant__ int xdim5_update_halo_kernel1_t2;
__constant__ int xdim6_update_halo_kernel1_t2;

#define OPS_ACC0(x,y) (x+xdim0_update_halo_kernel1_t2*(y))
#define OPS_ACC1(x,y) (x+xdim1_update_halo_kernel1_t2*(y))
#define OPS_ACC2(x,y) (x+xdim2_update_halo_kernel1_t2*(y))
#define OPS_ACC3(x,y) (x+xdim3_update_halo_kernel1_t2*(y))
#define OPS_ACC4(x,y) (x+xdim4_update_halo_kernel1_t2*(y))
#define OPS_ACC5(x,y) (x+xdim5_update_halo_kernel1_t2*(y))
#define OPS_ACC6(x,y) (x+xdim6_update_halo_kernel1_t2*(y))

//user function
__device__

inline void update_halo_kernel1_t2(double *density0, double *density1,
                          double *energy0, double *energy1,
                          double *pressure, double *viscosity,
                          double *soundspeed , const int* fields) {
  if(fields[FIELD_DENSITY0] == 1) density0[OPS_ACC0(0,0)] = density0[OPS_ACC0(0,-3)];
  if(fields[FIELD_DENSITY1] == 1) density1[OPS_ACC1(0,0)] = density1[OPS_ACC0(0,-3)];
  if(fields[FIELD_ENERGY0] == 1) energy0[OPS_ACC2(0,0)] = energy0[OPS_ACC0(0,-3)];
  if(fields[FIELD_ENERGY1] == 1) energy1[OPS_ACC3(0,0)] = energy1[OPS_ACC0(0,-3)];
  if(fields[FIELD_PRESSURE] == 1) pressure[OPS_ACC4(0,0)] = pressure[OPS_ACC0(0,-3)];
  if(fields[FIELD_VISCOSITY] == 1) viscosity[OPS_ACC5(0,0)] = viscosity[OPS_ACC0(0,-3)];
  if(fields[FIELD_SOUNDSPEED] == 1) soundspeed[OPS_ACC6(0,0)] = soundspeed[OPS_ACC0(0,-3)];

}



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6


__global__ void ops_update_halo_kernel1_t2(
double* __restrict arg0,
double* __restrict arg1,
double* __restrict arg2,
double* __restrict arg3,
double* __restrict arg4,
double* __restrict arg5,
double* __restrict arg6,
const int* __restrict arg7,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1 + idx_y * 1 * xdim0_update_halo_kernel1_t2;
  arg1 += idx_x * 1 + idx_y * 1 * xdim1_update_halo_kernel1_t2;
  arg2 += idx_x * 1 + idx_y * 1 * xdim2_update_halo_kernel1_t2;
  arg3 += idx_x * 1 + idx_y * 1 * xdim3_update_halo_kernel1_t2;
  arg4 += idx_x * 1 + idx_y * 1 * xdim4_update_halo_kernel1_t2;
  arg5 += idx_x * 1 + idx_y * 1 * xdim5_update_halo_kernel1_t2;
  arg6 += idx_x * 1 + idx_y * 1 * xdim6_update_halo_kernel1_t2;

  if (idx_x < size0 && idx_y < size1) {
    update_halo_kernel1_t2(arg0, arg1, arg2, arg3,
                   arg4, arg5, arg6, arg7);
  }

}

// host stub function
void ops_par_loop_update_halo_kernel1_t2(char const *name, ops_block Block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7) {

  ops_arg args[8] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7};


  int x_size = range[1]-range[0];
  int y_size = range[3]-range[2];

  int xdim0 = args[0].dat->block_size[0];
  int xdim1 = args[1].dat->block_size[0];
  int xdim2 = args[2].dat->block_size[0];
  int xdim3 = args[3].dat->block_size[0];
  int xdim4 = args[4].dat->block_size[0];
  int xdim5 = args[5].dat->block_size[0];
  int xdim6 = args[6].dat->block_size[0];


  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(44,"update_halo_kernel1_t2");
  ops_timers_core(&c1,&t1);

  if (OPS_kernels[44].count == 0) {
    cudaMemcpyToSymbol( xdim0_update_halo_kernel1_t2, &xdim0, sizeof(int) );
    cudaMemcpyToSymbol( xdim1_update_halo_kernel1_t2, &xdim1, sizeof(int) );
    cudaMemcpyToSymbol( xdim2_update_halo_kernel1_t2, &xdim2, sizeof(int) );
    cudaMemcpyToSymbol( xdim3_update_halo_kernel1_t2, &xdim3, sizeof(int) );
    cudaMemcpyToSymbol( xdim4_update_halo_kernel1_t2, &xdim4, sizeof(int) );
    cudaMemcpyToSymbol( xdim5_update_halo_kernel1_t2, &xdim5, sizeof(int) );
    cudaMemcpyToSymbol( xdim6_update_halo_kernel1_t2, &xdim6, sizeof(int) );
  }


  int *arg7h = (int *)arg7.data;

  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 block(OPS_block_size_x,OPS_block_size_y,1);

  int consts_bytes = 0;

  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));

  reallocConstArrays(consts_bytes);

  consts_bytes = 0;
  arg7.data = OPS_consts_h + consts_bytes;
  arg7.data_d = OPS_consts_d + consts_bytes;
  for (int d=0; d<NUM_FIELDS; d++) ((int *)arg7.data)[d] = arg7h[d];
  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));
  mvConstArraysToDevice(consts_bytes);

  char *p_a[8];


  //set up initial pointers
  p_a[0] = &args[0].data_d[
  + args[0].dat->size * args[0].dat->block_size[0] * ( range[2] * 1 - args[0].dat->offset[1] )
  + args[0].dat->size * ( range[0] * 1 - args[0].dat->offset[0] ) ];

  p_a[1] = &args[1].data_d[
  + args[1].dat->size * args[1].dat->block_size[0] * ( range[2] * 1 - args[1].dat->offset[1] )
  + args[1].dat->size * ( range[0] * 1 - args[1].dat->offset[0] ) ];

  p_a[2] = &args[2].data_d[
  + args[2].dat->size * args[2].dat->block_size[0] * ( range[2] * 1 - args[2].dat->offset[1] )
  + args[2].dat->size * ( range[0] * 1 - args[2].dat->offset[0] ) ];

  p_a[3] = &args[3].data_d[
  + args[3].dat->size * args[3].dat->block_size[0] * ( range[2] * 1 - args[3].dat->offset[1] )
  + args[3].dat->size * ( range[0] * 1 - args[3].dat->offset[0] ) ];

  p_a[4] = &args[4].data_d[
  + args[4].dat->size * args[4].dat->block_size[0] * ( range[2] * 1 - args[4].dat->offset[1] )
  + args[4].dat->size * ( range[0] * 1 - args[4].dat->offset[0] ) ];

  p_a[5] = &args[5].data_d[
  + args[5].dat->size * args[5].dat->block_size[0] * ( range[2] * 1 - args[5].dat->offset[1] )
  + args[5].dat->size * ( range[0] * 1 - args[5].dat->offset[0] ) ];

  p_a[6] = &args[6].data_d[
  + args[6].dat->size * args[6].dat->block_size[0] * ( range[2] * 1 - args[6].dat->offset[1] )
  + args[6].dat->size * ( range[0] * 1 - args[6].dat->offset[0] ) ];


  ops_H_D_exchanges_cuda(args, 8);


  //call kernel wrapper function, passing in pointers to data
  ops_update_halo_kernel1_t2<<<grid, block >>> (  (double *)p_a[0], (double *)p_a[1],
           (double *)p_a[2], (double *)p_a[3],
           (double *)p_a[4], (double *)p_a[5],
           (double *)p_a[6], (int *)arg7.data_d,x_size, y_size);

  if (OPS_diags>1) cudaDeviceSynchronize();
  ops_set_dirtybit_cuda(args, 8);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[44].count++;
  OPS_kernels[44].time += t2-t1;
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg4);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg5);
  OPS_kernels[44].transfer += ops_compute_transfer(dim, range, &arg6);
}
