//
// auto-generated by ops.py
//
__constant__ int xdim0_update_halo_kernel4_minus_4_b;
int xdim0_update_halo_kernel4_minus_4_b_h = -1;
int ydim0_update_halo_kernel4_minus_4_b_h = -1;
__constant__ int xdim1_update_halo_kernel4_minus_4_b;
int xdim1_update_halo_kernel4_minus_4_b_h = -1;
int ydim1_update_halo_kernel4_minus_4_b_h = -1;

#define OPS_ACC0(x,y) (x+xdim0_update_halo_kernel4_minus_4_b*(y))
#define OPS_ACC1(x,y) (x+xdim1_update_halo_kernel4_minus_4_b*(y))

//user function
__device__

inline void update_halo_kernel4_minus_4_b(double *vol_flux_y, double *mass_flux_y, const int* fields) {
  if(fields[FIELD_VOL_FLUX_Y] == 1) vol_flux_y[OPS_ACC0(0,0)] = -(vol_flux_y[OPS_ACC0(0,-4)]);
  if(fields[FIELD_MASS_FLUX_Y] == 1) mass_flux_y[OPS_ACC1(0,0)] = -(mass_flux_y[OPS_ACC1(0,-4)]);
}



#undef OPS_ACC0
#undef OPS_ACC1


__global__ void ops_update_halo_kernel4_minus_4_b(
double* __restrict arg0,
double* __restrict arg1,
const int* __restrict arg2,
int size0,
int size1 ){


  int idx_y = blockDim.y * blockIdx.y + threadIdx.y;
  int idx_x = blockDim.x * blockIdx.x + threadIdx.x;

  arg0 += idx_x * 1*1 + idx_y * 1*1 * xdim0_update_halo_kernel4_minus_4_b;
  arg1 += idx_x * 1*1 + idx_y * 1*1 * xdim1_update_halo_kernel4_minus_4_b;

  if (idx_x < size0 && idx_y < size1) {
    update_halo_kernel4_minus_4_b(arg0, arg1, arg2);
  }

}

// host stub function
void ops_par_loop_update_halo_kernel4_minus_4_b(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {

  //Timing
  double t1,t2,c1,c2;

  ops_arg args[3] = { arg0, arg1, arg2};


  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,3,range,76)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(76,"update_halo_kernel4_minus_4_b");
    OPS_kernels[76].count++;
    ops_timers_core(&c1,&t1);
  }

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<2; n++ ){
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
  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif //OPS_MPI

  int x_size = MAX(0,end[0]-start[0]);
  int y_size = MAX(0,end[1]-start[1]);

  int xdim0 = args[0].dat->size[0];
  int xdim1 = args[1].dat->size[0];

  if (xdim0 != xdim0_update_halo_kernel4_minus_4_b_h || xdim1 != xdim1_update_halo_kernel4_minus_4_b_h) {
    cudaMemcpyToSymbol( xdim0_update_halo_kernel4_minus_4_b, &xdim0, sizeof(int) );
    xdim0_update_halo_kernel4_minus_4_b_h = xdim0;
    cudaMemcpyToSymbol( xdim1_update_halo_kernel4_minus_4_b, &xdim1, sizeof(int) );
    xdim1_update_halo_kernel4_minus_4_b_h = xdim1;
  }


  int *arg2h = (int *)arg2.data;

  dim3 grid( (x_size-1)/OPS_block_size_x+ 1, (y_size-1)/OPS_block_size_y + 1, 1);
  dim3 tblock(OPS_block_size_x,OPS_block_size_y,1);

  int consts_bytes = 0;

  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));

  reallocConstArrays(consts_bytes);

  consts_bytes = 0;
  arg2.data = OPS_consts_h + consts_bytes;
  arg2.data_d = OPS_consts_d + consts_bytes;
  for (int d=0; d<NUM_FIELDS; d++) ((int *)arg2.data)[d] = arg2h[d];
  consts_bytes += ROUND_UP(NUM_FIELDS*sizeof(int));
  mvConstArraysToDevice(consts_bytes);
  int dat0 = args[0].dat->elem_size;
  int dat1 = args[1].dat->elem_size;

  char *p_a[3];

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
  p_a[0] = (char *)args[0].data_d + base0;

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
  p_a[1] = (char *)args[1].data_d + base1;


  ops_H_D_exchanges_device(args, 3);
  ops_halo_exchanges(args,3,range);

  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[76].mpi_time += t2-t1;
  }


  //call kernel wrapper function, passing in pointers to data
  ops_update_halo_kernel4_minus_4_b<<<grid, tblock >>> (  (double *)p_a[0], (double *)p_a[1],
           (int *)arg2.data_d,x_size, y_size);

  if (OPS_diags>1) {
    cutilSafeCall(cudaDeviceSynchronize());
    ops_timers_core(&c1,&t1);
    OPS_kernels[76].time += t1-t2;
  }

  ops_set_dirtybit_device(args, 3);
  ops_set_halo_dirtybit3(&args[0],range);
  ops_set_halo_dirtybit3(&args[1],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c2,&t2);
    OPS_kernels[76].mpi_time += t2-t1;
    OPS_kernels[76].transfer += ops_compute_transfer(dim, range, &arg0);
    OPS_kernels[76].transfer += ops_compute_transfer(dim, range, &arg1);
  }
}
