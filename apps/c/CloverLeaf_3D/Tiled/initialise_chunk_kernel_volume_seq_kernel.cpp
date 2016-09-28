//
// auto-generated by ops.py
//
#define OPS_ACC0(x,y,z) (n_x*1+n_y*xdim0_initialise_chunk_kernel_volume*1+n_z*xdim0_initialise_chunk_kernel_volume*ydim0_initialise_chunk_kernel_volume*1+x+xdim0_initialise_chunk_kernel_volume*(y)+xdim0_initialise_chunk_kernel_volume*ydim0_initialise_chunk_kernel_volume*(z))
#define OPS_ACC1(x,y,z) (n_x*0+n_y*xdim1_initialise_chunk_kernel_volume*1+n_z*xdim1_initialise_chunk_kernel_volume*ydim1_initialise_chunk_kernel_volume*0+x+xdim1_initialise_chunk_kernel_volume*(y)+xdim1_initialise_chunk_kernel_volume*ydim1_initialise_chunk_kernel_volume*(z))
#define OPS_ACC2(x,y,z) (n_x*1+n_y*xdim2_initialise_chunk_kernel_volume*1+n_z*xdim2_initialise_chunk_kernel_volume*ydim2_initialise_chunk_kernel_volume*1+x+xdim2_initialise_chunk_kernel_volume*(y)+xdim2_initialise_chunk_kernel_volume*ydim2_initialise_chunk_kernel_volume*(z))
#define OPS_ACC3(x,y,z) (n_x*1+n_y*xdim3_initialise_chunk_kernel_volume*0+n_z*xdim3_initialise_chunk_kernel_volume*ydim3_initialise_chunk_kernel_volume*0+x+xdim3_initialise_chunk_kernel_volume*(y)+xdim3_initialise_chunk_kernel_volume*ydim3_initialise_chunk_kernel_volume*(z))
#define OPS_ACC4(x,y,z) (n_x*1+n_y*xdim4_initialise_chunk_kernel_volume*1+n_z*xdim4_initialise_chunk_kernel_volume*ydim4_initialise_chunk_kernel_volume*1+x+xdim4_initialise_chunk_kernel_volume*(y)+xdim4_initialise_chunk_kernel_volume*ydim4_initialise_chunk_kernel_volume*(z))
#define OPS_ACC5(x,y,z) (n_x*0+n_y*xdim5_initialise_chunk_kernel_volume*0+n_z*xdim5_initialise_chunk_kernel_volume*ydim5_initialise_chunk_kernel_volume*1+x+xdim5_initialise_chunk_kernel_volume*(y)+xdim5_initialise_chunk_kernel_volume*ydim5_initialise_chunk_kernel_volume*(z))
#define OPS_ACC6(x,y,z) (n_x*1+n_y*xdim6_initialise_chunk_kernel_volume*1+n_z*xdim6_initialise_chunk_kernel_volume*ydim6_initialise_chunk_kernel_volume*1+x+xdim6_initialise_chunk_kernel_volume*(y)+xdim6_initialise_chunk_kernel_volume*ydim6_initialise_chunk_kernel_volume*(z))


//user function

// host stub function
void ops_par_loop_initialise_chunk_kernel_volume_execute(ops_kernel_descriptor *desc) {
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];
  ops_arg arg3 = desc->args[3];
  ops_arg arg4 = desc->args[4];
  ops_arg arg5 = desc->args[5];
  ops_arg arg6 = desc->args[6];

  //Timing
  double t1,t2,c1,c2;

  ops_arg args[7] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,7,range,55)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(55,"initialise_chunk_kernel_volume");
    OPS_kernels[55].count++;
    ops_timers_core(&c2,&t2);
  }

  //compute locally allocated range for the sub-block
  int start[3];
  int end[3];

  for ( int n=0; n<3; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "initialise_chunk_kernel_volume");
  #endif



  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  double * __restrict__ volume = (double *)(args[0].data + base0);

  int base1 = args[1].dat->base_offset;
  const double * __restrict__ celldy = (double *)(args[1].data + base1);

  int base2 = args[2].dat->base_offset;
  double * __restrict__ xarea = (double *)(args[2].data + base2);

  int base3 = args[3].dat->base_offset;
  const double * __restrict__ celldx = (double *)(args[3].data + base3);

  int base4 = args[4].dat->base_offset;
  double * __restrict__ yarea = (double *)(args[4].data + base4);

  int base5 = args[5].dat->base_offset;
  const double * __restrict__ celldz = (double *)(args[5].data + base5);

  int base6 = args[6].dat->base_offset;
  double * __restrict__ zarea = (double *)(args[6].data + base6);


  //initialize global variable with the dimension of dats
  int xdim0_initialise_chunk_kernel_volume = args[0].dat->size[0];
  int ydim0_initialise_chunk_kernel_volume = args[0].dat->size[1];
  int xdim1_initialise_chunk_kernel_volume = args[1].dat->size[0];
  int ydim1_initialise_chunk_kernel_volume = args[1].dat->size[1];
  int xdim2_initialise_chunk_kernel_volume = args[2].dat->size[0];
  int ydim2_initialise_chunk_kernel_volume = args[2].dat->size[1];
  int xdim3_initialise_chunk_kernel_volume = args[3].dat->size[0];
  int ydim3_initialise_chunk_kernel_volume = args[3].dat->size[1];
  int xdim4_initialise_chunk_kernel_volume = args[4].dat->size[0];
  int ydim4_initialise_chunk_kernel_volume = args[4].dat->size[1];
  int xdim5_initialise_chunk_kernel_volume = args[5].dat->size[0];
  int ydim5_initialise_chunk_kernel_volume = args[5].dat->size[1];
  int xdim6_initialise_chunk_kernel_volume = args[6].dat->size[0];
  int ydim6_initialise_chunk_kernel_volume = args[6].dat->size[1];

  //Halo Exchanges
  ops_H_D_exchanges_host(args, 7);
  ops_halo_exchanges(args,7,range);
  ops_H_D_exchanges_host(args, 7);

  if (OPS_diags > 1) {
    ops_timers_core(&c1,&t1);
    OPS_kernels[55].mpi_time += t1-t2;
  }

  #pragma omp parallel for collapse(2)
  for ( int n_z=start[2]; n_z<end[2]; n_z++ ){
    for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
      #ifdef intel
      #pragma omp simd
      #else
      #pragma simd
      #endif
      for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
        

  double d_x, d_y, d_z;

  d_x = (grid.xmax - grid.xmin)/(double)grid.x_cells;
  d_y = (grid.ymax - grid.ymin)/(double)grid.y_cells;
  d_z = (grid.zmax - grid.zmin)/(double)grid.z_cells;

  volume[OPS_ACC0(0,0,0)] = d_x*d_y*d_z;
  xarea[OPS_ACC2(0,0,0)] = celldy[OPS_ACC1(0,0,0)]*celldz[OPS_ACC5(0,0,0)];
  yarea[OPS_ACC4(0,0,0)] = celldx[OPS_ACC3(0,0,0)]*celldz[OPS_ACC5(0,0,0)];
  zarea[OPS_ACC6(0,0,0)] = celldx[OPS_ACC3(0,0,0)]*celldy[OPS_ACC1(0,0,0)];

      }
    }
  }
  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[55].time += t2-t1;
  }
  ops_set_dirtybit_host(args, 7);
  ops_set_halo_dirtybit3(&args[0],range);
  ops_set_halo_dirtybit3(&args[2],range);
  ops_set_halo_dirtybit3(&args[4],range);
  ops_set_halo_dirtybit3(&args[6],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c1,&t1);
    OPS_kernels[55].mpi_time += t1-t2;
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg1);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg2);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg3);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg4);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg5);
    OPS_kernels[55].transfer += ops_compute_transfer(dim, start, end, &arg6);
  }
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6


void ops_par_loop_initialise_chunk_kernel_volume(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6) {
  ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->index = 55;
  #ifdef OPS_MPI
  sub_block_list sb = OPS_sub_block_list[block->index];
  if (!sb->owned) return;
  for ( int n=0; n<3; n++ ){
    desc->range[2*n] = sb->decomp_disp[n];desc->range[2*n+1] = sb->decomp_disp[n]+sb->decomp_size[n];
    if (desc->range[2*n] >= range[2*n]) {
      desc->range[2*n] = 0;
    }
    else {
      desc->range[2*n] = range[2*n] - desc->range[2*n];
    }
    if (sb->id_m[n]==MPI_PROC_NULL && range[2*n] < 0) desc->range[2*n] = range[2*n];
    if (desc->range[2*n+1] >= range[2*n+1]) {
      desc->range[2*n+1] = range[2*n+1] - sb->decomp_disp[n];
    }
    else {
      desc->range[2*n+1] = sb->decomp_size[n];
    }
    if (sb->id_p[n]==MPI_PROC_NULL && (range[2*n+1] > sb->decomp_disp[n]+sb->decomp_size[n]))
      desc->range[2*n+1] += (range[2*n+1]-sb->decomp_disp[n]-sb->decomp_size[n]);
  }
  #else //OPS_MPI
  for ( int i=0; i<6; i++ ){
    desc->range[i] = range[i];
  }
  #endif //OPS_MPI
  desc->nargs = 7;
  desc->args = (ops_arg*)malloc(7*sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->args[1] = arg1;
  desc->args[2] = arg2;
  desc->args[3] = arg3;
  desc->args[4] = arg4;
  desc->args[5] = arg5;
  desc->args[6] = arg6;
  desc->function = ops_par_loop_initialise_chunk_kernel_volume_execute;
  ops_enqueue_kernel(desc);
  }
