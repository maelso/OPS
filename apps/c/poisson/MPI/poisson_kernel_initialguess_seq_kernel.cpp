//
// auto-generated by ops.py
//

//user function
inline void poisson_kernel_initialguess(double *u) {
  u[OPS_ACC0(0,0)] = 0.0;
}





// host stub function
void ops_par_loop_poisson_kernel_initialguess_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];

  //Timing
  double t1,t2,c1,c2;

  char *p_a[1];
  int  offs[1][2];
  ops_arg args[1] = { arg0};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,1,range,1)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(1,"poisson_kernel_initialguess");
    OPS_kernels[1].count++;
    ops_timers_core(&c2,&t2);
  }

  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];

  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "poisson_kernel_initialguess");
  #endif



  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  double *p_a0 = (double *)(args[0].data + base0);


  //initialize global variable with the dimension of dats
  xdim0 = args[0].dat->size[0];

  //Halo Exchanges
  ops_H_D_exchanges_host(args, 1);
  ops_halo_exchanges(args,1,range);
  ops_H_D_exchanges_host(args, 1);

  if (OPS_diags > 1) {
    ops_timers_core(&c1,&t1);
    OPS_kernels[1].mpi_time += t1-t2;
  }

  int n_x;
  for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
    for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
      poisson_kernel_initialguess(  p_a0 + n_x*1*1 + n_y*xdim0*1*1 );

    }
  }
  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[1].time += t2-t1;
  }
  ops_set_dirtybit_host(args, 1);
  ops_set_halo_dirtybit3(&args[0],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c1,&t1);
    OPS_kernels[1].mpi_time += t1-t2;
    OPS_kernels[1].transfer += ops_compute_transfer(dim, start, end, &arg0);
  }
}

void ops_par_loop_poisson_kernel_initialguess(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0) {
ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
desc->name = name;
desc->block = block;
desc->dim = dim;
#ifdef OPS_MPI
sub_block_list sb = OPS_sub_block_list[block->index];
if (!sb->owned) return;
for ( int n=0; n<2; n++ ){
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
for ( int i=0; i<4; i++ ){
  desc->range[i] = range[i];
}
#endif //OPS_MPI
desc->nargs = 1;
desc->args = (ops_arg*)malloc(1*sizeof(ops_arg));
desc->args[0] = arg0;
desc->function = ops_par_loop_poisson_kernel_initialguess_execute;
ops_enqueue_kernel(desc);
}
