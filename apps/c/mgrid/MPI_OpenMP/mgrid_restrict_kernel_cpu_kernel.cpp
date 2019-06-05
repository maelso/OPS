//
// auto-generated by ops.py
//
#define OPS_ACC0(x,y) (n_x*args[0].stencil->mgrid_stride[0] + x + (n_y*args[0].stencil->mgrid_stride[1]+(y))*xdim0_mgrid_restrict_kernel)
#define OPS_ACC1(x,y) (n_x*1 + x + (n_y*1+(y))*xdim1_mgrid_restrict_kernel)

//user function

// host stub function
#ifndef OPS_LAZY
void ops_par_loop_mgrid_restrict_kernel(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {
#else
void ops_par_loop_mgrid_restrict_kernel_execute(ops_kernel_descriptor *desc) {
  ops_block block = desc->block;
  int dim = desc->dim;
  int *range = desc->range;
  ops_arg arg0 = desc->args[0];
  ops_arg arg1 = desc->args[1];
  ops_arg arg2 = desc->args[2];
  #endif

  //Timing
  double __t1,__t2,__c1,__c2;

  ops_arg args[3] = { arg0, arg1, arg2};



  #if defined(CHECKPOINTING) && !defined(OPS_LAZY)
  if (!ops_checkpointing_before(args,3,range,6)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(6,"mgrid_restrict_kernel");
    OPS_kernels[6].count++;
    ops_timers_core(&__c2,&__t2);
  }

  #ifdef OPS_DEBUG
  ops_register_args(args, "mgrid_restrict_kernel");
  #endif


  //compute locally allocated range for the sub-block
  int start[2];
  int end[2];
  int arg_idx[2];
  #if defined(OPS_LAZY) || !defined(OPS_MPI)
  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #else
  if (compute_ranges(args, 3,block, range, start, end, arg_idx) < 0) return;
  #endif

  #ifdef OPS_MPI
  arg_idx[0] -= start[0];
  arg_idx[1] -= start[1];
#else
  arg_idx[0] = 0;
  arg_idx[1] = 0;
#endif // OPS_MPI

  //initialize global variable with the dimension of dats
  int xdim0_mgrid_restrict_kernel = args[0].dat->size[0];
  int xdim1_mgrid_restrict_kernel = args[1].dat->size[0];

  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset;
  const double * __restrict__ fine = (double *)(args[0].data + base0);
  #ifdef OPS_MPI
  sub_dat_list sd0 = OPS_sub_dat_list[args[0].dat->index];
  fine += arg_idx[0]*args[0].stencil->mgrid_stride[0] - sd0->decomp_disp[0] + args[0].dat->d_m[0];
  fine += (arg_idx[1]*args[0].stencil->mgrid_stride[1] - sd0->decomp_disp[1] + args[0].dat->d_m[1])*xdim0_mgrid_restrict_kernel;
  #endif

  int base1 = args[1].dat->base_offset;
  double * __restrict__ coarse = (double *)(args[1].data + base1);




  #ifndef OPS_LAZY
  //Halo Exchanges
  ops_H_D_exchanges_host(args, 3);
  ops_halo_exchanges(args,3,range);
  ops_H_D_exchanges_host(args, 3);
  #endif

  if (OPS_diags > 1) {
    ops_timers_core(&__c1,&__t1);
    OPS_kernels[6].mpi_time += __t1-__t2;
  }

  #pragma omp parallel for
  for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
    #ifdef __INTEL_COMPILER
    #pragma loop_count(10000)
    #pragma omp simd aligned(fine,coarse)
    #else
    #pragma simd
    #endif
    for ( int n_x=start[0]; n_x<end[0]; n_x++ ){
      int idx[] = {arg_idx[0]+n_x, arg_idx[1]+n_y};
      

  coarse[OPS_ACC1(0,0)] = fine[OPS_ACC0(0,0)];

    }
  }
  if (OPS_diags > 1) {
    ops_timers_core(&__c2,&__t2);
    OPS_kernels[6].time += __t2-__t1;
  }
  #ifndef OPS_LAZY
  ops_set_dirtybit_host(args, 3);
  ops_set_halo_dirtybit3(&args[1],range);
  #endif

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&__c1,&__t1);
    OPS_kernels[6].mpi_time += __t1-__t2;
    OPS_kernels[6].transfer += ops_compute_transfer(dim, start, end, &arg0);
    OPS_kernels[6].transfer += ops_compute_transfer(dim, start, end, &arg1);
  }
}
#undef OPS_ACC0
#undef OPS_ACC1


#ifdef OPS_LAZY
void ops_par_loop_mgrid_restrict_kernel(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2) {
  ops_kernel_descriptor *desc = (ops_kernel_descriptor *)malloc(sizeof(ops_kernel_descriptor));
  desc->name = name;
  desc->block = block;
  desc->dim = dim;
  desc->device = 1;
  desc->index = 6;
  desc->hash = 5381;
  desc->hash = ((desc->hash << 5) + desc->hash) + 6;
  for ( int i=0; i<4; i++ ){
    desc->range[i] = range[i];
    desc->orig_range[i] = range[i];
    desc->hash = ((desc->hash << 5) + desc->hash) + range[i];
  }
  desc->nargs = 3;
  desc->args = (ops_arg*)malloc(3*sizeof(ops_arg));
  desc->args[0] = arg0;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg0.dat->index;
  desc->args[1] = arg1;
  desc->hash = ((desc->hash << 5) + desc->hash) + arg1.dat->index;
  desc->args[2] = arg2;
  desc->function = ops_par_loop_mgrid_restrict_kernel_execute;
  if (OPS_diags > 1) {
    ops_timing_realloc(6,"mgrid_restrict_kernel");
  }
  ops_enqueue_kernel(desc);
}
#endif