//
// auto-generated by ops.py
//

//user function
inline void poisson_kernel_initialguess(double *u) {
  u[OPS_ACC0(0,0)] = 0.0;
}





// host stub function
void ops_par_loop_poisson_kernel_initialguess(char const *name, ops_block block, int dim, int* range,
 ops_arg arg0) {

  //Timing
  double t1,t2,c1,c2;

  char *p_a[1];
  int  offs[1][2];
  ops_arg args[1] = { arg0};



  #ifdef CHECKPOINTING
  if (!ops_checkpointing_before(args,1,range,2)) return;
  #endif

  if (OPS_diags > 1) {
    ops_timing_realloc(2,"poisson_kernel_initialguess");
    OPS_kernels[2].count++;
    ops_timers_core(&c2,&t2);
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
  #else
  for ( int n=0; n<2; n++ ){
    start[n] = range[2*n];end[n] = range[2*n+1];
  }
  #endif
  #ifdef OPS_DEBUG
  ops_register_args(args, "poisson_kernel_initialguess");
  #endif

  offs[0][0] = args[0].stencil->stride[0]*1;  //unit step in x dimension
  offs[0][1] = off2D(1, &start[0],
      &end[0],args[0].dat->size, args[0].stencil->stride) - offs[0][0];



  int off0_0 = offs[0][0];
  int off0_1 = offs[0][1];
  int dat0 = args[0].dat->elem_size;

  //set up initial pointers and exchange halos if necessary
  int base0 = args[0].dat->base_offset + args[0].dat->elem_size * start[0] * args[0].stencil->stride[0];
  base0 = base0+ args[0].dat->elem_size *
    args[0].dat->size[0] *
    start[1] * args[0].stencil->stride[1];
  p_a[0] = (char *)args[0].data + base0;


  //initialize global variable with the dimension of dats
  xdim0 = args[0].dat->size[0];

  //Halo Exchanges
  ops_H_D_exchanges_host(args, 1);
  ops_halo_exchanges(args,1,range);
  ops_H_D_exchanges_host(args, 1);

  if (OPS_diags > 1) {
    ops_timers_core(&c1,&t1);
    OPS_kernels[2].mpi_time += t1-t2;
  }

  int n_x;
  for ( int n_y=start[1]; n_y<end[1]; n_y++ ){
    #pragma novector
    for( n_x=start[0]; n_x<start[0]+((end[0]-start[0])/SIMD_VEC)*SIMD_VEC; n_x+=SIMD_VEC ) {
      //call kernel function, passing in pointers to data -vectorised
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        poisson_kernel_initialguess(  (double *)p_a[0]+ i*1*1 );

      }

      //shift pointers to data x direction
      p_a[0]= p_a[0] + (dat0 * off0_0)*SIMD_VEC;
    }

    for ( int n_x=start[0]+((end[0]-start[0])/SIMD_VEC)*SIMD_VEC; n_x<end[0]; n_x++ ){
      //call kernel function, passing in pointers to data - remainder
      poisson_kernel_initialguess(  (double *)p_a[0] );


      //shift pointers to data x direction
      p_a[0]= p_a[0] + (dat0 * off0_0);
    }

    //shift pointers to data y direction
    p_a[0]= p_a[0] + (dat0 * off0_1);
  }
  if (OPS_diags > 1) {
    ops_timers_core(&c2,&t2);
    OPS_kernels[2].time += t2-t1;
  }
  ops_set_dirtybit_host(args, 1);
  ops_set_halo_dirtybit3(&args[0],range);

  if (OPS_diags > 1) {
    //Update kernel record
    ops_timers_core(&c1,&t1);
    OPS_kernels[2].mpi_time += t1-t2;
    OPS_kernels[2].transfer += ops_compute_transfer(dim, start, end, &arg0);
  }
}
