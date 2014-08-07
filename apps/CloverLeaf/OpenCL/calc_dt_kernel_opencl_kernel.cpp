//
// auto-generated by ops.py on 2014-08-07 15:58
//


#ifdef OCL_FMA_SWITCH_ON
#define OCL_FMA 1
#else
#define OCL_FMA 0
#endif


static bool isbuilt_calc_dt_kernel = false;

void buildOpenCLKernels_calc_dt_kernel(int xdim0,
 int xdim1,
 int xdim2,
 int xdim3,
 int xdim4,
 int xdim5,
 int xdim6,
 int xdim7,
 int xdim8,
 int xdim9,
 int xdim10) {

  //int ocl_fma = OCL_FMA;
  if(!isbuilt_calc_dt_kernel) {
    buildOpenCLKernels();
    //clSafeCall( clUnloadCompiler() );
    cl_int ret;
    char* source_filename[1] = {"./OpenCL/calc_dt_kernel.cl"};

    // Load the kernel source code into the array source_str
    FILE *fid;
    char *source_str[1];
    size_t source_size[1];

    for(int i=0; i<1; i++) {
      fid = fopen(source_filename[i], "r");
      if (!fid) {
        fprintf(stderr, "Can't open the kernel source file!\n");
        exit(1);
      }

      source_str[i] = (char*)malloc(4*0x1000000);
      source_size[i] = fread(source_str[i], 1, 4*0x1000000, fid);
      if(source_size[i] != 4*0x1000000) {
        if (ferror(fid)) {
          printf ("Error while reading kernel source file %s\n", source_filename[i]);
          exit(-1);
        }
        if (feof(fid))
          printf ("Kernel source file %s succesfuly read.\n", source_filename[i]);
          //printf("%s\n",source_str[i]);
      }
      fclose(fid);
    }
    
    printf("Compiling calc_dt_kernel %d source -- start \n",OCL_FMA);

      // Create a program from the source
      OPS_opencl_core.program = clCreateProgramWithSource(OPS_opencl_core.context, 1, (const char **) &source_str, (const size_t *) &source_size, &ret);
      clSafeCall( ret );

      // Build the program
      char buildOpts[255*11];
      char* pPath = NULL;
      pPath = getenv ("OPS_INSTALL_PATH");
      if (pPath!=NULL)
        if(OCL_FMA)
          sprintf(buildOpts,"-cl-mad-enable -DOCL_FMA -I%s/include -DOPS_WARPSIZE=%d  -Dxdim0_calc_dt_kernel=%d -Dxdim1_calc_dt_kernel=%d -Dxdim2_calc_dt_kernel=%d -Dxdim3_calc_dt_kernel=%d -Dxdim4_calc_dt_kernel=%d -Dxdim5_calc_dt_kernel=%d -Dxdim6_calc_dt_kernel=%d -Dxdim7_calc_dt_kernel=%d -Dxdim8_calc_dt_kernel=%d -Dxdim9_calc_dt_kernel=%d -Dxdim10_calc_dt_kernel=%d", pPath, 32,xdim0,xdim1,xdim2,xdim3,xdim4,xdim5,xdim6,xdim7,xdim8,xdim9,xdim10);
        else
          sprintf(buildOpts,"-cl-mad-enable -I%s/include -DOPS_WARPSIZE=%d  -Dxdim0_calc_dt_kernel=%d -Dxdim1_calc_dt_kernel=%d -Dxdim2_calc_dt_kernel=%d -Dxdim3_calc_dt_kernel=%d -Dxdim4_calc_dt_kernel=%d -Dxdim5_calc_dt_kernel=%d -Dxdim6_calc_dt_kernel=%d -Dxdim7_calc_dt_kernel=%d -Dxdim8_calc_dt_kernel=%d -Dxdim9_calc_dt_kernel=%d -Dxdim10_calc_dt_kernel=%d", pPath, 32,xdim0,xdim1,xdim2,xdim3,xdim4,xdim5,xdim6,xdim7,xdim8,xdim9,xdim10);
      else {
        sprintf("Incorrect OPS_INSTALL_PATH %s\n",pPath);
        exit(EXIT_FAILURE);
      }

      ret = clBuildProgram(OPS_opencl_core.program, 1, &OPS_opencl_core.device_id, buildOpts, NULL, NULL);

      if(ret != CL_SUCCESS) {
        char* build_log;
        size_t log_size;
        clSafeCall( clGetProgramBuildInfo(OPS_opencl_core.program, OPS_opencl_core.device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &log_size) );
        build_log = (char*) malloc(log_size+1);
        clSafeCall( clGetProgramBuildInfo(OPS_opencl_core.program, OPS_opencl_core.device_id, CL_PROGRAM_BUILD_LOG, log_size, build_log, NULL) );
        build_log[log_size] = '\0';
        fprintf(stderr, "=============== OpenCL Program Build Info ================\n\n%s", build_log);
        fprintf(stderr, "\n========================================================= \n");
        free(build_log);
        exit(EXIT_FAILURE);
      }
      printf("compiling calc_dt_kernel -- done\n");

    // Create the OpenCL kernel
    OPS_opencl_core.kernel[71] = clCreateKernel(OPS_opencl_core.program, "ops_calc_dt_kernel", &ret);
    clSafeCall( ret );

    isbuilt_calc_dt_kernel = true;
  }

}


// host stub function
void ops_par_loop_calc_dt_kernel(char const *name, ops_block Block, int dim, int* range,
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

  //build opencl kernel if not already built

  buildOpenCLKernels_calc_dt_kernel(
  xdim0,
  xdim1,
  xdim2,
  xdim3,
  xdim4,
  xdim5,
  xdim6,
  xdim7,
  xdim8,
  xdim9,
  xdim10);

  //Timing
  double t1,t2,c1,c2;
  ops_timing_realloc(71,"calc_dt_kernel");
  ops_timers_core(&c2,&t2);

  //set up OpenCL thread blocks
  size_t globalWorkSize[3] = {((x_size-1)/OPS_block_size_x+ 1)*OPS_block_size_x, ((y_size-1)/OPS_block_size_y + 1)*OPS_block_size_y, 1};
  size_t localWorkSize[3] =  {OPS_block_size_x,OPS_block_size_y,1};




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

  //set up initial pointers
  int base0 = 1 * 
  (start[0] * args[0].stencil->stride[0] - args[0].dat->base[0] - args[0].dat->d_m[0]);
  base0 = base0 + args[0].dat->size[0] *
  (start[1] * args[0].stencil->stride[1] - args[0].dat->base[1] - args[0].dat->d_m[1]);

  int base1 = 1 * 
  (start[0] * args[1].stencil->stride[0] - args[1].dat->base[0] - args[1].dat->d_m[0]);
  base1 = base1 + args[1].dat->size[0] *
  (start[1] * args[1].stencil->stride[1] - args[1].dat->base[1] - args[1].dat->d_m[1]);

  int base2 = 1 * 
  (start[0] * args[2].stencil->stride[0] - args[2].dat->base[0] - args[2].dat->d_m[0]);
  base2 = base2 + args[2].dat->size[0] *
  (start[1] * args[2].stencil->stride[1] - args[2].dat->base[1] - args[2].dat->d_m[1]);

  int base3 = 1 * 
  (start[0] * args[3].stencil->stride[0] - args[3].dat->base[0] - args[3].dat->d_m[0]);
  base3 = base3 + args[3].dat->size[0] *
  (start[1] * args[3].stencil->stride[1] - args[3].dat->base[1] - args[3].dat->d_m[1]);

  int base4 = 1 * 
  (start[0] * args[4].stencil->stride[0] - args[4].dat->base[0] - args[4].dat->d_m[0]);
  base4 = base4 + args[4].dat->size[0] *
  (start[1] * args[4].stencil->stride[1] - args[4].dat->base[1] - args[4].dat->d_m[1]);

  int base5 = 1 * 
  (start[0] * args[5].stencil->stride[0] - args[5].dat->base[0] - args[5].dat->d_m[0]);
  base5 = base5 + args[5].dat->size[0] *
  (start[1] * args[5].stencil->stride[1] - args[5].dat->base[1] - args[5].dat->d_m[1]);

  int base6 = 1 * 
  (start[0] * args[6].stencil->stride[0] - args[6].dat->base[0] - args[6].dat->d_m[0]);
  base6 = base6 + args[6].dat->size[0] *
  (start[1] * args[6].stencil->stride[1] - args[6].dat->base[1] - args[6].dat->d_m[1]);

  int base7 = 1 * 
  (start[0] * args[7].stencil->stride[0] - args[7].dat->base[0] - args[7].dat->d_m[0]);
  base7 = base7 + args[7].dat->size[0] *
  (start[1] * args[7].stencil->stride[1] - args[7].dat->base[1] - args[7].dat->d_m[1]);

  int base8 = 1 * 
  (start[0] * args[8].stencil->stride[0] - args[8].dat->base[0] - args[8].dat->d_m[0]);
  base8 = base8 + args[8].dat->size[0] *
  (start[1] * args[8].stencil->stride[1] - args[8].dat->base[1] - args[8].dat->d_m[1]);

  int base9 = 1 * 
  (start[0] * args[9].stencil->stride[0] - args[9].dat->base[0] - args[9].dat->d_m[0]);
  base9 = base9 + args[9].dat->size[0] *
  (start[1] * args[9].stencil->stride[1] - args[9].dat->base[1] - args[9].dat->d_m[1]);

  int base10 = 1 * 
  (start[0] * args[10].stencil->stride[0] - args[10].dat->base[0] - args[10].dat->d_m[0]);
  base10 = base10 + args[10].dat->size[0] *
  (start[1] * args[10].stencil->stride[1] - args[10].dat->base[1] - args[10].dat->d_m[1]);


  ops_H_D_exchanges_device(args, 11);

  ops_timers_core(&c1,&t1);
  OPS_kernels[71].mpi_time += t1-t2;


  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 0, sizeof(cl_mem), (void*) &arg0.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 1, sizeof(cl_mem), (void*) &arg1.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 2, sizeof(cl_mem), (void*) &arg2.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 3, sizeof(cl_mem), (void*) &arg3.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 4, sizeof(cl_mem), (void*) &arg4.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 5, sizeof(cl_mem), (void*) &arg5.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 6, sizeof(cl_mem), (void*) &arg6.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 7, sizeof(cl_mem), (void*) &arg7.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 8, sizeof(cl_mem), (void*) &arg8.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 9, sizeof(cl_mem), (void*) &arg9.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 10, sizeof(cl_mem), (void*) &arg10.data_d ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 11, sizeof(cl_double), (void*) &g_small ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 12, sizeof(cl_double), (void*) &g_big ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 13, sizeof(cl_double), (void*) &dtc_safe ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 14, sizeof(cl_double), (void*) &dtu_safe ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 15, sizeof(cl_double), (void*) &dtv_safe ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 16, sizeof(cl_double), (void*) &dtdiv_safe ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 17, sizeof(cl_int), (void*) &base0 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 18, sizeof(cl_int), (void*) &base1 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 19, sizeof(cl_int), (void*) &base2 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 20, sizeof(cl_int), (void*) &base3 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 21, sizeof(cl_int), (void*) &base4 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 22, sizeof(cl_int), (void*) &base5 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 23, sizeof(cl_int), (void*) &base6 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 24, sizeof(cl_int), (void*) &base7 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 25, sizeof(cl_int), (void*) &base8 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 26, sizeof(cl_int), (void*) &base9 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 27, sizeof(cl_int), (void*) &base10 ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 28, sizeof(cl_int), (void*) &x_size ));
  clSafeCall( clSetKernelArg(OPS_opencl_core.kernel[71], 29, sizeof(cl_int), (void*) &y_size ));

  //call/enque opencl kernel wrapper function
  clSafeCall( clEnqueueNDRangeKernel(OPS_opencl_core.command_queue, OPS_opencl_core.kernel[71], 3, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL) );
  if (OPS_diags>1) {
    clSafeCall( clFinish(OPS_opencl_core.command_queue) );
  }

  ops_set_dirtybit_device(args, 11);

  //Update kernel record
  ops_timers_core(&c2,&t2);
  OPS_kernels[71].count++;
  OPS_kernels[71].time += t2-t1;
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg0);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg1);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg2);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg3);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg4);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg5);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg6);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg7);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg8);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg9);
  OPS_kernels[71].transfer += ops_compute_transfer(dim, range, &arg10);
}
