//
// auto-generated by ops.py on 2013-11-04 16:04
//

#include "lib.h"
//user function
#include "calc_dt_kernel.h"

// host stub function
void ops_par_loop_calc_dt_kernel(char const *name, int dim, int* range,
 ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3,
 ops_arg arg4, ops_arg arg5, ops_arg arg6, ops_arg arg7,
 ops_arg arg8, ops_arg arg9, ops_arg arg10) {

  char *p_a[11];
  int  offs[11][2];
  int  count[dim];

  ops_arg args[11] = { arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9, arg10};


  for ( int i=0; i<11; i++ ){
    if (args[i].stencil!=NULL) {
      offs[i][0] = 1;  //unit step in x dimension
      offs[i][1] = ops_offs_set(range[0],range[2]+1, args[i]) - ops_offs_set(range[1],range[2], args[i]) +1;
      //stride in y as x stride is 0
      if (args[i].stencil->stride[0] == 0) {
        offs[i][0] = 0;
        offs[i][1] = args[i].dat->block_size[0];
      }
      //stride in x as y stride is 0
      else if (args[i].stencil->stride[1] == 0) {
        offs[i][0] = 1;
        offs[i][1] = -( range[1] - range[0] ) +1;
      }
    }
  }
  //set up initial pointers
  for ( int i=0; i<11; i++ ){
    if (args[i].argtype == OPS_ARG_DAT) {
      p_a[i] = (char *)args[i].data //base of 2D array
      +
      //y dimension -- get to the correct y line
      args[i].dat->size * args[i].dat->block_size[0] * ( range[2] * args[i].stencil->stride[1] - args[i].dat->offset[1] )
      +
      //x dimension - get to the correct x point on the y line
      args[i].dat->size * ( range[0] * args[i].stencil->stride[0] - args[i].dat->offset[0] );
    }
    else if (args[i].argtype == OPS_ARG_GBL) {
      p_a[i] = (char *)args[i].data;
    }
  }

  int total_range = 1;
  for ( int m=0; m<dim; m++ ){
    //number in each dimension
    count[m] = range[2*m+1]-range[2*m];
    total_range *= count[m];
  }
  //extra in last to ensure correct termination
  count[dim-1]++;


  xdim0 = args[0].dat->block_size[0];
  xdim1 = args[1].dat->block_size[0];
  xdim2 = args[2].dat->block_size[0];
  xdim3 = args[3].dat->block_size[0];
  xdim4 = args[4].dat->block_size[0];
  xdim5 = args[5].dat->block_size[0];
  xdim6 = args[6].dat->block_size[0];
  xdim7 = args[7].dat->block_size[0];
  xdim8 = args[8].dat->block_size[0];
  xdim9 = args[9].dat->block_size[0];
  xdim10 = args[10].dat->block_size[0];

  for ( int nt=0; nt<total_range; nt++ ){
    //call kernel function, passing in pointers to data

    calc_dt_kernel(  (double *)p_a[0], (double *)p_a[1], (double *)p_a[2],
           (double *)p_a[3], (double *)p_a[4], (double *)p_a[5], (double *)p_a[6],
           (double *)p_a[7], (double *)p_a[8], (double *)p_a[9], (double *)p_a[10] );

    //decrement counter
    count[0]--;

    //max dimension with changed index
    int m = 0;
    while ( (count[m]==0) ){
      count[m] = range[2*m+1]-range[2*m]; // reset counter
      m++;                                // next dimension
      count[m]--;                         // decrement counter
    }

    //shift pointers to data
    p_a[0]= p_a[0] + (args[0].dat->size * offs[0][m]);
    p_a[1]= p_a[1] + (args[1].dat->size * offs[1][m]);
    p_a[2]= p_a[2] + (args[2].dat->size * offs[2][m]);
    p_a[3]= p_a[3] + (args[3].dat->size * offs[3][m]);
    p_a[4]= p_a[4] + (args[4].dat->size * offs[4][m]);
    p_a[5]= p_a[5] + (args[5].dat->size * offs[5][m]);
    p_a[6]= p_a[6] + (args[6].dat->size * offs[6][m]);
    p_a[7]= p_a[7] + (args[7].dat->size * offs[7][m]);
    p_a[8]= p_a[8] + (args[8].dat->size * offs[8][m]);
    p_a[9]= p_a[9] + (args[9].dat->size * offs[9][m]);
    p_a[10]= p_a[10] + (args[10].dat->size * offs[10][m]);
  }

}