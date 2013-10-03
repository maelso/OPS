/*
* Open source copyright declaration based on BSD open source template:
* http://www.opensource.org/licenses/bsd-license.php
*
* This file is part of the OPS distribution.
*
* Copyright (c) 2013, Mike Giles and others. Please see the AUTHORS file in
* the main source directory for a full list of copyright holders.
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* * Redistributions of source code must retain the above copyright
* notice, this list of conditions and the following disclaimer.
* * Redistributions in binary form must reproduce the above copyright
* notice, this list of conditions and the following disclaimer in the
* documentation and/or other materials provided with the distribution.
* * The name of Mike Giles may not be used to endorse or promote products
* derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY Mike Giles ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL Mike Giles BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
* ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
* (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/** @brief headder file declaring the functions for the ops sequential backend
  * @author Gihan Mudalige, Istvan Reguly
  * @details Declares the OPS API calls for the sequential backend
  */

#include "ops_lib_cpp.h"

inline void ops_arg_set(int n_x, ops_arg arg, char **p_arg){
  if (arg.stencil!=NULL) {

    for (int i = 0; i < arg.stencil->points; i++){
      p_arg[i] =
         arg.data + //base of 1D array
         arg.dat->size * ( //multiply by the number of bytes per element
         (n_x - arg.dat->offset[0]) * //calculate the offset from index 0
         arg.stencil->stride[0]  + // jump in strides ??
         arg.stencil->stencil[i * arg.stencil->dims + 0] //get the value at the ith stencil point
         );
    }
  } else {
    *p_arg = arg.data;
  }
}

inline void ops_arg_set(int n_x,
                        int n_y, ops_arg arg, char **p_arg){
  if (arg.stencil!=NULL) {
    for (int i = 0; i < arg.stencil->points; i++){
      p_arg[i] =
        arg.data + //base of 2D array
        //y dimension -- get to the correct y line
        arg.dat->size * arg.dat->block_size[0] * ( //multiply by the number of
                                                    //bytes per element and xdim block size
        (n_y - arg.dat->offset[1]) * // calculate the offset from index 0 for y dim
        arg.stencil->stride[1] + // jump in strides in y dim ??
        arg.stencil->stencil[i*arg.stencil->dims + 1]) //get the value at the ith
                                                       //stencil point "+ 1" is the y dim
        +
        //x dimension - get to the correct x point on the y line
        arg.dat->size * ( //multiply by the number of bytes per element
        (n_x - arg.dat->offset[0]) * //calculate the offset from index 0 for x dim
        arg.stencil->stride[0] + // jump in strides in x dim ??
        arg.stencil->stencil[i*arg.stencil->dims + 0] //get the value at the ith
                                                      //stencil point "+ 0" is the x dim
      );
    }
  } else {
    *p_arg = arg.data;
  }
}

inline void ops_arg_set(int n_x,
                        int n_y,
                        int n_z, ops_arg arg, char **p_arg){
  if (arg.stencil!=NULL) {
    for (int i = 0; i < arg.stencil->points; i++)
      p_arg[i] =
      arg.data +
      //z dimension - get to the correct z plane
      arg.dat->size * arg.dat->block_size[1] * arg.dat->block_size[0] * (
      (n_z - arg.dat->offset[2]) *
      arg.stencil->stride[2] +
      arg.stencil->stencil[i*arg.stencil->dims + 2])
      +
      //y dimension -- get to the correct y line on the z plane
      arg.dat->size * arg.dat->block_size[0] * (
      (n_y - arg.dat->offset[1]) *
      arg.stencil->stride[1] +
      arg.stencil->stencil[i*arg.stencil->dims + 1])
      +
      //x dimension - get to the correct x point on the y line
      arg.dat->size * (
      (n_x - arg.dat->offset[0]) *
      arg.stencil->stride[0] +
      arg.stencil->stencil[i*arg.stencil->dims + 0]
      );
  } else {
    *p_arg = arg.data;
  }
}



inline void ops_args_set(int iter_x,
                         int nargs, ops_arg *args, char ***p_a){
  for (int n=0; n<nargs; n++) {
    ops_arg_set(iter_x, args[n], p_a[n]);
  }
}

inline void ops_args_set(int iter_x,
                         int iter_y,
                         int nargs, ops_arg *args, char ***p_a){
  for (int n=0; n<nargs; n++) {
    ops_arg_set(iter_x, iter_y, args[n], p_a[n]);
  }
}

inline void ops_args_set(int iter_x,
                         int iter_y,
                         int iter_z, int nargs, ops_arg *args, char ***p_a){
  for (int n=0; n<nargs; n++) {
    ops_arg_set(iter_x, iter_y, iter_z, args[n], p_a[n]);
  }
}



template < class T0 >
void ops_par_loop_opt(void (*kernel)( T0*),
  char const * name, int dim, int *range,
  ops_arg arg0 ) {

  char  **p_a[1];
  int   *offs[1] = {arg0.stencil->stencil};
  int   count[dim];

  ops_arg args[1] = {arg0};



  if (args[0].argtype == OPS_ARG_DAT){
    p_a[0] = (char **)malloc(args[0].stencil->points * sizeof(char *));
  }
  else if (args[0].argtype == OPS_ARG_GBL)
      p_a[0] = (char **)malloc(args[0].dim * sizeof(char *));

  int total_range = 1;
  for (int m=0; m<dim; m++) {
    count[m] = range[2*m+1]-range[2*m];  // number in each dimension
    total_range *= count[m];
  }
  count[dim-1]++;     // extra in last to ensure correct termination


  // loop over set elements

  //p_a[0][0] = args[0].data;
  ops_args_set(range[0], range[2],1,args,p_a);

  for (int nt=0; nt<total_range; nt++) {


    // call kernel function, passing in pointers to data
    kernel( (T0 *)p_a[0] );

    count[0]--;   // decrement counter
    int m = 0;    // max dimension with changed index

    while (count[m]==0) {
      count[m] = range[2*m+1]-range[2*m]; // reset counter
      m++;                                // next dimension
      count[m]--;                         // decrement counter
    }

    // shift pointers to data
    for (int np=0; np<args[0].stencil->points; np++) {
      //p_a[0][np] += args[0].dat->size * m;
      p_a[0][np] += args[0].dat->size * offs[0][m];
    }
  }

  for (int i = 0; i < 1; i++) {
  if (args[i].argtype == OPS_ARG_DAT)
    free(p_a[i]);
  }
}


template < class T0, class T1 >
void ops_par_loop_opt(void (*kernel)( T0*, T1* ),
                  char const * name, int dim, int *range,
                  ops_arg arg0, ops_arg arg1 ) {

  char  **p_a[2];
  int   *offs[2] = {arg0.stencil->stencil, arg1.stencil->stencil};
  int   count[dim];

  ops_arg args[2] = {arg0, arg1};

  for (int i = 0; i < 2; i++) {
    if (args[i].argtype == OPS_ARG_DAT)
      p_a[i] = (char **)malloc(args[i].stencil->points * sizeof(char *));
    else if (args[i].argtype == OPS_ARG_GBL)
      p_a[i] = (char **)malloc(args[i].dim * sizeof(char *));
  }

  int total_range = 1;
  for (int m=0; m<dim; m++) {
    count[m] = range[2*m+1]-range[2*m];  // number in each dimension
    total_range *= count[m];
  }
  count[dim-1]++;     // extra in last to ensure correct termination

  ops_args_set(range[0], range[2],1,args,p_a);

  for (int nt=0; nt<total_range; nt++) {

    // call kernel function, passing in pointers to data
    kernel( (T0 *)p_a[0], (T0 *)p_a[1] );

    count[0]--;   // decrement counter
    int m = 0;    // max dimension with changed index

    while (count[m]==0) {
      count[m] = range[2*m+1]-range[2*m]; // reset counter
      m++;                                // next dimension
      count[m]--;                         // decrement counter
    }

    // shift pointers to data
    for (int i=0; i<2; i++) {
      for (int np=0; np<args[i].stencil->points; np++) {
        p_a[i][np] = p_a[i][np] + (args[i].dat->size * offs[i][m]);
      }
    }
  }

  for (int i = 0; i < 2; i++) {
  if (args[i].argtype == OPS_ARG_DAT)
    free(p_a[i]);
  }
}



template < class T0, class T1, class T2, class T3  >
void ops_par_loop_opt(void (*kernel)( T0*, T1*, T2*, T3*),
                  char const * name, int dim, int *range,
                  ops_arg arg0, ops_arg arg1, ops_arg arg2, ops_arg arg3  ) {

  char  **p_a[4];
  int   *offs[4] = {arg0.stencil->stencil, arg1.stencil->stencil,
                    arg2.stencil->stencil, arg3.stencil->stencil};
  int   count[dim];

  ops_arg args[4] = {arg0, arg1, arg2, arg3};

  for (int i = 0; i < 4; i++) {
    if (args[i].argtype == OPS_ARG_DAT)
      p_a[i] = (char **)malloc(args[i].stencil->points * sizeof(char *));
    else if (args[i].argtype == OPS_ARG_GBL)
      p_a[i] = (char **)malloc(args[i].dim * sizeof(char *));
  }

  int total_range = 1;
  for (int m=0; m<dim; m++) {
    count[m] = range[2*m+1]-range[2*m];  // number in each dimension
    total_range *= count[m];
  }
  count[dim-1]++;     // extra in last to ensure correct termination

  ops_args_set(range[0], range[2],1,args,p_a);

  for (int nt=0; nt<total_range; nt++) {

    // call kernel function, passing in pointers to data
    kernel( (T0 *)p_a[0], (T0 *)p_a[1], (T0 *)p_a[2], (T0 *)p_a[3] );

    count[0]--;   // decrement counter
    int m = 0;    // max dimension with changed index

    while (count[m]==0) {
      count[m] = range[2*m+1]-range[2*m]; // reset counter
      m++;                                // next dimension
      count[m]--;                         // decrement counter
    }

    // shift pointers to data
    for (int i=0; i<4; i++) {
      for (int np=0; np<args[i].stencil->points; np++) {
        p_a[i][np] += args[i].dat->size * offs[i][m];
      }
    }
  }

  for (int i = 0; i < 4; i++) {
  if (args[i].argtype == OPS_ARG_DAT)
    free(p_a[i]);
  }
}