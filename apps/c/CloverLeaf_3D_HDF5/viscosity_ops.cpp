//
// auto-generated by ops.py
//

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define OPS_3D
#include "ops_lib_cpp.h"

//
// ops_par_loop declarations
//

void ops_par_loop_viscosity_kernel(char const *, ops_block, int, int *, ops_arg,
                                   ops_arg, ops_arg, ops_arg, ops_arg, ops_arg,
                                   ops_arg, ops_arg, ops_arg, ops_arg, ops_arg,
                                   ops_arg);

#include "data.h"
#include "definitions.h"

//#include "viscosity_kernel.h"

void viscosity_func() {

  int x_min = field.x_min;
  int x_max = field.x_max;
  int y_min = field.y_min;
  int y_max = field.y_max;
  int z_min = field.z_min;
  int z_max = field.z_max;

  int rangexyz_inner[] = {x_min, x_max, y_min, y_max, z_min, z_max};

  ops_par_loop_viscosity_kernel(
      "viscosity_kernel", clover_grid, 3, rangexyz_inner,
      ops_arg_dat(xvel0, 1, S3D_000_fP1P1P1, "double", OPS_READ),
      ops_arg_dat(yvel0, 1, S3D_000_fP1P1P1, "double", OPS_READ),
      ops_arg_dat(celldx, 1, S3D_000_P100_STRID3D_X, "double", OPS_READ),
      ops_arg_dat(celldy, 1, S3D_000_0P10_STRID3D_Y, "double", OPS_READ),
      ops_arg_dat(pressure, 1, S3D_P100_M100_0P10_0M10_00P1_00M1, "double",
                  OPS_READ),
      ops_arg_dat(density0, 1, S3D_000, "double", OPS_READ),
      ops_arg_dat(viscosity, 1, S3D_000, "double", OPS_WRITE),
      ops_arg_dat(zvel0, 1, S3D_000_fP1P1P1, "double", OPS_READ),
      ops_arg_dat(celldz, 1, S3D_000_00P1_STRID3D_Z, "double", OPS_READ),
      ops_arg_dat(xarea, 1, S3D_000, "double", OPS_READ),
      ops_arg_dat(yarea, 1, S3D_000, "double", OPS_READ),
      ops_arg_dat(zarea, 1, S3D_000, "double", OPS_READ));
}
