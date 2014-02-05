/* Crown Copyright 2012 AWE.

 This file is part of CloverLeaf.

 CloverLeaf is free software: you can redistribute it and/or modify it under
 the terms of the GNU General Public License as published by the
 Free Software Foundation, either version 3 of the License, or (at your option)
 any later version.

 CloverLeaf is distributed in the hope that it will be useful, but
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 details.

 You should have received a copy of the GNU General Public License along with
 CloverLeaf. If not, see http://www.gnu.org/licenses/. */

/** @brief Ideal gas kernel driver
 *  @author Wayne Gaudin, converted to OPS by Gihan Mudalige
 *  @details Invokes the user specified kernel for the ideal gas equation of
 *  state using the specified time level data.
**/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

// OPS header file

#include "ops_seq.h"
#include "ops_mpi_seq.h"

#include "data.h"
#include "definitions.h"

//Cloverleaf kernels
#include "ideal_gas_kernel.h"


void ideal_gas(int predict)
{
  //initialize sizes using global values
  int x_cells = grid->x_cells;
  int y_cells = grid->y_cells;
  int x_min = field->x_min;
  int x_max = field->x_max;
  int y_min = field->y_min;
  int y_max = field->y_max;

  int rangexy_inner[] = {x_min,x_max,y_min,y_max}; // inner range without border

  if(predict != TRUE) {
    ops_par_loop_mpi(ideal_gas_kernel, "ideal_gas_kernel", clover_grid, 2, rangexy_inner,
      ops_arg_dat(density0, S2D_00, "double", OPS_READ),
      ops_arg_dat(energy0, S2D_00, "double", OPS_READ),
      ops_arg_dat(pressure, S2D_00, "double", OPS_READ),
      ops_arg_dat(soundspeed, S2D_00, "double", OPS_READ));
  }
  else {
    ops_par_loop_mpi(ideal_gas_kernel, "ideal_gas_kernel", clover_grid, 2, rangexy_inner,
      ops_arg_dat(density1, S2D_00, "double", OPS_READ),
      ops_arg_dat(energy1, S2D_00, "double", OPS_READ),
      ops_arg_dat(pressure, S2D_00, "double", OPS_READ),
      ops_arg_dat(soundspeed, S2D_00, "double", OPS_READ));
  }
}
