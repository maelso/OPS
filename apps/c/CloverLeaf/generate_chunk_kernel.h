#ifndef GENERATE_CHUNK_KERNEL_H
#define GENERATE_CHUNK_KERNEL_H

#include "data.h"
#include "definitions.h"



void generate_chunk_kernel( const double *vertexx, const double *vertexy,
                     double *energy0, double *density0,
                     double *xvel0,  double *yvel0,
                     const double *cellx, const double *celly) {

  double radius, x_cent, y_cent;
  int is_in = 0;
  int is_in2 = 0;
  //State 1 is always the background state

  energy0[OPS_ACC2(0,0)]= states[0].energy;
  density0[OPS_ACC3(0,0)]= states[0].density;
  xvel0[OPS_ACC4(0,0)]=states[0].xvel;
  yvel0[OPS_ACC5(0,0)]=states[0].yvel;

  for(int i = 1; i<number_of_states; i++) {

    x_cent=states[i].xmin;
    y_cent=states[i].ymin;
    is_in = 0;
    is_in2 = 0;

    if (states[i].geometry == g_rect) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          if(vertexx[OPS_ACC0(1+i1,0)] >= states[i].xmin  && vertexx[OPS_ACC0(0+i1,0)] < states[i].xmax) {
            if(vertexy[OPS_ACC1(0,1+j1)] >= states[i].ymin && vertexy[OPS_ACC1(0,0+j1)] < states[i].ymax) {
              is_in = 1;
            }
          }
        }
      }
      if(vertexx[OPS_ACC0(1,0)] >= states[i].xmin  && vertexx[OPS_ACC0(0,0)] < states[i].xmax) {
        if(vertexy[OPS_ACC1(0,1)] >= states[i].ymin && vertexy[OPS_ACC1(0,0)] < states[i].ymax) {
          is_in2 = 1;
        }
      }
      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      }
      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
    else if(states[i].geometry == g_circ) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          radius = sqrt ((cellx[OPS_ACC6(i1,0)] - x_cent) * (cellx[OPS_ACC6(i1,0)] - x_cent) +
                     (celly[OPS_ACC7(0,j1)] - y_cent) * (celly[OPS_ACC7(0,j1)] - y_cent));
          if (radius <= states[i].radius) {
            is_in = 1;
          }
        }
      }
      if (radius <= states[i].radius) is_in2 = 1;

      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      }

      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
    else if(states[i].geometry == g_point) {
      for (int i1 = -1; i1 <= 0; i1++) {
        for (int j1 = -1; j1 <= 0; j1++) {
          if(vertexx[OPS_ACC0(i1,0)] == x_cent && vertexy[OPS_ACC1(0,j1)] == y_cent) {
            is_in = 1;
          }
        }
      }
      if(vertexx[OPS_ACC0(0,0)] == x_cent && vertexy[OPS_ACC1(0,0)] == y_cent) 
        is_in2 = 1;

      if (is_in2) {
        energy0[OPS_ACC2(0,0)] = states[i].energy;
        density0[OPS_ACC3(0,0)] = states[i].density;
      } 

      if (is_in) {
        xvel0[OPS_ACC4(0,0)] = states[i].xvel;
        yvel0[OPS_ACC5(0,0)] = states[i].yvel;
      }
    }
  }
}


#endif
