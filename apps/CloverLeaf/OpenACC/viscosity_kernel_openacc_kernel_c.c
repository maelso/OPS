//
// auto-generated by ops.py on 2014-07-11 14:02
//

#include "./OpenACC/clover_leaf_common.h"

#define OPS_GPU

int xdim0_viscosity_kernel;
int xdim1_viscosity_kernel;
int xdim2_viscosity_kernel;
int xdim3_viscosity_kernel;
int xdim4_viscosity_kernel;
int xdim5_viscosity_kernel;
int xdim6_viscosity_kernel;

#define OPS_ACC0(x,y) (x+xdim0_viscosity_kernel*(y))
#define OPS_ACC1(x,y) (x+xdim1_viscosity_kernel*(y))
#define OPS_ACC2(x,y) (x+xdim2_viscosity_kernel*(y))
#define OPS_ACC3(x,y) (x+xdim3_viscosity_kernel*(y))
#define OPS_ACC4(x,y) (x+xdim4_viscosity_kernel*(y))
#define OPS_ACC5(x,y) (x+xdim5_viscosity_kernel*(y))
#define OPS_ACC6(x,y) (x+xdim6_viscosity_kernel*(y))

//user function
inline 
void viscosity_kernel( const double *xvel0, const double *yvel0,
                       const double *celldx, const double *celldy,
                       const double *pressure, const double *density0,
                       double *viscosity) {

  double ugrad, vgrad,
         grad2,
         pgradx,pgrady,
         pgradx2,pgrady2,
         grad,
         ygrad, xgrad,
         div,
         strain2,
         limiter,
         pgrad;


  ugrad = (xvel0[OPS_ACC0(1,0)] + xvel0[OPS_ACC0(1,1)]) - (xvel0[OPS_ACC0(0,0)] + xvel0[OPS_ACC0(0,1)]);
  vgrad = (yvel0[OPS_ACC1(0,1)] + yvel0[OPS_ACC1(1,1)]) - (yvel0[OPS_ACC1(0,0)] + yvel0[OPS_ACC1(1,0)]);

  div = (celldx[OPS_ACC2(0,0)])*(ugrad) + (celldy[OPS_ACC3(0,0)])*(vgrad);

  strain2 = 0.5*(xvel0[OPS_ACC0(0,1)] + xvel0[OPS_ACC0(1,1)] - xvel0[OPS_ACC0(0,0)] - xvel0[OPS_ACC0(1,0)])/(celldy[OPS_ACC3(0,0)]) +
            0.5*(yvel0[OPS_ACC1(1,0)] + yvel0[OPS_ACC1(1,1)] - yvel0[OPS_ACC1(0,0)] - yvel0[OPS_ACC1(0,1)])/(celldx[OPS_ACC2(0,0)]);


  pgradx  = (pressure[OPS_ACC4(1,0)] - pressure[OPS_ACC4(-1,0)])/(celldx[OPS_ACC2(0,0)]+ celldx[OPS_ACC2(1,0)]);
  pgrady = (pressure[OPS_ACC4(0,1)] - pressure[OPS_ACC4(0,-1)])/(celldy[OPS_ACC3(0,0)]+ celldy[OPS_ACC3(0,1)]);

  pgradx2 = pgradx * pgradx;
  pgrady2 = pgrady * pgrady;

  limiter = ((0.5*(ugrad)/celldx[OPS_ACC2(0,0)]) * pgradx2 +
             (0.5*(vgrad)/celldy[OPS_ACC3(0,0)]) * pgrady2 +
              strain2 * pgradx * pgrady)/ MAX(pgradx2 + pgrady2 , 1.0e-16);

  if( (limiter > 0.0) || (div >= 0.0)) {
        viscosity[OPS_ACC6(0,0)] = 0.0;
  }
  else {
    pgradx = SIGN( MAX(1.0e-16, fabs(pgradx)), pgradx);
    pgrady = SIGN( MAX(1.0e-16, fabs(pgrady)), pgrady);
    pgrad = sqrt(pgradx*pgradx + pgrady*pgrady);
    xgrad = fabs(celldx[OPS_ACC2(0,0)] * pgrad/pgradx);
    ygrad = fabs(celldy[OPS_ACC3(0,0)] * pgrad/pgrady);
    grad  = MIN(xgrad,ygrad);
    grad2 = grad*grad;

    viscosity[OPS_ACC6(0,0)] = 2.0 * (density0[OPS_ACC5(0,0)]) * grad2 * limiter * limiter;
  }
}



#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5
#undef OPS_ACC6


void viscosity_kernel_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  double *p_a4,
  double *p_a5,
  double *p_a6,
  int x_size, int y_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1,p_a2,p_a3,p_a4,p_a5,p_a6)
  #pragma acc loop
  #endif
  for ( int n_y=0; n_y<y_size; n_y++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_x=0; n_x<x_size; n_x++ ){
      viscosity_kernel(  p_a0 + n_x*1 + n_y*xdim0_viscosity_kernel*1,
           p_a1 + n_x*1 + n_y*xdim1_viscosity_kernel*1, p_a2 + n_x*1 + n_y*xdim2_viscosity_kernel*0,
           p_a3 + n_x*0 + n_y*xdim3_viscosity_kernel*1, p_a4 + n_x*1 + n_y*xdim4_viscosity_kernel*1,
           p_a5 + n_x*1 + n_y*xdim5_viscosity_kernel*1, p_a6 + n_x*1 + n_y*xdim6_viscosity_kernel*1 );

    }
  }
}
