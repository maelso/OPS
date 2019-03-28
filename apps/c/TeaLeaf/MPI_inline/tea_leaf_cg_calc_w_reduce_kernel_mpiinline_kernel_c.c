//
// auto-generated by ops.py
//

int xdim0_tea_leaf_cg_calc_w_reduce_kernel;
int xdim1_tea_leaf_cg_calc_w_reduce_kernel;
int xdim2_tea_leaf_cg_calc_w_reduce_kernel;
int xdim3_tea_leaf_cg_calc_w_reduce_kernel;


#define OPS_ACC0(x,y) (n_x*1+n_y*xdim0_tea_leaf_cg_calc_w_reduce_kernel*1+x+xdim0_tea_leaf_cg_calc_w_reduce_kernel*(y))
#define OPS_ACC1(x,y) (n_x*1+n_y*xdim1_tea_leaf_cg_calc_w_reduce_kernel*1+x+xdim1_tea_leaf_cg_calc_w_reduce_kernel*(y))
#define OPS_ACC2(x,y) (n_x*1+n_y*xdim2_tea_leaf_cg_calc_w_reduce_kernel*1+x+xdim2_tea_leaf_cg_calc_w_reduce_kernel*(y))
#define OPS_ACC3(x,y) (n_x*1+n_y*xdim3_tea_leaf_cg_calc_w_reduce_kernel*1+x+xdim3_tea_leaf_cg_calc_w_reduce_kernel*(y))

//user function



void tea_leaf_cg_calc_w_reduce_kernel_c_wrapper(
  double * restrict w,
  const double * restrict Kx,
  const double * restrict Ky,
  const double * restrict p,
  const double * restrict rx,
  const double * restrict ry,
  double * restrict pw_g,
  int x_size, int y_size) {
  double pw_v = *pw_g;
  #pragma omp parallel for reduction(+:pw_v)
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      double * restrict pw = &pw_v;
      
  w[OPS_ACC0(0,0)] = (1.0
                + (*ry)*(Ky[OPS_ACC2(0,1)] + Ky[OPS_ACC2(0,0)])
                + (*rx)*(Kx[OPS_ACC1(1,0)] + Kx[OPS_ACC1(0,0)]))*p[OPS_ACC3(0,0)]
                - (*ry)*(Ky[OPS_ACC2(0,1)]*p[OPS_ACC3(0,1)] + Ky[OPS_ACC2(0,0)]*p[OPS_ACC3(0,-1)])
                - (*rx)*(Kx[OPS_ACC1(1,0)]*p[OPS_ACC3(1,0)] + Kx[OPS_ACC1(0,0)]*p[OPS_ACC3(-1,0)]);
  *pw = *pw + w[OPS_ACC0(0,0)]*p[OPS_ACC3(0,0)];

    }
  }
  *pw_g = pw_v;
}
#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3

