//
// auto-generated by ops.py
//

#define OPS_GPU

int xdim0_advec_cell_kernel1_xdir;
int ydim0_advec_cell_kernel1_xdir;
int xdim1_advec_cell_kernel1_xdir;
int ydim1_advec_cell_kernel1_xdir;
int xdim2_advec_cell_kernel1_xdir;
int ydim2_advec_cell_kernel1_xdir;
int xdim3_advec_cell_kernel1_xdir;
int ydim3_advec_cell_kernel1_xdir;
int xdim4_advec_cell_kernel1_xdir;
int ydim4_advec_cell_kernel1_xdir;
int xdim5_advec_cell_kernel1_xdir;
int ydim5_advec_cell_kernel1_xdir;


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5


#define OPS_ACC0(x,y,z) (x+xdim0_advec_cell_kernel1_xdir*(y)+xdim0_advec_cell_kernel1_xdir*ydim0_advec_cell_kernel1_xdir*(z))
#define OPS_ACC1(x,y,z) (x+xdim1_advec_cell_kernel1_xdir*(y)+xdim1_advec_cell_kernel1_xdir*ydim1_advec_cell_kernel1_xdir*(z))
#define OPS_ACC2(x,y,z) (x+xdim2_advec_cell_kernel1_xdir*(y)+xdim2_advec_cell_kernel1_xdir*ydim2_advec_cell_kernel1_xdir*(z))
#define OPS_ACC3(x,y,z) (x+xdim3_advec_cell_kernel1_xdir*(y)+xdim3_advec_cell_kernel1_xdir*ydim3_advec_cell_kernel1_xdir*(z))
#define OPS_ACC4(x,y,z) (x+xdim4_advec_cell_kernel1_xdir*(y)+xdim4_advec_cell_kernel1_xdir*ydim4_advec_cell_kernel1_xdir*(z))
#define OPS_ACC5(x,y,z) (x+xdim5_advec_cell_kernel1_xdir*(y)+xdim5_advec_cell_kernel1_xdir*ydim5_advec_cell_kernel1_xdir*(z))

//user function

inline void advec_cell_kernel1_xdir( double *pre_vol, double *post_vol, const double *volume,
                        const double *vol_flux_x, const double *vol_flux_y, const double *vol_flux_z) {

  pre_vol[OPS_ACC0(0,0,0)] = volume[OPS_ACC2(0,0,0)] +
                         ( vol_flux_x[OPS_ACC3(1,0,0)] - vol_flux_x[OPS_ACC3(0,0,0)] +
                           vol_flux_y[OPS_ACC4(0,1,0)] - vol_flux_y[OPS_ACC4(0,0,0)] +
                           vol_flux_z[OPS_ACC5(0,0,1)] - vol_flux_z[OPS_ACC5(0,0,0)]);
  post_vol[OPS_ACC1(0,0,0)] = pre_vol[OPS_ACC0(0,0,0)] - ( vol_flux_x[OPS_ACC3(1,0,0)] - vol_flux_x[OPS_ACC3(0,0,0)]);

}


#undef OPS_ACC0
#undef OPS_ACC1
#undef OPS_ACC2
#undef OPS_ACC3
#undef OPS_ACC4
#undef OPS_ACC5



void advec_cell_kernel1_xdir_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  double *p_a4,
  double *p_a5,
  int x_size, int y_size, int z_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1,p_a2,p_a3,p_a4,p_a5)
  #pragma acc loop
  #endif
  for ( int n_z=0; n_z<z_size; n_z++ ){
    #ifdef OPS_GPU
    #pragma acc loop
    #endif
    for ( int n_y=0; n_y<y_size; n_y++ ){
      #ifdef OPS_GPU
      #pragma acc loop
      #endif
      for ( int n_x=0; n_x<x_size; n_x++ ){
        advec_cell_kernel1_xdir(  p_a0 + n_x*1*1 + n_y*xdim0_advec_cell_kernel1_xdir*1*1 + n_z*xdim0_advec_cell_kernel1_xdir*ydim0_advec_cell_kernel1_xdir*1*1,
           p_a1 + n_x*1*1 + n_y*xdim1_advec_cell_kernel1_xdir*1*1 + n_z*xdim1_advec_cell_kernel1_xdir*ydim1_advec_cell_kernel1_xdir*1*1,
           p_a2 + n_x*1*1 + n_y*xdim2_advec_cell_kernel1_xdir*1*1 + n_z*xdim2_advec_cell_kernel1_xdir*ydim2_advec_cell_kernel1_xdir*1*1,
           p_a3 + n_x*1*1 + n_y*xdim3_advec_cell_kernel1_xdir*1*1 + n_z*xdim3_advec_cell_kernel1_xdir*ydim3_advec_cell_kernel1_xdir*1*1,
           p_a4 + n_x*1*1 + n_y*xdim4_advec_cell_kernel1_xdir*1*1 + n_z*xdim4_advec_cell_kernel1_xdir*ydim4_advec_cell_kernel1_xdir*1*1,
           p_a5 + n_x*1*1 + n_y*xdim5_advec_cell_kernel1_xdir*1*1 + n_z*xdim5_advec_cell_kernel1_xdir*ydim5_advec_cell_kernel1_xdir*1*1 );

      }
    }
  }
}
