//
// auto-generated by ops.py
//

#define OPS_GPU

int xdim0_advec_cell_kernel2_zdir;
int ydim0_advec_cell_kernel2_zdir;
int xdim1_advec_cell_kernel2_zdir;
int ydim1_advec_cell_kernel2_zdir;
int xdim2_advec_cell_kernel2_zdir;
int ydim2_advec_cell_kernel2_zdir;
int xdim3_advec_cell_kernel2_zdir;
int ydim3_advec_cell_kernel2_zdir;

//user function

inline void advec_cell_kernel2_zdir(ptr_double pre_vol,
  ptr_double post_vol,
  const ptr_double volume,
  const ptr_double vol_flux_z) {

  OPS_ACC(pre_vol, 0,0,0) = OPS_ACC(volume, 0,0,0) + OPS_ACC(vol_flux_z, 0,0,1) - OPS_ACC(vol_flux_z, 0,0,0);
  OPS_ACC(post_vol, 0,0,0) = OPS_ACC(volume, 0,0,0);

}


void advec_cell_kernel2_zdir_c_wrapper(
  double *p_a0,
  double *p_a1,
  double *p_a2,
  double *p_a3,
  int x_size, int y_size, int z_size) {
  #ifdef OPS_GPU
  #pragma acc parallel deviceptr(p_a0,p_a1,p_a2,p_a3)
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
        ptr_double ptr0 = {  p_a0 + n_x*1*1 + n_y*xdim0_advec_cell_kernel2_zdir*1*1 + n_z*xdim0_advec_cell_kernel2_zdir*ydim0_advec_cell_kernel2_zdir*1*1, xdim0_advec_cell_kernel2_zdir, ydim0_advec_cell_kernel2_zdir};
        ptr_double ptr1 = {  p_a1 + n_x*1*1 + n_y*xdim1_advec_cell_kernel2_zdir*1*1 + n_z*xdim1_advec_cell_kernel2_zdir*ydim1_advec_cell_kernel2_zdir*1*1, xdim1_advec_cell_kernel2_zdir, ydim1_advec_cell_kernel2_zdir};
        const ptr_double ptr2 = {  p_a2 + n_x*1*1 + n_y*xdim2_advec_cell_kernel2_zdir*1*1 + n_z*xdim2_advec_cell_kernel2_zdir*ydim2_advec_cell_kernel2_zdir*1*1, xdim2_advec_cell_kernel2_zdir, ydim2_advec_cell_kernel2_zdir};
        const ptr_double ptr3 = {  p_a3 + n_x*1*1 + n_y*xdim3_advec_cell_kernel2_zdir*1*1 + n_z*xdim3_advec_cell_kernel2_zdir*ydim3_advec_cell_kernel2_zdir*1*1, xdim3_advec_cell_kernel2_zdir, ydim3_advec_cell_kernel2_zdir};
        advec_cell_kernel2_zdir( ptr0,
          ptr1,
          ptr2,
          ptr3 );

      }
    }
  }
}
