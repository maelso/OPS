//
// auto-generated by ops.py
//

int xdim0_tea_leaf_norm2_kernel;


//user function

void tea_leaf_norm2_kernel_c_wrapper(double *restrict x_p,
                                     double *restrict norm_g, int x_size,
                                     int y_size) {
  double norm_0 = norm_g[0];
  #pragma omp parallel for reduction(+:norm_0)
  for ( int n_y=0; n_y<y_size; n_y++ ){
    for ( int n_x=0; n_x<x_size; n_x++ ){
      double norm[1];
      norm[0] = ZERO_double;
      const ptr_double x = {x_p + n_x * 1 +
                                n_y * xdim0_tea_leaf_norm2_kernel * 1,
                            xdim0_tea_leaf_norm2_kernel};

      *norm = *norm + OPS_ACC(x, 0, 0) * OPS_ACC(x, 0, 0);

      norm_0 +=norm[0];
    }
  }
  norm_g[0] = norm_0;
}
