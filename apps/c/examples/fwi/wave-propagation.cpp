#define _POSIX_C_SOURCE 200809L
#define OPS_2D
#include "stdlib.h"
#include "math.h"
#include "sys/time.h"
#include "ops_seq.h"

#include "wave-propagation-kernels.h"

struct dataobj
{
  void * data;
  int * size;
  int * npsize;
  int * dsize;
  int * hsize;
  int * hofs;
  int * oofs;
} ;

struct profiler
{
  double section0;
} ;

int main(const int x_M, const int x_m, const int y_M, const int y_m, const int time_M, const int time_m, struct profiler * timers, struct dataobj * u_vec)
{
  float (* u)[u_vec->size[1]][u_vec->size[2]] __attribute__ ((aligned (64))) = (float (*)[u_vec->size[1]][u_vec->size[2]]) u_vec->data;
  ops_init(0,0,6);
  int range_0[4] = {x_m, x_M, y_m, y_M};
  ops_block block_0 = ops_decl_block(2,"block_0");
  int s2d_ut1_1pt[2] = {0, 0};
  ops_stencil S2D_UT1_1PT = ops_decl_stencil(2,1,(int *)s2d_ut1_1pt,"S2D_UT1_1PT");
  int s2d_ut0_1pt[2] = {0, 0};
  ops_stencil S2D_UT0_1PT = ops_decl_stencil(2,1,(int *)s2d_ut0_1pt,"S2D_UT0_1PT");
  int u_dim[2] = {3, 3};
  int u_base[2] = {0, 0};
  int u_d_p[2] = {0, 0};
  int u_d_m[2] = {0, 0};
  ops_dat u_dat[2];
  u_dat[0] = ops_decl_dat(block_0,1,(int *)u_dim,(int *)u_base,(int *)u_d_m,(int *)u_d_p,(float *)&u[0],"float","ut0");
  u_dat[1] = ops_decl_dat(block_0,1,(int *)u_dim,(int *)u_base,(int *)u_d_m,(int *)u_d_p,(float *)&u[1],"float","ut1");
  ops_partition("");
  for (int time = time_m, t0 = (time)%(2), t1 = (time + 1)%(2); time <= time_M; time += 1, t0 = (time)%(2), t1 = (time + 1)%(2))
  {
    /* Begin section0 */
    ops_par_loop(Kernel0,"Kernel0",block_0,2,(int *)range_0,ops_arg_dat(u_dat[t0],1,S2D_UT0_1PT,"float",OPS_READ),ops_arg_dat(u_dat[t1],1,S2D_UT1_1PT,"float",OPS_WRITE));
    /* End section0 */
  }
  ops_timing_output(stdout);
  ops_exit();
  return 0;
}
