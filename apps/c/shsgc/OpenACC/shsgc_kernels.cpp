//
// auto-generated by ops.py//

#include "./OpenACC/shsgc_common.h"

#include <openacc.h>

void ops_init_backend() {acc_set_device_num(ops_get_proc()%acc_get_num_devices(acc_device_nvidia),acc_device_nvidia);}

void ops_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!strcmp(name,"nxp")) {
    nxp = *(int*)dat;
  }
  else
  if (!strcmp(name,"nyp")) {
    nyp = *(int*)dat;
  }
  else
  if (!strcmp(name,"xhalo")) {
    xhalo = *(int*)dat;
  }
  else
  if (!strcmp(name,"yhalo")) {
    yhalo = *(int*)dat;
  }
  else
  if (!strcmp(name,"xmin")) {
    xmin = *(double*)dat;
  }
  else
  if (!strcmp(name,"ymin")) {
    ymin = *(double*)dat;
  }
  else
  if (!strcmp(name,"xmax")) {
    xmax = *(double*)dat;
  }
  else
  if (!strcmp(name,"ymax")) {
    ymax = *(double*)dat;
  }
  else
  if (!strcmp(name,"dx")) {
    dx = *(double*)dat;
  }
  else
  if (!strcmp(name,"dy")) {
    dy = *(double*)dat;
  }
  else
  if (!strcmp(name,"pl")) {
    pl = *(double*)dat;
  }
  else
  if (!strcmp(name,"pr")) {
    pr = *(double*)dat;
  }
  else
  if (!strcmp(name,"rhol")) {
    rhol = *(double*)dat;
  }
  else
  if (!strcmp(name,"rhor")) {
    rhor = *(double*)dat;
  }
  else
  if (!strcmp(name,"ul")) {
    ul = *(double*)dat;
  }
  else
  if (!strcmp(name,"ur")) {
    ur = *(double*)dat;
  }
  else
  if (!strcmp(name,"gam")) {
    gam = *(double*)dat;
  }
  else
  if (!strcmp(name,"gam1")) {
    gam1 = *(double*)dat;
  }
  else
  if (!strcmp(name,"eps")) {
    eps = *(double*)dat;
  }
  else
  if (!strcmp(name,"lambda")) {
    lambda = *(double*)dat;
  }
  else
  if (!strcmp(name,"dt")) {
    dt = *(double*)dat;
  }
  else
  if (!strcmp(name,"del2")) {
    del2 = *(double*)dat;
  }
  else
  if (!strcmp(name,"akap2")) {
    akap2 = *(double*)dat;
  }
  else
  if (!strcmp(name,"tvdsmu")) {
    tvdsmu = *(double*)dat;
  }
  else
  if (!strcmp(name,"con")) {
    con = *(double*)dat;
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "initialize_kernel_openacc_kernel.cpp"
#include "save_kernel_openacc_kernel.cpp"
#include "zerores_kernel_openacc_kernel.cpp"
#include "drhoudx_kernel_openacc_kernel.cpp"
#include "drhouupdx_kernel_openacc_kernel.cpp"
#include "drhoEpudx_kernel_openacc_kernel.cpp"
#include "updateRK3_kernel_openacc_kernel.cpp"
#include "Riemann_kernel_openacc_kernel.cpp"
#include "limiter_kernel_openacc_kernel.cpp"
#include "tvd_kernel_openacc_kernel.cpp"
#include "vars_kernel_openacc_kernel.cpp"
#include "calupwindeff_kernel_openacc_kernel.cpp"
#include "fact_kernel_openacc_kernel.cpp"
#include "update_kernel_openacc_kernel.cpp"
#include "test_kernel_openacc_kernel.cpp"
