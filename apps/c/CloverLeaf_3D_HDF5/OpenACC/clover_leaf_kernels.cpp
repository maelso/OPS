//
// auto-generated by ops.py//

#include "./OpenACC/clover_leaf_common.h"


void ops_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!strcmp(name,"g_small")) {
    g_small = *(double*)dat;
  }
  else
  if (!strcmp(name,"g_big")) {
    g_big = *(double*)dat;
  }
  else
  if (!strcmp(name,"dtc_safe")) {
    dtc_safe = *(double*)dat;
  }
  else
  if (!strcmp(name,"dtu_safe")) {
    dtu_safe = *(double*)dat;
  }
  else
  if (!strcmp(name,"dtv_safe")) {
    dtv_safe = *(double*)dat;
  }
  else
  if (!strcmp(name,"dtw_safe")) {
    dtw_safe = *(double*)dat;
  }
  else
  if (!strcmp(name,"dtdiv_safe")) {
    dtdiv_safe = *(double*)dat;
  }
  else
  if (!strcmp(name,"field")) {
    field = *(field_type*)dat;
  }
  else
  if (!strcmp(name,"grid")) {
    grid = *(grid_type*)dat;
  }
  else
  if (!strcmp(name,"states")) {
    states = (state_type*)dat;
  }
  else
  if (!strcmp(name,"number_of_states")) {
    number_of_states = *(int*)dat;
  }
  else
  if (!strcmp(name,"g_sphe")) {
    g_sphe = *(int*)dat;
  }
  else
  if (!strcmp(name,"g_point")) {
    g_point = *(int*)dat;
  }
  else
  if (!strcmp(name,"g_cube")) {
    g_cube = *(int*)dat;
  }
  else
  if (!strcmp(name,"dt")) {
    dt = *(double*)dat;
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}
