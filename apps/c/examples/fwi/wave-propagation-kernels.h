void Kernel0(const float * ut0, float * ut1)
{
  ut1[OPS_ACC1(0,0)] = ut0[OPS_ACC0(0,0)] + 2;
}