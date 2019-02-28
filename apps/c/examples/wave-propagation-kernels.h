#ifndef WAVE_PROPAGATION_KERNELS_H
#define WAVE_PROPAGATION_KERNELS_H

void wave_propagation_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    // ut10[OPS_ACC2(0, 0, 0)] = 1.0;
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(-1, -1, -1)];
    
    ut10[OPS_ACC2(0, 0, 0)] = 2.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 +
                              1.0 * (-dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(0, 0, -1)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(0, -1, 0)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(-1, 0, 0)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(1, 0, 0)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(0, 1, 0)] / r0 +
                                     (dt * dt * dt) * ut00[OPS_ACC3(0, 0, 1)] / r0) +
                              0.5 * (dt * dt) * dampx0[OPS_ACC0(-1, -1, -1)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              6.0 * dt * dt * dt * ut00[OPS_ACC3(0, 0, 0)] / r0;

    // if (ut10[OPS_ACC2(0, 0, 0)] != 0)
    // {
    //        printf("Injecting u[%d][%d][%d]=%f ", idx[2], idx[1], idx[0], ut10[OPS_ACC2(0, 0, 0)]);
    //        printf("r0=%f | dt=%f | m=%f | damp=%f ", r0, dt, mx0[OPS_ACC1(0, 0, 0)], dampx0[OPS_ACC0(-1, -1, -1)]);
    //        printf("\t ut2=%f | ut0=%f %f %f %f %f %f %f\n", ut20[OPS_ACC4(0, 0, 0)],
    //                                                         ut00[OPS_ACC3(0, 0, 0)],
    //                                                         ut00[OPS_ACC3(0, 0, -1)],
    //                                                         ut00[OPS_ACC3(0, -1, 0)],
    //                                                         ut00[OPS_ACC3(-1, 0, 0)],
    //                                                         ut00[OPS_ACC3(1, 0, 0)],
    //                                                         ut00[OPS_ACC3(0, 1, 0)],
    //                                                         ut00[OPS_ACC3(0, 0, 1)]);
    // }
}

#endif
