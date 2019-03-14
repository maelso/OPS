#ifndef WAVE_PROPAGATION_KERNELS_H
#define WAVE_PROPAGATION_KERNELS_H

void initialize_velocity_model_kernel(float *m){
    m[OPS_ACC0(0, 0, 0)] = 0.25;
}

void initialize_damp_kernel(float *damp, int *idx){
    float dampcoeff;
    float pos;
    float val;
    int pivot;

    int damp_left_initial = space_order/2;
    int damp_left_final = space_order/2 + border_size;
    int damp_right_initial = damp_left_final + X_size;
    int damp_right_final = damp_right_initial + border_size;

    dampcoeff = 0.2590408229618301; // 1.5 * log(1.0 / 0.001) / (40.)
    damp[OPS_ACC0(0, 0, 0)] = 0.0;

    if ( (idx[0] >= damp_left_initial && idx[1] >= damp_left_initial && idx[2] >= damp_left_initial) && 
         (idx[0] < damp_right_final && idx[1] < damp_right_final && idx[2] < damp_right_final)
    ){
        //IDX 2
        if(idx[2] >= damp_left_initial && idx[2] < damp_left_final){
            pivot = idx[2];
            pos = abs((border_size - pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
        else if(idx[2] > damp_right_initial && idx[2] < damp_right_final){
            pivot = idx[2] % (border_size + 1);
            pos = abs((pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
        //IDX 1
        else if(idx[1] >= damp_left_initial && idx[1] < damp_left_final){
            pivot = idx[1];
            pos = abs((border_size - pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
        else if(idx[1] > damp_right_initial && idx[1] < damp_right_final){
            pivot = idx[1] % (border_size + 1);
            pos = abs((pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
        //IDX 0
        else if(idx[0] >= damp_left_initial && idx[0] < damp_left_final){
            pivot = idx[0];
            pos = abs((border_size - pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
        else if(idx[0] > damp_right_initial && idx[0] < damp_right_final){
            pivot = idx[0] % (border_size + 1);
            pos = abs((pivot + 2) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));
            damp[OPS_ACC0(0, 0, 0)] = pos;
        }
    }
}

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
