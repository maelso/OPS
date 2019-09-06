#ifndef WAVE_PROPAGATION_KERNELS_H
#define WAVE_PROPAGATION_KERNELS_H

void initialize_velocity_model_kernel(float *m)
{
    m[OPS_ACC0(0, 0, 0)] = 0.25;
}

void set_zero_kernel(float *damp)
{
    damp[OPS_ACC0(0, 0, 0)] = 0.0;
}

void wave_propagation_2th_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    // ut10[OPS_ACC2(0, 0, 0)] = 1.0;
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)];
    float dt3r0 = dt * dt * dt / r0;
    ut10[OPS_ACC2(0, 0, 0)] = 2.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 +
                              1.0 * (-dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                                     dt3r0 *
                                         (ut00[OPS_ACC3(0, 0, -1)] +
                                          ut00[OPS_ACC3(0, -1, 0)] +
                                          ut00[OPS_ACC3(-1, 0, 0)] +
                                          ut00[OPS_ACC3(1, 0, 0)] +
                                          ut00[OPS_ACC3(0, 1, 0)] +
                                          ut00[OPS_ACC3(0, 0, 1)])) +
                              0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              6.0 * dt3r0 * ut00[OPS_ACC3(0, 0, 0)];
}

void wave_propagation_4th_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)];
    float dt3r0 = dt * dt * dt / r0;
    ut10[OPS_ACC2(0, 0, 0)] = 2.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 -
                              1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                              1.33333333 * dt3r0 *
                                  (ut00[OPS_ACC3(-1, 0, 0)] +
                                   ut00[OPS_ACC3(0, -1, 0)] +
                                   ut00[OPS_ACC3(0, 0, -1)] +
                                   ut00[OPS_ACC3(0, 0, 1)] +
                                   ut00[OPS_ACC3(0, 1, 0)] +
                                   ut00[OPS_ACC3(1, 0, 0)]) -
                              0.0833333333 * dt3r0 *
                                  (ut00[OPS_ACC3(-2, 0, 0)] +
                                   ut00[OPS_ACC3(0, -2, 0)] +
                                   ut00[OPS_ACC3(0, 0, -2)] +
                                   ut00[OPS_ACC3(0, 0, 2)] +
                                   ut00[OPS_ACC3(0, 2, 0)] +
                                   ut00[OPS_ACC3(2, 0, 0)]) +
                              0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              7.5 * dt3r0 * ut00[OPS_ACC3(0, 0, 0)];
}

void wave_propagation_8th_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)];
    float dt3r0 = dt * dt * dt / r0;
    ut10[OPS_ACC2(0, 0, 0)] = 2.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 -
                              1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                              1.6 * dt3r0 *
                                  (ut00[OPS_ACC3(-1, 0, 0)] +
                                   ut00[OPS_ACC3(0, -1, 0)] +
                                   ut00[OPS_ACC3(0, 0, -1)] +
                                   ut00[OPS_ACC3(0, 0, 1)] +
                                   ut00[OPS_ACC3(0, 1, 0)] +
                                   ut00[OPS_ACC3(1, 0, 0)]) -
                              0.2 * dt3r0 *
                                  (ut00[OPS_ACC3(-2, 0, 0)] +
                                   ut00[OPS_ACC3(0, -2, 0)] +
                                   ut00[OPS_ACC3(0, 0, -2)] +
                                   ut00[OPS_ACC3(0, 0, 2)] +
                                   ut00[OPS_ACC3(0, 2, 0)] +
                                   ut00[OPS_ACC3(2, 0, 0)]) +
                              0.025 * dt3r0 *
                                  (ut00[OPS_ACC3(-3, 0, 0)] +
                                   ut00[OPS_ACC3(0, -3, 0)] +
                                   ut00[OPS_ACC3(0, 0, -3)] +
                                   ut00[OPS_ACC3(0, 0, 3)] +
                                   ut00[OPS_ACC3(0, 3, 0)] +
                                   ut00[OPS_ACC3(3, 0, 0)]) -
                              0.0017 * dt3r0 *
                                  (ut00[OPS_ACC3(-4, 0, 0)] +
                                   ut00[OPS_ACC3(0, -4, 0)] +
                                   ut00[OPS_ACC3(0, 0, -4)] +
                                   ut00[OPS_ACC3(0, 0, 4)] +
                                   ut00[OPS_ACC3(0, 4, 0)] +
                                   ut00[OPS_ACC3(4, 0, 0)]) +
                              0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              8.5 * dt3r0 * ut00[OPS_ACC3(0, 0, 0)];
}

void wave_propagation_16th_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)];
    float dt3r0 = dt * dt * dt / r0;
    ut10[OPS_ACC2(0, 0, 0)] = 1.777 * dt3r0 *
                                  (ut00[OPS_ACC3(-1, 0, 0)] +
                                   ut00[OPS_ACC3(0, -1, 0)] +
                                   ut00[OPS_ACC3(0, 0, -1)] +
                                   ut00[OPS_ACC3(0, 0, 1)] +
                                   ut00[OPS_ACC3(0, 1, 0)] +
                                   ut00[OPS_ACC3(1, 0, 0)]) -
                              0.3111 * dt3r0 *
                                  (ut00[OPS_ACC3(-2, 0, 0)] +
                                   ut00[OPS_ACC3(0, -2, 0)] +
                                   ut00[OPS_ACC3(0, 0, -2)] +
                                   ut00[OPS_ACC3(0, 0, 2)] +
                                   ut00[OPS_ACC3(0, 2, 0)] +
                                   ut00[OPS_ACC3(2, 0, 0)]) +
                              0.0754 * dt3r0 *
                                  (ut00[OPS_ACC3(-3, 0, 0)] +
                                   ut00[OPS_ACC3(0, -3, 0)] +
                                   ut00[OPS_ACC3(0, 0, -3)] +
                                   ut00[OPS_ACC3(0, 0, 3)] +
                                   ut00[OPS_ACC3(0, 3, 0)] +
                                   ut00[OPS_ACC3(3, 0, 0)]) -
                              0.017676 * dt3r0 *
                                  (ut00[OPS_ACC3(-4, 0, 0)] +
                                   ut00[OPS_ACC3(0, -4, 0)] +
                                   ut00[OPS_ACC3(0, 0, -4)] +
                                   ut00[OPS_ACC3(0, 0, 4)] +
                                   ut00[OPS_ACC3(0, 4, 0)] +
                                   ut00[OPS_ACC3(4, 0, 0)]) +
                              0.00348 * dt3r0 *
                                  (ut00[OPS_ACC3(-5, 0, 0)] +
                                   ut00[OPS_ACC3(0, -5, 0)] +
                                   ut00[OPS_ACC3(0, 0, -5)] +
                                   ut00[OPS_ACC3(0, 0, 5)] +
                                   ut00[OPS_ACC3(0, 5, 0)] +
                                   ut00[OPS_ACC3(5, 0, 0)]) -
                              5.18000518e-4F * dt3r0 *
                                  (ut00[OPS_ACC3(-6, 0, 0)] +
                                   ut00[OPS_ACC3(0, -6, 0)] +
                                   ut00[OPS_ACC3(0, 0, -6)] +
                                   ut00[OPS_ACC3(0, 0, 6)] +
                                   ut00[OPS_ACC3(0, 6, 0)] +
                                   ut00[OPS_ACC3(6, 0, 0)]) +
                              5.07429079e-5F * dt3r0 *
                                  (ut00[OPS_ACC3(-7, 0, 0)] +
                                   ut00[OPS_ACC3(0, -7, 0)] +
                                   ut00[OPS_ACC3(0, 0, -7)] +
                                   ut00[OPS_ACC3(0, 0, 7)] +
                                   ut00[OPS_ACC3(0, 7, 0)] +
                                   ut00[OPS_ACC3(7, 0, 0)]) -
                              2.42812743e-6F * dt3r0 *
                                  (ut00[OPS_ACC3(-8, 0, 0)] +
                                   ut00[OPS_ACC3(0, -8, 0)] +
                                   ut00[OPS_ACC3(0, 0, -8)] +
                                   ut00[OPS_ACC3(0, 0, 8)] +
                                   ut00[OPS_ACC3(0, 8, 0)] +
                                   ut00[OPS_ACC3(8, 0, 0)]) +
                              2.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 -
                              1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                              0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              9.164 * dt3r0 * ut00[OPS_ACC3(0, 0, 0)];
}

void wave_propagation_32th_kernel(const float *dampx0, const float *mx0, float *ut10, const float *ut00, const float *ut20, int *idx)
{
    float r0 = 1.0 * dt * mx0[OPS_ACC1(0, 0, 0)] + 0.5 * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)];
    float dt3r0 = dt * dt * dt / r0;
    ut10[OPS_ACC2(0, 0, 0)] = 1.77777778F * dt3r0 *
                                  (ut00[OPS_ACC3(-1, 0, 0)] +
                                   ut00[OPS_ACC3(0, -1, 0)] +
                                   ut00[OPS_ACC3(0, 0, -1)] +
                                   ut00[OPS_ACC3(0, 0, 1)] +
                                   ut00[OPS_ACC3(0, 1, 0)] +
                                   ut00[OPS_ACC3(1, 0, 0)]) -
                              3.11111111e-1F * dt3r0 *
                                  (ut00[OPS_ACC3(-2, 0, 0)] +
                                   ut00[OPS_ACC3(0, -2, 0)] +
                                   ut00[OPS_ACC3(0, 0, -2)] +
                                   ut00[OPS_ACC3(0, 0, 2)] +
                                   ut00[OPS_ACC3(0, 2, 0)] +
                                   ut00[OPS_ACC3(2, 0, 0)]) +
                              7.54208754e-2F * dt3r0 *
                                  (ut00[OPS_ACC3(-3, 0, 0)] +
                                   ut00[OPS_ACC3(0, -3, 0)] +
                                   ut00[OPS_ACC3(0, 0, -3)] +
                                   ut00[OPS_ACC3(0, 0, 3)] +
                                   ut00[OPS_ACC3(0, 3, 0)] +
                                   ut00[OPS_ACC3(3, 0, 0)]) -
                              1.76767677e-2F * dt3r0 *
                                  (ut00[OPS_ACC3(-4, 0, 0)] +
                                   ut00[OPS_ACC3(0, -4, 0)] +
                                   ut00[OPS_ACC3(0, 0, -4)] +
                                   ut00[OPS_ACC3(0, 0, 4)] +
                                   ut00[OPS_ACC3(0, 4, 0)] +
                                   ut00[OPS_ACC3(4, 0, 0)]) +
                              3.48096348e-3F * dt3r0 *
                                  (ut00[OPS_ACC3(-5, 0, 0)] +
                                   ut00[OPS_ACC3(0, -5, 0)] +
                                   ut00[OPS_ACC3(0, 0, -5)] +
                                   ut00[OPS_ACC3(0, 0, 5)] +
                                   ut00[OPS_ACC3(0, 5, 0)] +
                                   ut00[OPS_ACC3(5, 0, 0)]) -
                              5.18000518e-4F * dt3r0 *
                                  (ut00[OPS_ACC3(-6, 0, 0)] +
                                   ut00[OPS_ACC3(0, -6, 0)] +
                                   ut00[OPS_ACC3(0, 0, -6)] +
                                   ut00[OPS_ACC3(0, 0, 6)] +
                                   ut00[OPS_ACC3(0, 6, 0)] +
                                   ut00[OPS_ACC3(6, 0, 0)]) +
                              5.07429079e-5F * dt3r0 *
                                  (ut00[OPS_ACC3(-7, 0, 0)] +
                                   ut00[OPS_ACC3(0, -7, 0)] +
                                   ut00[OPS_ACC3(0, 0, -7)] +
                                   ut00[OPS_ACC3(0, 0, 7)] +
                                   ut00[OPS_ACC3(0, 7, 0)] +
                                   ut00[OPS_ACC3(7, 0, 0)]) -
                              2.42812743e-6F * dt3r0 *
                                  (ut00[OPS_ACC3(-8, 0, 0)] +
                                   ut00[OPS_ACC3(0, -8, 0)] +
                                   ut00[OPS_ACC3(0, 0, -8)] +
                                   ut00[OPS_ACC3(0, 0, 8)] +
                                   ut00[OPS_ACC3(0, 8, 0)] +
                                   ut00[OPS_ACC3(8, 0, 0)]) +
                              2.0F * dt * mx0[OPS_ACC1(0, 0, 0)] * ut00[OPS_ACC3(0, 0, 0)] / r0 -
                              1.0F * dt * mx0[OPS_ACC1(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 +
                              5.0e-1F * (dt * dt) * dampx0[OPS_ACC0(0, 0, 0)] * ut20[OPS_ACC4(0, 0, 0)] / r0 -
                              9.16453231F * dt3r0 * ut00[OPS_ACC3(0, 0, 0)];
}
void initialize_damp_kernel(float *damp, int *idx)
{
    float dampcoeff;
    float pos;
    float val;
    int pivot;

    int damp_left_initial = space_order / 2;
    int damp_left_final = damp_left_initial + border_size;
    int damp_right_initial = damp_left_final + X_size;
    int damp_right_final = damp_right_initial + border_size;

    dampcoeff = 0.2590408229618301; // 1.5 * log(1.0 / 0.001) / (40.)

    for (int dim = 0; dim < 3; dim++)
    {
        if ((idx[0] >= damp_left_initial && idx[1] >= damp_left_initial && idx[2] >= damp_left_initial) &&
            (idx[0] < damp_right_final && idx[1] < damp_right_final && idx[2] < damp_right_final))
        {

            if (idx[dim] >= damp_left_initial && idx[dim] < damp_left_final)
            {
                pivot = idx[dim];
                pos = abs((border_size - pivot + 2) / (float)border_size);
                val = dampcoeff * (pos - sin(2 * 3.14159265358979323846 * pos) / (2 * 3.14159265358979323846));
                damp[OPS_ACC0(0, 0, 0)] += val;
            }
            else if (idx[dim] > damp_right_initial && idx[dim] < damp_right_final)
            {
                pivot = idx[dim] - damp_right_initial - 1;
                pos = abs((pivot + 2) / (float)border_size);
                val = dampcoeff * (pos - sin(2 * 3.14159265358979323846 * pos) / (2 * 3.14159265358979323846));
                damp[OPS_ACC0(0, 0, 0)] += val;
            }
        }
    }
}

void set_space_order_border_kernel(float *damp, const float *damp2, int *idx)
{
    int so_border_left = space_order / 2;
    int so_border_right = space_order / 2 + border_size * 2 + X_size;
    int x, y, z;

    if (idx[0] < so_border_left)
        x = so_border_left;
    else if (idx[0] >= so_border_right)
        x = -so_border_left;
    else
        x = 0;
    if (idx[1] < so_border_left)
        y = so_border_left;
    else if (idx[1] >= so_border_right)
        y = -so_border_left;
    else
        y = 0;
    if (idx[2] < so_border_left)
        z = so_border_left;
    else if (idx[2] >= so_border_right)
        z = -so_border_left;
    else
        z = 0;

    damp[OPS_ACC0(0, 0, 0)] = damp2[OPS_ACC1(x, y, z)];
}

void source_injection_kernel(float *u, const float *m, const float *src_value, const int *idx)
{
    if (idx[2] == ii_src[0] + 2 && idx[1] == ii_src[1] + 2 && idx[0] == ii_src[2] + 2)
    {
        float r7 = 2.2801e-2F *
                   (-1.0F * p[0] * p[1] * p[2] + 1.0F * p[0] * p[1] + 1.0F * p[0] * p[2] - 1.0F * p[0] +
                    1.0F * p[1] * p[2] - 1.0F * p[1] - 1.0F * p[2] + 1) *
                   src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r7);
        u[OPS_ACC0(0, 0, 0)] += r7;
    }
    else if (idx[2] == ii_src[0] + 2 && idx[1] == ii_src[1] + 2 && idx[0] == ii_src[3] + 2)
    {
        float r11 = 2.2801e-2F * (1.0F * p[0] * p[1] * p[2] - 1.0F * p[0] * p[2] - 1.0F * p[1] * p[2] + 1.0F * p[2]) *
                    src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r11);
        u[OPS_ACC0(0, 0, 0)] += r11;
    }
    else if (idx[2] == ii_src[0] + 2 && idx[1] == ii_src[4] + 2 && idx[0] == ii_src[2] + 2)
    {
        float r15 = 2.2801e-2F * (1.0F * p[0] * p[1] * p[2] - 1.0F * p[0] * p[1] - 1.0F * p[1] * p[2] + 1.0F * p[1]) *
                    src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r15);
        u[OPS_ACC0(0, 0, 0)] += r15;
    }
    else if (idx[2] == ii_src[0] + 2 && idx[1] == ii_src[4] + 2 && idx[0] == ii_src[3] + 2)
    {
        float r19 =
            2.2801e-2F * (-1.0F * p[0] * p[1] * p[2] + 1.0F * p[1] * p[2]) * src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r19);
        u[OPS_ACC0(0, 0, 0)] += r19;
    }
    else if (idx[2] == ii_src[5] + 2 && idx[1] == ii_src[1] + 2 && idx[0] == ii_src[2] + 2)
    {
        float r23 = 2.2801e-2F * (1.0F * p[0] * p[1] * p[2] - 1.0F * p[0] * p[1] - 1.0F * p[0] * p[2] + 1.0F * p[0]) *
                    src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r23);
        u[OPS_ACC0(0, 0, 0)] += r23;
    }
    else if (idx[2] == ii_src[5] + 2 && idx[1] == ii_src[1] + 2 && idx[0] == ii_src[3] + 2)
    {
        float r27 =
            2.2801e-2F * (-1.0F * p[0] * p[1] * p[2] + 1.0F * p[0] * p[2]) * src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r27);
        u[OPS_ACC0(0, 0, 0)] += r27;
    }
    else if (idx[2] == ii_src[5] + 2 && idx[1] == ii_src[4] + 2 && idx[0] == ii_src[2] + 2)
    {
        float r31 =
            2.2801e-2F * (-1.0F * p[0] * p[1] * p[2] + 1.0F * p[0] * p[1]) * src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r31);
        u[OPS_ACC0(0, 0, 0)] += r31;
    }
    else if (idx[2] == ii_src[5] + 2 && idx[1] == ii_src[4] + 2 && idx[0] == ii_src[3] + 2)
    {
        float r35 = 2.2801e-2F * p[0] * p[1] * p[2] * src_value[0] / m[OPS_ACC1(0, 0, 0)];
        // printf("u[%d][%d][%d]=%f\n", idx[0], idx[1], idx[2], r35);
        u[OPS_ACC0(0, 0, 0)] += r35;
    }
}

#endif
