#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define OPS_3D
#include "ops_seq.h"

int X_size = 40;
int Y_size = 40;
int Z_size = 40;
int padding = 2;
double dt = 0.001, start = 0, stop = 30; // time variables
int border_size = 10;                    // Abosrbent border
int space_order = 0;                     // Space order

#include "wave-propagation-kernels.h"

void initialize_velocity_model(float *m, int x_size, int y_size, int z_size);
void initialize_damp(float *damp, int x_size, int y_size, int z_size);
void initialize_source(float *src, int T_intervals);
void print_vector(float *vec, int x_limit, int y_limit, int z_limit);

int main(int argc, char *argv[])
{
    int t = 1; // no devito, ele comeca no intervalo 1.
    int size[3];
    int damp_size[3];
    int base[] = {0, 0, 0};
    int d_m[] = {-padding, -padding, -padding};
    int d_p[] = {padding, padding, padding};
    // int d_m_0[] = {0, 0, 0};
    // int d_p_0[] = {0, 0, 0};
    float **u, *m, *damp, *src, src_coords[3];
    int T_intervals;
    char title[25];

    // Esse size leva em consideracao a borda absorvente???
    size[0] = X_size + 2 * border_size + 2 * space_order;
    size[1] = Y_size + 2 * border_size + 2 * space_order;
    size[2] = Z_size + 2 * border_size + 2 * space_order;
    damp_size[0] = X_size + 2 * border_size + 2;
    damp_size[1] = Y_size + 2 * border_size + 2;
    damp_size[2] = Z_size + 2 * border_size + 2;

    T_intervals = ceil((stop - start + dt) / dt);
    // ops_printf("T_intervals = %d\n", T_intervals);

    ops_init(argc, NULL, 1);

    // Declare global constant
    ops_decl_const("X_size", 1, "int", &X_size);
    ops_decl_const("Y_size", 1, "int", &Y_size);
    ops_decl_const("Z_size", 1, "int", &Z_size);
    ops_decl_const("border_size", 1, "int", &border_size);
    ops_decl_const("dt", 1, "double", &dt);

    // Alocates and initialize grid
    u = (float **)malloc(3 * sizeof(float *));
    u[0] = (float *)calloc((X_size + 2 * border_size + 2 * padding) * (Y_size + 2 * border_size + 2 * padding) * (Z_size + 2 * border_size + 2 * padding), sizeof(float));
    u[1] = (float *)calloc((X_size + 2 * border_size + 2 * padding) * (Y_size + 2 * border_size + 2 * padding) * (Z_size + 2 * border_size + 2 * padding), sizeof(float));
    u[2] = (float *)calloc((X_size + 2 * border_size + 2 * padding) * (Y_size + 2 * border_size + 2 * padding) * (Z_size + 2 * border_size + 2 * padding), sizeof(float));
    m = (float *)malloc((X_size + 2 * border_size + 2 * padding) * (Y_size + 2 * border_size + 2 * padding) * (Z_size + 2 * border_size + 2 * padding) * sizeof(float));
    damp = (float *)calloc((X_size + 2 * border_size + 1 + 1) * (Y_size + 2 * border_size + 1 + 1) * (Z_size + 2 * border_size + 1 + 1), sizeof(float));
    src = (float *)malloc(T_intervals * sizeof(float));

    // Initialize velocity model
    initialize_velocity_model(m, X_size + 2 * border_size + 2 * padding, Y_size + 2 * border_size + 2 * padding, Z_size + 2 * border_size + 2 * padding);
    // print_vector(m, X_size + 2 * border_size + 2 * padding, Y_size + 2 * border_size + 2 * padding, Z_size + 2 * border_size + 2 * padding);

    // Initialize Damp
    initialize_damp(damp, X_size + 2 * border_size + 1 + 1, Y_size + 2 * border_size + 1 + 1, Z_size + 2 * border_size + 1 + 1);
    // print_vector(damp, X_size + 2 * border_size + 1 + 1, Y_size + 2 * border_size + 1 + 1, Z_size + 2 * border_size + 1 + 1);

    // Initialize source
    initialize_source(src, T_intervals);
    // print_vector(src, T_intervals, 0, 0);

    // Initialize source position
    src_coords[0] = (float)(X_size - 1) / (float)2;
    src_coords[1] = (float)(Y_size - 1) / (float)2;
    src_coords[2] = (float)1.;

    ops_printf("Src coordinates:%f %f %f\n", src_coords[0], src_coords[1], src_coords[2]);

    // Declare ops_block
    ops_block grid = ops_decl_block(3, "grid");

    // Declare ops_dat objects
    ops_dat dat_ut[3];
    ops_dat dat_m;
    ops_dat dat_damp;
    dat_ut[0] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[0], "float", "ut0");
    dat_ut[1] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[1], "float", "ut1");
    dat_ut[2] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[2], "float", "ut2");
    dat_m = ops_decl_dat(grid, 1, size, base, d_m, d_p, m, "float", "m");
    dat_damp = ops_decl_dat(grid, 1, damp_size, base, d_m, d_p, damp, "float", "damp");

    int s3d_000[] = {0, 0, 0};
    int s3d_1pt[] = {-1, -1, -1};
    int s3d_7pts[] = {0, 0, 0,
                      -1, 0, 0,
                      0, -1, 0,
                      0, 0, -1,
                      0, 0, 1,
                      0, 1, 0,
                      1, 0, 0};

    // Declare stencil
    ops_stencil S3D_000 = ops_decl_stencil(3, 1, s3d_000, "0,0,0");
    ops_stencil S3D_1PT = ops_decl_stencil(3, 1, s3d_1pt, "-1,-1,-1");
    ops_stencil S3D_7PTS = ops_decl_stencil(3, 7, s3d_7pts, "7pts");

    int shift = space_order;
    int whole_range[] = {space_order, X_size + 2 * border_size, space_order, Y_size + 2 * border_size, shift, Z_size + 2 * border_size};

    ops_printf("Starts time propagation.\n");

    int disp[3], sizes[3], strides[3];

    // ops_printf("time=%d\n", t);
    // print_vector(u[(t + 1) % 3], 34, 34, 34);
    // ops_print_dat_to_txtfile(dat_ut[t % 3], "output/u_ops000.txt");
    
    ops_partition("");

    ops_printf("\n***********************************************************\n");
    do
    {
        // ops_printf("%d\n", t);
        ops_par_loop(wave_propagation_kernel, "wave_propagation_kernel", grid, 3, whole_range,
                     ops_arg_dat(dat_damp, 1, S3D_1PT, "float", OPS_READ),
                     ops_arg_dat(dat_m, 1, S3D_000, "float", OPS_READ),
                     ops_arg_dat(dat_ut[(t + 1) % 3], 1, S3D_000, "float", OPS_WRITE), // t1
                     ops_arg_dat(dat_ut[t % 3], 1, S3D_7PTS, "float", OPS_READ),       // t0
                     ops_arg_dat(dat_ut[(t + 2) % 3], 1, S3D_000, "float", OPS_READ),  // t2
                     ops_arg_idx());

        ops_dat_get_extents(dat_ut[(t + 1) % 3], 0, disp, sizes);
        u[(t + 1) % 3] = (float *)ops_dat_get_raw_pointer(dat_ut[(t + 1) % 3], 0, S3D_000, strides);

        // ops_printf("disp=(%d,%d,%d)\n", disp[0], disp[1], disp[2]);
        // ops_printf("size=(%d,%d,%d)\n", sizes[0], sizes[1], sizes[2]);
        // ops_printf("strides=(%d,%d,%d)\n", strides[0], strides[1], strides[2]);

        // source injection interpolation
        float r1 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[0]));
        int ii_src_0 = r1;
        float r2 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[1]));
        int ii_src_1 = r2;
        float r3 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[2]));
        int ii_src_2 = r3;
        int ii_src_3 = r3 + 1;
        int ii_src_4 = r2 + 1;
        int ii_src_5 = r1 + 1;
        float px = (float)(border_size - 1.0F * r1 + src_coords[0]);
        float py = (float)(border_size - 1.0F * r2 + src_coords[1]);
        float pz = (float)(border_size - 1.0F * r3 + src_coords[2]);

        int max_index = X_size + 2 * border_size + 2 * padding;
        int u_dim_size = X_size + 2 * border_size + 2 * padding;

        // ops_printf("src=%f\n", src[t]);
        // ops_printf("r1=%f r2=%f r3=%f\n", r1, r2, r3);
        if (ii_src_0 >= -1 && ii_src_1 >= -1 && ii_src_2 >= -1 && ii_src_0 <= max_index && ii_src_1 <= max_index && ii_src_2 <= max_index)
        {
            int r4 = ii_src_0 + 2;
            int r5 = ii_src_1 + 2;
            int r6 = ii_src_2 + 2;
            float r7 = 2.2801e-2F * (-1.0F * px * py * pz + 1.0F * px * py + 1.0F * px * pz - 1.0F * px + 1.0F * py * pz - 1.0F * py - 1.0F * pz + 1) * src[t] / m[r4 * max_index * max_index + r5 * max_index + r6];
            u[((t + 1) % 3)][r4 * u_dim_size * u_dim_size +
                             r5 * u_dim_size +
                             r6] += r7;
            // ops_printf("u[%d][%d][%d]=%f | m=%f | src[%d]=%f\n", r4, r5, r6, r7, m[r4 * max_index * max_index + r5 * max_index + r6], t, src[t]);
        }
        if (ii_src_0 >= -1 && ii_src_1 >= -1 && ii_src_3 >= -1 && ii_src_0 <= max_index && ii_src_1 <= max_index && ii_src_3 <= max_index)
        {
            int r8 = ii_src_0 + 2;
            int r9 = ii_src_1 + 2;
            int r10 = ii_src_3 + 2;
            float r11 = 2.2801e-2F * (1.0F * px * py * pz - 1.0F * px * pz - 1.0F * py * pz + 1.0F * pz) * src[t] / m[r8 * max_index * max_index + r9 * max_index + r10];
            u[(t + 1) % 3][r8 * u_dim_size * u_dim_size +
                           r9 * u_dim_size +
                           r10] += r11;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r8, r9, r10, r11, m[r8 * max_index * max_index + r9 * max_index + r10]);
        }

        if (ii_src_0 >= -1 && ii_src_2 >= -1 && ii_src_4 >= -1 && ii_src_0 <= max_index && ii_src_2 <= max_index && ii_src_4 <= max_index)
        {
            int r12 = ii_src_0 + 2;
            int r13 = ii_src_4 + 2;
            int r14 = ii_src_2 + 2;
            float r15 = 2.2801e-2F * (1.0F * px * py * pz - 1.0F * px * py - 1.0F * py * pz + 1.0F * py) * src[t] / m[r12 * max_index * max_index + r13 * max_index + r14];
            u[(t + 1) % 3][r12 * u_dim_size * u_dim_size +
                           r13 * u_dim_size +
                           r14] += r15;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r12, r13, r14, r15, m[r12 * max_index * max_index + r13 * max_index + r14]);
        }

        if (ii_src_0 >= -1 && ii_src_3 >= -1 && ii_src_4 >= -1 && ii_src_0 <= max_index && ii_src_3 <= max_index && ii_src_4 <= max_index)
        {
            int r16 = ii_src_0 + 2;
            int r17 = ii_src_4 + 2;
            int r18 = ii_src_3 + 2;
            float r19 = 2.2801e-2F * (-1.0F * px * py * pz + 1.0F * py * pz) * src[t] / m[r16 * max_index * max_index + r17 * max_index + r18];
            u[(t + 1) % 3][r16 * u_dim_size * u_dim_size +
                           r17 * u_dim_size +
                           r18] += r19;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r16, r17, r18, r19, m[r16 * max_index * max_index + r17 * max_index + r18]);
        }

        if (ii_src_1 >= -1 && ii_src_2 >= -1 && ii_src_5 >= -1 && ii_src_1 <= max_index && ii_src_2 <= max_index && ii_src_5 <= max_index)
        {
            int r20 = ii_src_5 + 2;
            int r21 = ii_src_1 + 2;
            int r22 = ii_src_2 + 2;
            float r23 = 2.2801e-2F * (1.0F * px * py * pz - 1.0F * px * py - 1.0F * px * pz + 1.0F * px) * src[t] / m[r20 * max_index * max_index + r21 * max_index + r22];
            u[(t + 1) % 3][r20 * u_dim_size * u_dim_size +
                           r21 * u_dim_size +
                           r22] += r23;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r20, r21, r22, r23, m[r20 * max_index * max_index + r21 * max_index + r22]);
        }

        if (ii_src_1 >= -1 && ii_src_3 >= -1 && ii_src_5 >= -1 && ii_src_1 <= max_index && ii_src_3 <= max_index && ii_src_5 <= max_index)
        {
            int r24 = ii_src_5 + 2;
            int r25 = ii_src_1 + 2;
            int r26 = ii_src_3 + 2;
            float r27 = 2.2801e-2F * (-1.0F * px * py * pz + 1.0F * px * pz) * src[t] / m[r24 * max_index * max_index + r25 * max_index + r26];
            u[(t + 1) % 3][r24 * u_dim_size * u_dim_size +
                           r25 * u_dim_size +
                           r26] += r27;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r24, r25, r26, r27, m[r24 * max_index * max_index + r25 * max_index + r26]);
        }

        if (ii_src_2 >= -1 && ii_src_4 >= -1 && ii_src_5 >= -1 && ii_src_2 <= max_index && ii_src_4 <= max_index && ii_src_5 <= max_index)
        {
            int r28 = ii_src_5 + 2;
            int r29 = ii_src_4 + 2;
            int r30 = ii_src_2 + 2;
            float r31 = 2.2801e-2F * (-1.0F * px * py * pz + 1.0F * px * py) * src[t] / m[r28 * max_index * max_index + r29 * max_index + r30];
            u[(t + 1) % 3][r28 * u_dim_size * u_dim_size +
                           r29 * u_dim_size +
                           r30] += r31;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r28, r29, r30, r31, m[r28 * max_index * max_index + r29 * max_index + r30]);
        }

        if (ii_src_3 >= -1 && ii_src_4 >= -1 && ii_src_5 >= -1 && ii_src_3 <= max_index && ii_src_4 <= max_index && ii_src_5 <= max_index)
        {
            int r32 = ii_src_5 + 2;
            int r33 = ii_src_4 + 2;
            int r34 = ii_src_3 + 2;
            float r35 = 2.2801e-2F * px * py * pz * src[t] / m[r32 * max_index * max_index + r33 * max_index + r34];
            u[(t + 1) % 3][r32 * u_dim_size * u_dim_size +
                           r33 * u_dim_size +
                           r34] += r35;
            // ops_printf("u[%d][%d][%d]=%f | m=%f\n", r32, r33, r34, r35, m[r32 * max_index * max_index + r33 * max_index + r34]);
        }

        ops_dat_release_raw_data(dat_ut[(t + 1) % 3], 0, OPS_WRITE);
        t++;

        // ops_printf(title, "output/u_ops%03d.txt", t);

        // if (t % 250 == 0)
            // ops_print_dat_to_txtfile(dat_ut[t % 3], title);

    } while (t < T_intervals);

    ops_printf("\n---------------------------\n");
    ops_print_dat_to_txtfile(dat_ut[t % 3], "title.txt");

    ops_exit();

    // free(u);
    // free(m);
    // free(damp);

    return 0;
}

void initialize_velocity_model(float *m, int x_size, int y_size, int z_size)
{
    // Initialize velocoty model. 2 layered model.
    for (int z = 0; z < z_size; z++)
    {
        for (int j = 0; j < y_size; j++)
        {
            for (int i = 0; i < x_size; i++)
            {

                m[i + j * x_size + z * x_size * y_size] = 0.25;
                // if (j < (y_size / 2))
                // {
                //     m[i + j * x_size + z * x_size * y_size] = 1 / (1.5 * 1.5);
                // }
                // else
                // {
                //     m[i + j * x_size + z * x_size * y_size] = 1 / (2.5 * 2.5);
                // }
            }
        }
    }
}

void initialize_damp(float *damp, int x_size, int y_size, int z_size)
{
    float dampcoeff;
    float pos;
    float val;

    dampcoeff = 0.2590408229618301; // 1.5 * log(1.0 / 0.001) / (40.)

    for (int dim = 0; dim < 3; dim++)
    {
        for (int border = 0; border < border_size; border++)
        {
            pos = abs((border_size - border + 1) / (float)border_size);
            val = dampcoeff * (pos - sin(2 * M_PI * pos) / (2 * M_PI));

            for (int k = 1; k < z_size - 1; k++)
            {
                for (int j = 1; j < y_size - 1; j++)
                {
                    for (int i = 1; i < x_size - 1; i++)
                    {
                        // Left slice
                        if ((dim == 0 && i == border + 1) ||
                            (dim == 1 && j == border + 1) ||
                            (dim == 2 && k == border + 1))
                        {
                            damp[i + j * x_size + k * x_size * y_size] += val;
                        }

                        // Right slice
                        if ((dim == 0 && i == x_size - border - 1) ||
                            (dim == 1 && j == y_size - border - 1) ||
                            (dim == 2 && k == z_size - border - 1))
                        {
                            damp[i + j * x_size + k * x_size * y_size] += val;
                        }
                    }
                }
            }
        }
    }

    // Pad edge
    int new_i, new_j, new_k;
    for (int k = 0; k < z_size; k++)
    {
        for (int j = 0; j < y_size; j++)
        {
            for (int i = 0; i < x_size; i++)
            {
                if ((i < 1 || i > X_size + 2 * border_size) ||
                    (j < 1 || j > Y_size + 2 * border_size) ||
                    (k < 1 || k > Z_size + 2 * border_size))
                {
                    new_i = i < 1 ? i + 1 : i > x_size - 2 ? i - 1 : i;
                    new_j = j < 1 ? j + 1 : j > y_size - 2 ? j - 1 : j;
                    new_k = k < 1 ? k + 1 : k > z_size - 2 ? k - 1 : k;

                    // ops_printf("(%d %d %d) -> (%d %d %d)\n", new_i, new_j, new_k, i, j, k);
                    damp[i + j * x_size + k * x_size * y_size] = damp[new_i + new_j * x_size + new_k * x_size * y_size];
                }
            }
        }
    }
}

void initialize_source(float *src, int total_time)
{
    float f0 = 0.1; // Peak frequency

    ops_printf("Total time: %d\n", total_time);
    for (int t = 0; t < total_time; t++)
    {
        float r = (M_PI * f0 * (dt * t - 1. / f0));
        src[t] = (1 - 2. * pow(r, 2)) * exp(-pow(r, 2));
        // ops_printf("src[%d] = %f\n", t, src[t]);
    }
}

void print_vector(float *vec, int x_limit, int y_limit, int z_limit)
{
    FILE * pFile;
    pFile = fopen ("velocity1.txt","w");

    fprintf(pFile, "Vector of dimensions: %d %d %d\n", x_limit, y_limit, z_limit);

    for (int k = 0; k < (z_limit == 0 ? 1 : z_limit); k++)
    {
        for (int j = 0; j < (y_limit == 0 ? 1 : y_limit); j++)
        {
            for (int i = 0; i < (x_limit == 0 ? 1 : x_limit); i++)
            {
                fprintf(pFile, "%.6e ", vec[i + j * x_limit + k * x_limit * y_limit]);
            }
            fprintf(pFile, "\n");
        }
        fprintf(pFile, "\n");
    }

    fprintf(pFile, "\n\n");
    fclose (pFile);
}
