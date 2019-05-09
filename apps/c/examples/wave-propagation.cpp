#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define OPS_3D
#include "ops_seq.h"
int dimensions_number = 3;
int X_size = 1024;
int Y_size = 1024;
int Z_size = 1024;
int padding = 1;
double dt = 0.001, start = 0, stop = 30; // time variables
int border_size = 10;                    // Abosrbent border
int space_order = 2;                     // Space order
int ii_src[6], u_dim_size;
float p[3];

#include "wave-propagation-kernels.h"

void initialize_velocity_model(float *m, int x_size, int y_size, int z_size);
void initialize_damp(float *damp, int x_size, int y_size, int z_size);
void initialize_source(float *src, int T_intervals);
void print_vector(float *vec, int x_limit, int y_limit, int z_limit);
void initialize_source_coordinates(float* src_coords);
void calculate_source_interpolation_position(float* src_coords, float* p, int* ii_src);

int main(int argc, char *argv[])
{
    int t;
    int size[3];
    int damp_size[3];
    int base[] = {0, 0, 0};
    int d_m[] = {-padding, -padding, -padding};
    int d_p[] = {padding, padding, padding};
    int d_m_0[] = {0, 0, 0};
    int d_p_0[] = {0, 0, 0};
    float **u, *m, *damp, *src, src_coords[3];
    int T_intervals;
    char title[25];
    double ct0, ct1, et0, et1; 
    X_size -= 2*border_size + 2*space_order;
    Y_size -= 2*border_size + 2*space_order;
    Z_size -= 2*border_size + 2*space_order;

    ops_printf("Size = %d\n", X_size);

    // Esse size leva em consideracao a borda absorvente???
    size[0] = X_size + 2 * border_size + 2 * space_order;
    size[1] = Y_size + 2 * border_size + 2 * space_order;
    size[2] = Z_size + 2 * border_size + 2 * space_order;
    damp_size[0] = X_size + 2 * border_size + space_order;
    damp_size[1] = Y_size + 2 * border_size + space_order;
    damp_size[2] = Z_size + 2 * border_size + space_order;
    T_intervals = ceil((stop - start + dt) / dt);
    ops_printf("T_intervals = %d\n", T_intervals);

    ops_init(argc, NULL, 1);

    // Alocates and initialize grid
    u = (float **)malloc(3 * sizeof(float *));
    u[0] = (float *)malloc((X_size + 2 * border_size + 2 * space_order) * (Y_size + 2 * border_size + 2 * space_order) * (Z_size + 2 * border_size + 2 * space_order) * sizeof(float));
    u[1] = (float *)malloc((X_size + 2 * border_size + 2 * space_order) * (Y_size + 2 * border_size + 2 * space_order) * (Z_size + 2 * border_size + 2 * space_order) * sizeof(float));
    u[2] = (float *)malloc((X_size + 2 * border_size + 2 * space_order) * (Y_size + 2 * border_size + 2 * space_order) * (Z_size + 2 * border_size + 2 * space_order) * sizeof(float));
    m = (float *)malloc((X_size + 2 * border_size + 2 * space_order) * (Y_size + 2 * border_size + 2 * space_order) * (Z_size + 2 * border_size + 2 * space_order) * sizeof(float));
    damp = (float *)malloc((X_size + 2 * border_size + space_order) * (Y_size + 2 * border_size + space_order) * (Z_size + 2 * border_size + space_order) * sizeof(float));
    src = (float *)malloc(T_intervals * sizeof(float));

    // Initialize velocity model
    // initialize_velocity_model(m, X_size + 2 * border_size + 2 * space_order, Y_size + 2 * border_size + 2 * space_order, Z_size + 2 * border_size + 2 * space_order);
    // print_vector(m, X_size + 2 * border_size + 2 * space_order, Y_size + 2 * border_size + 2 * space_order, Z_size + 2 * border_size + 2 * space_order);

    // Initialize Damp
    // initialize_damp(damp, X_size + 2 * border_size + 1 + 1, Y_size + 2 * border_size + 1 + 1, Z_size + 2 * border_size + 1 + 1);
    // print_vector(damp, X_size + 2 * border_size + 1 + 1, Y_size + 2 * border_size + 1 + 1, Z_size + 2 * border_size + 1 + 1);

    // Initialize source
    initialize_source(src, T_intervals);
    // print_vector(src, T_intervals, 0, 0);
    // As we have only one source, this is simpler. For multiple sources, this has to be changed later.
    calculate_source_interpolation_position(src_coords, p, ii_src);
    // printf("ii_src: %d %d %d %d %d %d\n", ii_src[0], ii_src[1], ii_src[2], ii_src[3], ii_src[4], ii_src[5]);

    // Initialize source position
    initialize_source_coordinates(src_coords);

    // Declare global constant
    ops_decl_const("X_size", 1, "int", &X_size);
    ops_decl_const("Y_size", 1, "int", &Y_size);
    ops_decl_const("Z_size", 1, "int", &Z_size);
    ops_decl_const("border_size", 1, "int", &border_size);
    ops_decl_const("dt", 1, "double", &dt);
    ops_decl_const("space_order", 1, "int", &space_order);
    ops_decl_const("u_dim_size", 1, "int", &u_dim_size);
    ops_decl_const("ii_src", 6, "int", ii_src);
    ops_decl_const("p", 3, "float", p);

    ops_printf("Src coordinates:%f %f %f\n", src_coords[0], src_coords[1], src_coords[2]);

    // Declare ops_block
    ops_block grid = ops_decl_block(3, "grid");

    ops_block srcB = ops_decl_block(1, "src");

    // Declare ops_dat objects
    ops_dat dat_ut[3];
    ops_dat dat_m;
    ops_dat dat_damp, dat_damp2;
    ops_dat dat_x_size, dat_y_size, dat_z_size;
    dat_ut[0] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[0], "float", "ut0");
    dat_ut[1] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[1], "float", "ut1");
    dat_ut[2] = ops_decl_dat(grid, 1, size, base, d_m, d_p, u[2], "float", "ut2");
    dat_m = ops_decl_dat(grid, 1, size, base, d_m, d_p, m, "float", "m");
    dat_damp = ops_decl_dat(grid, 1, damp_size, base, d_m, d_p, damp, "float", "damp");
    dat_damp2 = ops_decl_dat(grid, 1, damp_size, base, d_m, d_p, damp, "float", "damp");

    int s3d_000[] = {0, 0, 0};
    int s3d_1pt[] = {-1, -1, -1};
    int s3d_7pts[] = {0, 0, 0,
                      -1, 0, 0,
                      0, -1, 0,
                      0, 0, -1,
                      0, 0, 1,
                      0, 1, 0,
                      1, 0, 0};
    int s3d_15pts[] = {0, 0, 0,
                       0, 0, 1,
                       0, 1, 0,
                       0, 1, 1,
                       1, 0, 0,
                       1, 0, 1,
                       1, 1, 0,
                       1, 1, 1,
                       0, 0, -1,
                       0, -1, 0,
                       0, -1, -1,
                      -1, 0, 0,
                      -1, 0, -1,
                      -1, -1, 0,
                      -1, -1, -1};

    // Declare stencil
    ops_stencil S3D_000 = ops_decl_stencil(3, 1, s3d_000, "0,0,0");
    ops_stencil S3D_1PT = ops_decl_stencil(3, 1, s3d_1pt, "-1,-1,-1");
    ops_stencil S3D_7PTS = ops_decl_stencil(3, 7, s3d_7pts, "7pts");
    ops_stencil S3D_15PTS = ops_decl_stencil(3, 15, s3d_15pts, "15pts");
    int max_index = X_size + 2 * border_size + 2 * space_order;
    int u_dim_size = X_size + 2 * border_size;

    int shift = space_order;
    int whole_range[] = {space_order, X_size + 2 * border_size, space_order, Y_size + 2 * border_size, space_order, Z_size + 2 * border_size};
    int velocity_model_range[] = {
        0, X_size + 2 * border_size,
        0, X_size + 2 * border_size,
        0, X_size + 2 * border_size};
    int damp_range[] = {
        0, X_size + 2 * border_size + 1 + 1,
        0, Y_size + 2 * border_size + 1 + 1,
        0, Z_size + 2 * border_size + 1 + 1
        };
    // Define range of injection kernel. That's the interpolation result.
    int injection_range[] = {
        ii_src[2] + 2,
        ii_src[3] + 2 + 1,
        ii_src[1] + 2,
        ii_src[4] + 2 + 1,
        ii_src[0] + 2,
        ii_src[5] + 2 + 1
    };
    // ops_printf("Starts time propagation.\n");

    int disp[3], sizes[3], strides[3];

    // ops_printf("time=%d\n", t);
    // print_vector(u[(t + 1) % 3], 34, 34, 34);
    // ops_print_dat_to_txtfile(dat_ut[t % 3], "output/u_ops000.txt");

    ops_partition("");
    ops_diagnostic_output();

    // Initialize velocity model
    ops_par_loop(initialize_velocity_model_kernel, "initialize_velocity_model_kernel", grid, 3, velocity_model_range,
                    ops_arg_dat(dat_m, 1, S3D_000, "float", OPS_WRITE));
    // ops_print_dat_to_txtfile(dat_m, "model.txt");
    //Initialize u1, u2 and u3
    ops_par_loop(set_zero_kernel, "set_zero_kernel", grid, 3, damp_range,
                    ops_arg_dat(dat_ut[0], 1, S3D_000, "float", OPS_WRITE));
    ops_par_loop(set_zero_kernel, "set_zero_kernel", grid, 3, damp_range,
                    ops_arg_dat(dat_ut[1], 1, S3D_000, "float", OPS_WRITE));
    ops_par_loop(set_zero_kernel, "set_zero_kernel", grid, 3, damp_range,
                    ops_arg_dat(dat_ut[2], 1, S3D_000, "float", OPS_WRITE));

    // Initialize Initialize Damp
    ops_par_loop(set_zero_kernel, "set_zero_kernel", grid, 3, damp_range,
                    ops_arg_dat(dat_damp2, 1, S3D_000, "float", OPS_WRITE));
    ops_par_loop(initialize_damp_kernel, "initialize_damp_kernel", grid, 3, damp_range,
                    ops_arg_dat(dat_damp2, 1, S3D_000, "float", OPS_WRITE),
                    ops_arg_idx());

    ops_par_loop(set_space_order_border_kernel, "set_space_order_border_kernel", grid, 3, damp_range,
                ops_arg_dat(dat_damp, 1, S3D_000, "float", OPS_WRITE),
                ops_arg_dat(dat_damp2, 1, S3D_15PTS, "float", OPS_READ),
                ops_arg_idx());
    // ops_print_dat_to_txtfile(dat_damp, "damp_par.txt");
    // damp = (float *)ops_dat_get_raw_pointer(dat_damp, 0, S3D_000, strides);
    // print_vector(damp, X_size + 2 * border_size + 1 + 1, Y_size + 2 * border_size + 1 + 1, Z_size + 2 * border_size + 1 + 1);

    ops_printf("\n***********************************************************\n");
    // Start time
    ops_timers(&ct0, &et0);
    for(t = 1; t < T_intervals; t++)
    {
        // ops_printf("%d\n", t);
        ops_par_loop(wave_propagation_kernel, "wave_propagation_kernel", grid, 3, whole_range,
                     ops_arg_dat(dat_damp, 1, S3D_1PT, "float", OPS_READ),
                     ops_arg_dat(dat_m, 1, S3D_000, "float", OPS_READ),
                     ops_arg_dat(dat_ut[(t + 1) % 3], 1, S3D_000, "float", OPS_WRITE), // t1
                     ops_arg_dat(dat_ut[t % 3], 1, S3D_7PTS, "float", OPS_READ),       // t0
                     ops_arg_dat(dat_ut[(t + 2) % 3], 1, S3D_000, "float", OPS_READ),  // t2
                     ops_arg_idx());

        ops_par_loop(source_injection_kernel, "source_injection_kernel", grid, 3, injection_range,
                     ops_arg_dat(dat_ut[(t + 1) % 3], 1, S3D_000, "float", OPS_WRITE),
                     ops_arg_dat(dat_m, 1, S3D_000, "float", OPS_READ),
                     ops_arg_gbl(&src[t], 1, "float", OPS_READ),
                     ops_arg_idx());
    }
    // End time
    ops_timers(&ct1, &et1);
    ops_printf("\nTotal Wall time %lf\n", et1 - et0);

    ops_printf("\n---------------------------\n");
    // ops_print_dat_to_txtfile(dat_ut[t % 3], "domain/title.txt");

    ops_exit();

    free(u);
    free(m);
    free(damp);
    free(src);

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


/**
 * Initializes the sources coordinates.
 * This example only has one source, so the initialization is simple, but in the
 *future it will handle multiple sources. The current source is intialize in the
 *middle top of the volume.
 **/
void initialize_source_coordinates(float* src_coords) {
    src_coords[0] = (float)(X_size - 1) / (float)2;
    src_coords[1] = (float)(Y_size - 1) / (float)2;
    src_coords[2] = (float)1.;
}



void initialize_source(float* src, int total_time) {
    float f0 = 0.1;  // Peak frequency

    double sigma = 1.0 / (M_PI * f0 * sqrt(2.0));
    double t0 = 6.0 * sigma;

    for (int t = 0; t < total_time; t++) {
        float r = (M_PI * f0 * (t * dt - t0)) * (M_PI * f0 * (t * dt - t0));
        src[t] = (1 - 2.0 * r) * exp(-r);
    }
}

/**
 * @brief  Calculates the source interpolation.
 * @note
 * @param  *src_coords: Source coordinates defined in other system.
 * @param  *p: [output] source center position.
 * @param  *ii_src: [output] Interpolation position.
 * @retval None
 */
void calculate_source_interpolation_position(float* src_coords, float* p, int* ii_src) {
    // source injection interpolation
    float r1 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[0]));
    ii_src[0] = r1;
    float r2 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[1]));
    ii_src[1] = r2;
    float r3 = (int)(floor(-1.0F * -border_size + 1.0F * src_coords[2]));
    ii_src[2] = r3;
    ii_src[3] = r3 + 1;
    ii_src[4] = r2 + 1;
    ii_src[5] = r1 + 1;
    p[0] = (float)(border_size - 1.0F * r1 + src_coords[0]);
    p[1] = (float)(border_size - 1.0F * r2 + src_coords[1]);
    p[2] = (float)(border_size - 1.0F * r3 + src_coords[2]);
}

void print_vector(float *vec, int x_limit, int y_limit, int z_limit)
{
    FILE * pFile;
    pFile = fopen ("damp_sequential.txt", "w");

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
