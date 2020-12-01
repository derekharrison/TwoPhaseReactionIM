/*
 * user_types.hpp
 *
 *  Created on: Dec 1, 2020
 *      Author: d-w-h
 */

#ifndef USER_TYPES_HPP_
#define USER_TYPES_HPP_

typedef struct grid_params {
    int N;
    double H;
    double W;
    double L;
} g_params;

typedef struct physical_params {
    double Ug;
    double Ul;
    double kg;
    double kl;
    double K;
    double Cag0;
    double Cal0;
} p_params;

typedef struct time_data {
    int Nt;
    double to;
    double tf;
} t_data;

typedef struct solver_data {
    double* Cag;
    double Cal;
    double* x_c;
} s_data;

#endif /* USER_TYPES_HPP_ */
