/*
 * main.cpp
 *
 *  Created on: Nov 25, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "solver.hpp"
#include "user_types.hpp"

double ra(double Ca) {
    /* Reaction rate law */
    double kr = 1.0;
    return -kr*Ca*Ca/(1+Ca);
}

int main(int argc, char* argv[]) {
    g_params grid_params;
    p_params physical_params;
    t_data time_data;
    s_data solver_data;

    /* Simulation parameters */
    grid_params.N = 10;                    //Number of nodes
    grid_params.H = 0.5;                   //Height of liquid/gas section
    grid_params.W = 1.0;                   //Width of reactor
    grid_params.L = 1.0;                   //Length of reactor

    physical_params.Ug = 1.0;              //Gas phase velocity
    physical_params.Ul = 1.0;              //Liquid phase velocity
    physical_params.kg = 1.0;              //Mass transfer coefficient gas phase
    physical_params.kl = 1.0;              //Mass transfer coefficient liquid phase
    physical_params.K = 2.0;               //Equilibrium coefficient
    physical_params.Cag0 = 1.0;            //Inlet concentration of component A in gas phase
    physical_params.Cal0 = 0.0;            //Inlet concentration of component A in liquid phase

    time_data.Nt = 20;                     //Number of timesteps
    time_data.to = 0.0;                    //Initial time
    time_data.tf = 1.0;                    //Final time

    /* Allocate memory for solver results */
    solver_data.Cag = new double[grid_params.N];
    solver_data.x_c = new double[grid_params.N];

    /* Execute solver */
    solver(grid_params, physical_params, time_data, &solver_data);

    /* Print results */
    for(int i = 0; i < grid_params.N; ++i) {
        printf("x: %f, Cag: %f, Cal: %f\n", solver_data.x_c[i], solver_data.Cag[i], solver_data.Cal);
    }

    return 0;
}
