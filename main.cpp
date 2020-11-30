/*
 * main.cpp
 *
 *  Created on: Nov 25, 2020
 *      Author: d-w-h
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

double ra(double Ca) {
    double kr = 1.0;
    return -kr*Ca*Ca/(1+Ca);
}

double dra_dCa(double Ca) {
    double dCa = 0.00001;
    return (ra(Ca + dCa) - ra(Ca)) / dCa;
}

int main(int argc, char* argv[]) {

    double H, W, L, Ug, Ul, kg, kl, K, Cag0, Cal0;
    double to, tf;
    int N, Nt;

    /* Parameters */
    N = 10;                      //Number of nodes
    Nt = 20;                     //Number of timesteps
    H = 0.5;                     //Height of liquid/gas section
    W = 1.0;                     //Width of reactor
    L = 1.0;                     //Length of reactor
    Ug = 1.0;                    //Gas phase velocity
    Ul = 1.0;                    //Liquid phase velocity
    kg = 1.0;                    //Mass transfer coefficient gas phase
    kl = 1.0;                    //Mass transfer coefficient liquid phase
    K = 2.0;                     //Equilibrium coefficient
    Cag0 = 1.0;                  //Inlet concentration of component A in gas phase
    Cal0 = 0.0;                  //Inlet concentration of component A in liquid phase

    to = 0.0;                    //Initial time
    tf = 1.0;                    //Final time

    /* Allocate memory for solver results */
    double* Cag = new double[N];
    double* Cago = new double[N];
    double* Cag_prev = new double[N];
    double Cal = 0.0;
    double Calo = 0.0;
    double Cal_prev = 0.0;
    double* x_c = new double[N];

    /* Start calculations */
    double del_x = L / N;
    double del_t = (tf - to) / Nt;
    double B = kg + kl / K;
    double t = to;
    double Calin = Cal0;
    int max_outer_iter = 300;
    int max_inner_iter = 300;

    /* Initialize data */
    for(int i = 0; i < N; ++i) {
        Cag[i] = 0.0;
        Cago[i] = 0.0;
        Cal = 0.0;
        Calo = 0.0;
        Cag_prev[i] = 0.0;
        Cal_prev = 0.0;
        x_c[i] = i*del_x + 0.5*del_x;
    }

    /* Start simulation */
    while(t < tf) {
        /* Start Gauss-Seidel iterations */
        int outer_iter = 0;
        while(outer_iter < max_outer_iter) {
            int inner_iter = 0;
            for(int i = 0; i < N; ++i) {
                Cag_prev[i] = Cag[i];
                Cal_prev = Cal;
            }
            /* Inner iterations */
            while(inner_iter < max_inner_iter) {
                /*First node gas phase*/
                Cag[0] = (Ug*H*W*Cag0 + kg*W*del_x*kl*Cal/B + H*W*del_x*Cago[0]/del_t) /
                         (Ug*H*W + kg*W*del_x + H*W*del_x/del_t - kg*W*del_x*kg/(B+1e-20));

                /* 2nd node to last node gas phase */
                for(int i = 1; i < N; ++i) {
                    Cag[i] = (Ug*H*W*Cag[i-1] + kg*W*del_x*kl*Cal/B + H*W*del_x*Cago[i]/del_t) /
                             (Ug*H*W + kg*W*del_x + H*W*del_x/del_t - kg*W*del_x*kg/(B+1e-20));
                }

                /* Liquid phase node */
                double klWdxkgCag_KB = 0.0;
                for(int i = 0; i < N; ++i) {
                    klWdxkgCag_KB += kl*W*del_x*kg*Cag[i]/(K*B);
                }
                Cal = (Ul*W*H*Calin + klWdxkgCag_KB + ra(Cal_prev)*W*H*L - dra_dCa(Cal_prev)*Cal_prev*W*H*L + W*H*L*Calo/del_t) /
                      (Ul*W*H + N*kl*W*del_x + W*H*L/del_t - N*kl*W*del_x*kl/(K*B) - dra_dCa(Cal_prev)*W*H*L);

                inner_iter++;
            }

            outer_iter++;
        }

        /* Set 'old' timestep values */
        for(int i = 0; i < N; ++i) {
            Cago[i] = Cag[i];
            Calo = Cal;
        }

        t += del_t;
    }

    /* Print results */
    for(int i = 0; i < N; ++i) {
        printf("x: %f, Cag: %f, Cal: %f\n", x_c[i], Cag[i], Cal);
    }

    printf("done\n");

    return 0;
}
