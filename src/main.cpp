/* The following script corresponds to the algororithm Cech.scale,
 * which was established in the paper "A numerical approach for the filtered generalized Cech complex",
 * by Jesus F. Espinoza, Rosalia Hernandez-Amador, Hector Alfredo Hernandez Hernandez and Beatriz Ramonetti Valencia.
 */

/* Cech.scale is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. This program is distributed without any warranty. See the GNU General Public License for more details.
 */

#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <fstream>
#include <cstring>
#include <limits> //for std::numeric_limits<double>::lowest
#include <vector> //for std::vector
#include <string> //for std::string
#include <algorithm> //for std::min_element

#include "circle_circle_intersection.h"
#include "auxiliary_functions.h"

using namespace std;

#define maxm 1000


//////////////////////////////////////
//       AUXILIARY FUNCTIONS        //
//////////////////////////////////////


//////////////////////////////////////
//////          MAIN            //////
//////////////////////////////////////

int main(int argc, char *argv[]){
    double M[maxm][3];
    int m, d;
    int i, j, D, q, k, k_index;
    int ijk_rho_sign, edge_max, double_intersection, ij_rho_nonnegative;
    int Rips_one_simplex[2], two_simplex[2], k_values[m-2];

    double Cech_scale, Rips_scale, Cech_scale_sqrt2, ijk_cech_scale, ij_scale, two_simplex_cech_scale;
    double dist_ij, r_ij, a, b, distdkcenter;
    double sum, sum21,sum22;
    double i_radious_reescaled, j_radious_reescaled, k_radious_reescaled, k_radious_reescaled2, radious1_reescaled, radious2_reescaled;

    double edges_scale[maxm][maxm];
    double i_center[2], j_center[2], k_center[2];
    double c_ij[2], d_ij[2], d_ji[2], n_ij[2];
    double D1[3],D2[3],D3[3];
    double center1[2],center2[2];
    double c_star[2];
    double M_ijk[3][3];
    double aa, ab, ac, ad;

    std::string filename = "textfiles/Disks-system.txt";

    if(argc > 1)
        filename = std::string(argv[1]);

    std::cout << "Reading from file: " << filename << std::endl;

    double p_M[d];
    std::vector< std::vector<double> > N;

    if(!read_file(m, d, M, N, filename) )
        return 1;


    Rips_scale=0;
    for(i=0;i<m-1;i++){
        for(j=i+1;j<m;j++){
            dist_ij=0; for(D=0; D < 2; D++) dist_ij+=pow(M[i][D]-M[j][D],2); dist_ij=sqrt(dist_ij);
            edges_scale[i][j] = dist_ij/(M[i][2]+M[j][2]);
            if(edges_scale[i][j] > Rips_scale){
                Rips_scale = edges_scale[i][j];
                Rips_one_simplex[0] = i;
                Rips_one_simplex[1] = j;
            }
        }
    }

    // Verify if rho_M(rips_scale) >= 0.
    i = Rips_one_simplex[0]; j = Rips_one_simplex[1];

    i_radious_reescaled = Rips_scale*M[i][2]; j_radious_reescaled = Rips_scale*M[j][2];
    for(D=0; D < 2; D++) i_center[D] = M[i][D]; for(D=0; D < 2; D++) j_center[D] = M[j][D];

    for(D=0; D < 2; D++){ d_ij[D] = (j_radious_reescaled*i_center[D] + i_radious_reescaled*j_center[D])/(i_radious_reescaled + j_radious_reescaled);}

    int count = 0;
    for(D=0; D<m; D++)
        if(D!=i && D!=j)
            k_values[count++] = D;

    // Suppose nonnegativity:
    ij_rho_nonnegative = 1, k_index = 0;
    while(ij_rho_nonnegative == 1 && k_index != m-2){
        k = k_values[k_index];
        for(D=0; D < 2; D++) k_center[D]=M[k][D];

        k_radious_reescaled2 = pow(Rips_scale*M[k][2],2);
        distdkcenter=0; for(D=0; D < 2; D++) distdkcenter+=pow(d_ij[D] - k_center[D],2);
        if(k_radious_reescaled2 < distdkcenter){
            ij_rho_nonnegative = 0;
        }
        k_index++;
    }

    // If rho_M(Rips_scale) >=0, then the Rips scale agree with the Cech scale
    if(ij_rho_nonnegative==1){
        ofstream output_file;
        output_file.open("Cech_scale.txt", ios::trunc );
        output_file << "The Cech scale agrees with the Rips scale.\n" << endl;
        output_file << "Cech scale:\n" << scientific << Rips_scale << endl;
        output_file << "\nIntersection point of the rescaled disks system:\n";
        output_file << scientific << d_ij[0] << endl;
        output_file << scientific << d_ij[1] << endl;
        output_file.close();

//         cout << scientific << Rips_scale << endl;
//         cout << scientific << d_ij[0] << endl;
//         cout << scientific << d_ij[1] << endl;
        return 0;
    }


    // The following is under the value: verify_nonneg == 1

    Cech_scale = Rips_scale;
    // We propose as Cech scale the value of the Rips scale.
    // We will update the Cech scale in each iteration:

    // Check every triplet [i,j,k]:
    for(i=0;i<m-2;i++){
        for(j=i+1;j<m-1;j++){
            ij_scale = edges_scale[i][j];

            for(k=j+1;k<m;k++){
                // Calculate the Rips scale of the 2-simplex[i,j,k] and assume it agree with the Cech scale: ijk_cech_scale

                edge_max=3;
                if(ij_scale>=edges_scale[i][k] && ij_scale>=edges_scale[j][k] ) edge_max=1;
                if(edges_scale[i][k]>=ij_scale && edges_scale[i][k]>=edges_scale[j][k]) edge_max=2;

                if(edge_max == 1){
                    for(D=0;D<3;D++){
                        D1[D]=M[i][D];
                        D2[D]=M[j][D];
                        D3[D]=M[k][D];
                    }
                    ijk_cech_scale = ij_scale;
                }else{
                    if(edge_max == 2){
                        for(D=0;D<3;D++){
                            D1[D]=M[i][D];
                            D2[D]=M[k][D];
                            D3[D]=M[j][D];
                        }
                        ijk_cech_scale = edges_scale[i][k];
                    }else{
                        for(D=0;D<3;D++){
                            D1[D]=M[j][D];
                            D2[D]=M[k][D];
                            D3[D]=M[i][D];
                        }
                        ijk_cech_scale = edges_scale[j][k];
                    }
                }

                // The value of ijk_cech_scale must be at least Cech_scale/sqrt(4/3):
                if(ijk_cech_scale*sqrt(4.0/3) >= Cech_scale){

                    // Verify is the Rips scale of [i,j,k] agree with the Cech scale of [i,j,k].
                    for(D=0; D < 2; D++) center1[D]  = D1[D]; radious1_reescaled = ijk_cech_scale*D1[2];
                    for(D=0; D < 2; D++) center2[D]  = D2[D]; radious2_reescaled = ijk_cech_scale*D2[2];
                    for(D=0; D < 2; D++) d_ij[D] = (radious2_reescaled*center1[D] + radious1_reescaled*center2[D])/(radious1_reescaled + radious2_reescaled);

                    sum = sqrt( pow(D3[0]-d_ij[0], 2) + pow(D3[1]-d_ij[1], 2) );
                    // Confirm if the Rips scale of the 2-simplex [i,j,k] agrees with the Cech scale.
                    // If it does not agree, update the Cech scale using the bisection numerical method.
                    if( ijk_cech_scale*D3[2] < sum){
                        for(q=0;q<3;q++){
                            M_ijk[0][q] = M[i][q];
                            M_ijk[1][q] = M[j][q];
                            M_ijk[2][q] = M[k][q];
                        }
                        ijk_cech_scale = bisection(M_ijk, 3, ijk_cech_scale, ijk_cech_scale*sqrt(4.0/3), 12, 1e-12);
                    }


                    // With the Cech scale of the 2-simplex [i,j,k], update the Cech scale of the disks system.
                    if(ijk_cech_scale >= Cech_scale){
                        Cech_scale = ijk_cech_scale;
                        // If inter_point==1, then save the 2-simplex with the maximal Cech scale,
                        // to compute the intersection point.
                        two_simplex[0]=i; two_simplex[1]=j; two_simplex[2]=k;
                        two_simplex_cech_scale = Cech_scale;
                    }
                }
            }
        }
    }

    i = two_simplex[0]; j = two_simplex[1]; k = two_simplex[2];
    for(D=0; D < 2; D++) i_center[D] = M[i][D]; i_radious_reescaled = two_simplex_cech_scale*M[i][2];
    for(D=0; D < 2; D++) j_center[D] = M[j][D]; j_radious_reescaled = two_simplex_cech_scale*M[j][2];
    for(D=0; D < 2; D++) k_center[D] = M[k][D]; k_radious_reescaled = two_simplex_cech_scale*M[k][2];

    for(D=0; D < 2; D++) c_ij[D] = j_center[D] - i_center[D];
    dist_ij=0; for(D=0; D < 2; D++){ dist_ij+=pow(i_center[D]-j_center[D],2);} dist_ij=sqrt(dist_ij);
    n_ij[0]=-c_ij[1]; n_ij[1]=c_ij[0];

    // Calculate the coefficients "a" and "b" such that
    // d_ij = i_center + a * c_ij + b * n_ij and d_ji = i_center + a * c_ij - b * n_ij.
    a = 0.5*(1 + (i_radious_reescaled - j_radious_reescaled)*(i_radious_reescaled + j_radious_reescaled)/pow(dist_ij,2));
    b = sqrt( pow(i_radious_reescaled,2) - pow(a*dist_ij,2) )/dist_ij;

    for(D=0; D < 2; D++) d_ij[D] = i_center[D] + a * c_ij[D] + b * n_ij[D];
    for(D=0; D < 2; D++) d_ji[D] = i_center[D] + a * c_ij[D] - b * n_ij[D];

    sum21=0; for(D=0; D < 2; D++) sum21+=pow(k_center[D] - d_ij[D],2);
    sum22=0; for(D=0; D < 2; D++) sum22+=pow(k_center[D] - d_ji[D],2);
    if(sum21 < sum22){
        for(D=0; D < 2; D++) p_M[D] = d_ij[D];
    }else{
        for(D=0; D < 2; D++) p_M[D] = d_ji[D];
    }

    ofstream output_file;
    output_file.open("Cech_scale.txt", ios::trunc);
    output_file << "The Cech scale is greater than the Rips scale.\n" << endl;
    output_file << "Cech scale:\n" << scientific << Cech_scale << endl;
    output_file << "\nIntersection point of the rescaled disks system:\n";
    output_file << scientific << p_M[0] << endl;
    output_file << scientific << p_M[1] << endl;
    output_file.close();

    return 0;
}

