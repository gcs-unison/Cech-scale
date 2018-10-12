/* The following script corresponds to the algororithm Cech.scale,
 * which was established in the paper "A numerical approach for the filtered generalized Cech complex",
 * by Jesus F. Espinoza, Rosalia Hernandez-Amador, Hector Alfredo Hernandez Hernandez and Beatriz Ramonetti Valencia.
 *
 * The code was rewriten by Luis Sotomayor.
 */

/* Cech.scale is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. This program is distributed without any warranty. See the GNU General Public License for more details.
 */

#include <iostream> //for std::cout, std::endl
#include <tuple> // for std::tie

#include "cech_scale.h"
#include "auxiliary_functions.h"


int main(int argc, char *argv[])
{
    std::string input_file = "textfiles/disks_system.txt";
    std::string output_file = "textfiles/cech_results.txt";
    if(argc >= 2)
        input_file = std::string(argv[1]);
    if(argc >= 3)
        output_file = std::string(argv[2]);

    //read and validate disk system from text file
    std::cout << "Reading from file: " << input_file << std::endl;
    std::vector< std::vector<double> > disk_system;
    if(!read_file(disk_system, input_file)){
        return 1;
    }

    double cech_scale;
    double vietori_rips;
    std::vector<double> intersection;
    std::tie(cech_scale, vietori_rips, intersection) = calculate_cech_scale(disk_system);

    write_file(cech_scale, vietori_rips, intersection, output_file);
    std::cout << "Wrote results to file: " << output_file << std::endl;

    return 0;
}

