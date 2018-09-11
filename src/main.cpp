/* The following script corresponds to the algororithm Cech.scale,
 * which was established in the paper "A numerical approach for the filtered generalized Cech complex",
 * by Jesus F. Espinoza, Rosalia Hernandez-Amador, Hector Alfredo Hernandez Hernandez and Beatriz Ramonetti Valencia.
 */

/* Cech.scale is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation. This program is distributed without any warranty. See the GNU General Public License for more details.
 */

#include "cech_scale.h"


int main(int argc, char *argv[])
{
    std::string input_file = "textfiles/Disks-system.txt";
    std::string output_file = "textfiles/new_cech_scale.txt";
    if(argc >= 2)
        input_file = std::string(argv[1]);
    if(argc >= 3)
        output_file = std::string(argv[2]);

    calculate_cech_scale(input_file, output_file);

    return 0;
}

