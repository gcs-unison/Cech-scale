#ifndef CECH_SCALE_H
#define CECH_SCALE_H

#include <tuple> // for std::tuple, std::make_tuple

#include "auxiliary_functions.h"

/*
 * Calculates the cech scale of a disks system.
 * This function corresponds to algorithm 5 of the paper:
 * "A numerical approach for the filtered generalized cech complex".
 * @param disk_system System of disks.
 *
 * @return A tuple of the results which are, in order: the cech scale of the,
 *  system, the vietoris-rips scale of the system and the intesrection point of
 *  the system.
 */
std::tuple<double, double, std::vector<double>>
calculate_cech_scale(std::vector< std::vector<double> > disk_system);

#endif //CECH_SCALE_H

