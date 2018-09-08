#ifndef CECH_SCALE_H
#define CECH_SCALE_H

#include "auxiliary_functions.h"

/**
 * Checks whether for all disks in the disk system and the lambda value
 * specified, the rho function always returns a positive number or not.
 * @param disk_system System of disks in R^n.
 * @param lambda_val Value needed to evaluate rho.
 *
 * @return True if for all disks in the disk system, evaluating rho using
 *          the lambda value return a non negative number. False otherwise.
 */
bool rho_nonnegative(const std::vector< std::vector<double> >& disk_system,
                       double lambda_val);


bool cech_scale(std::string input_file = "", std::string output_file = "");

#endif //CECH_SCALE_H

