#ifndef CECH_SCALE_H
#define CECH_SCALE_H

#include <iostream> //for std::cout, std::endl

#include "auxiliary_functions.h"

typedef VietoriRipsScale CechScale;

bool calculate_cech_scale(std::string input_file = "", std::string output_file = "");

#endif //CECH_SCALE_H

