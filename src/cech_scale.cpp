#include "cech_scale.h"

//#############################################################################

std::tuple<double, double, std::vector<double>>
calculate_cech_scale(std::vector< std::vector<double> > input_system)
{
    int dimensions = input_system[0].size() - 1;

    //transforms the system from more than 3 dimensions to only 2 and 3 disks or
    //just makes a copy of the disk system with 2 dimensions
    std::vector< std::vector<double> > disk_system;
    if(dimensions >= 3){
        disk_system = transform_disk_system(input_system);
    }else{
        disk_system = input_system;
    }

    double vietoris_rips_system;
    std::vector<double> intersection;
    //calculate max vietoris rips scale of the disk system
    std::tie(vietoris_rips_system, intersection) = max_vietoris_rips_intersection(disk_system);
    double cech_scale = vietoris_rips_system;

    //verify if rho is nonnegative
    if(!rho_nonnegative(disk_system, cech_scale)){

        unsigned d1_idx = 0;
        unsigned d2_idx = 1;
        unsigned number_disks = disk_system.size();
        for(unsigned i = 0; i < number_disks - 2; ++i){
            for(unsigned j = i + 1; j < number_disks - 1; ++j){
                for(unsigned k = j + 1; k < number_disks; ++k){

                    //calculate vietoris rips of disk_i, disk_j, disk_k
                    std::vector<std::vector<double>> disk_systemijk = {disk_system[i],
                                                                       disk_system[j],
                                                                       disk_system[k]};
                    double vietoris_rips_dijk = max_vietoris_rips(disk_systemijk);


                    if(sqrt(4.0/3.0)*vietoris_rips_dijk >= cech_scale){
                        if(rho(disk_systemijk, vietoris_rips_dijk) >= 0){
                            cech_scale = vietoris_rips_dijk;
                        }else{
                            double cech_bisection = bisection(disk_systemijk,
                                                              vietoris_rips_dijk,
                                                              sqrt(4.0/3.0)*vietoris_rips_dijk,
                                                              12);
                            if(cech_bisection > cech_scale){
                                cech_scale = cech_bisection;
                            }
                        }

                        d1_idx = i;
                        d2_idx = j;
                    }
                }
            }
        }

        if(dimensions < 3)
            intersection = left_intersection_scaled(disk_system[d1_idx], disk_system[d2_idx], cech_scale);
        else
            intersection = right_intersection_scaled(disk_system[d1_idx], disk_system[d2_idx], cech_scale);
    }

    if(dimensions >= 3){
        //transform the intersecion back to a greater dimension
        intersection = transform_intersection(disk_system, input_system, intersection);
    }

    return std::make_tuple(cech_scale, vietoris_rips_system, intersection);
}

//#############################################################################

