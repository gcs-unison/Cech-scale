#include "cech_scale.h"

//#############################################################################

bool calculate_cech_scale(std::string input_file /*= ""*/, std::string output_file /*= ""*/)
{
    //read and validate disk system from text file
    std::cout << "Reading from file: " << input_file << std::endl;
    std::vector< std::vector<double> > read_system;
    if(!read_file(read_system, input_file)){
        return false;
    }
    int dimentions = read_system[0].size() - 1;

    //transforms the system from more than 3 dimentions to only 2 and 3 disks or
    //just makes a copy of the disk system with 2 dimentions
    std::vector< std::vector<double> > disk_system;
    if(dimentions >= 3){
        disk_system = transform_disk_system(read_system);
    }else{
        disk_system = read_system;
    }

    double vietori_rips_system;
    std::vector<double> intersection;
    //calculate max vietori rips scale of the disk system
    std::tie(vietori_rips_system, intersection) = max_vietori_rips_intersection(disk_system);
    double cech_scale = vietori_rips_system;

    //verify if rho is nonnegative
    if(!rho_nonnegative(disk_system, cech_scale)){

        unsigned d1_idx = 0;
        unsigned d2_idx = 1;
        unsigned number_disks = disk_system.size();
        for(unsigned i = 0; i < number_disks - 2; ++i){
            for(unsigned j = i + 1; j < number_disks - 1; ++j){
                for(unsigned k = j + 1; k < number_disks; ++k){

                    //calculate vietori rips of disk_i, disk_j, disk_k
                    std::vector<std::vector<double>> disk_systemijk = {disk_system[i],
                                                                       disk_system[j],
                                                                       disk_system[k]};
                    double vietori_rips_dijk = max_vietori_rips(disk_systemijk);


                    if(sqrt(4.0/3.0)*vietori_rips_dijk >= cech_scale){
                        if(rho(disk_systemijk, vietori_rips_dijk) >= 0){
                            cech_scale = vietori_rips_dijk;
                        }else{
                            double cech_bisection = bisection(disk_systemijk,
                                                              vietori_rips_dijk,
                                                              sqrt(4.0/3.0)*vietori_rips_dijk,
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

        if(dimentions < 3)
            intersection = left_intersection_scaled(disk_system[d1_idx], disk_system[d2_idx], cech_scale);
        else
            intersection = right_intersection_scaled(disk_system[d1_idx], disk_system[d2_idx], cech_scale);
    }

    if(dimentions >= 3){
        //transform the intersecion back to a greater dimention
        intersection = transform_intersection(disk_system, read_system, intersection);
    }

    write_file(cech_scale, vietori_rips_system, intersection, output_file);
    std::cout << "Wrote results to file: " << output_file << std::endl;

    return true;
}

//#############################################################################

