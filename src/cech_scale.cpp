#include "cech_scale.h"

//#############################################################################

bool cech_scale(std::string input_file /*= ""*/, std::string output_file /*= ""*/)
{
    //read and validate disk system from text file
    std::cout << "Reading from file: " << input_file << std::endl;

    std::vector< std::vector<double> > disk_system;
    if(!read_file(disk_system, input_file)){
        return false;
    }

    //calculate max vietori rips scale of the disk system
    double vietori_rips_system = max_vietori_rips(disk_system);
    double cech_scale_val = vietori_rips_system;


    //verify if rho is nonnegative
    if(!rho_nonnegative(disk_system, vietori_rips_system)){
        unsigned number_disks = disk_system.size();
        for(unsigned i = 0; i < number_disks - 2; ++i){
            for(unsigned j = i + 1; j < number_disks - 1; ++j){
                for(unsigned k = j + 1; k < number_disks; ++k){

                    //calculate vietori rips of disk1, disk2, disk3
                    double vietori_rips_d123 = max_vietori_rips(std::vector< std::vector<double> >({disk_system[i],
                                                                                                    disk_system[j],
                                                                                                    disk_system[k]}));

                    if(sqrt(4.0/3.0)*vietori_rips_d123 >= cech_scale_val){
                        cech_scale_val = std::max(cech_scale_val, vietori_rips_d123);
                    }else{
                        cech_scale_val = std::max(cech_scale_val,
                                     bisection(disk_system, vietori_rips_d123, sqrt(4.0/3.0)*vietori_rips_d123, 12));
                    }
                }
            }
        }
    }

    write_file(cech_scale_val, vietori_rips_system, std::vector<double>({0,0}), output_file);
    std::cout << "Wrote results to file: " << output_file << std::endl;

    return true;
    //validation: if 2d then apply all above steps with any number of disks.
    //            if Rd, R>2 then must be only 3 disks and must apply a projection plane using the 3 centers of the disks, apply above steps and redimention the intersection points to Rd.
    //
    //output cech scale, point of intersection and wheter cech and vieti rips scale coincide
}

//#############################################################################

