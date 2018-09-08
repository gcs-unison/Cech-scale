#include "cech_scale.h"


//temporals
#include <iterator>
//temporals


bool rho_nonnegative(const std::vector< std::vector<double> >& disk_system,
                       double lambda_val)
{
    for(std::vector<double> disk1 : disk_system){
        for(std::vector<double> disk2 : disk_system){
            if(disk1 == disk2)
                continue;

            bool nonnegative = true;
            for(std::vector<double> disk3: disk_system){
                if(disk3 == disk1 || disk3 == disk2)
                    continue;

                if(lambda(disk1, disk2, disk3, lambda_val) < 0){
                    nonnegative = false;
                }
            }
            if(nonnegative)
                return true;

        }
    }

    return false;
}

//#############################################################################

bool cech_scale(std::string input_file /*= ""*/, std::string output_file /*= ""*/)
{
    //steps for the algorithm
    //read and validate disk system from text file

    std::cout << "Reading from file: " << input_file << std::endl;

    std::vector< std::vector<double> > disk_system;
    if(!read_file(disk_system, input_file)){
        return false;
    }

    //calculate max vietori rips scale of the disk system
    double vietori_rips_system = max_vietori_rips(disk_system);
    double cech_scale_val = vietori_rips_system;


    //verify if rho is positive (if the minimum output is positive)
    if(!rho_nonnegative(disk_system, vietori_rips_system)){
        unsigned number_disks = disk_system.size();
        for(unsigned i = 0; i < number_disks - 2; ++i){
            for(unsigned j = i + 1; j < number_disks - 1; ++j){
                for(unsigned k = j + 1; k < number_disks; ++k){

                    //calculate vietori rips of d1, d2, d3
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
    //  if so, then output the vietori rips scale and point of intersection and end
    //if is not completely positive, apply bisection to a 3-disk system rho to approximate cech scale. Repeat for all possible 3-disk systems
    //take the max cech scale calculated
    //
    //validation: if 2d then apply all above steps with any number of disks.
    //            if Rd, R>2 then must be only 3 disks and must apply a projection plane using the 3 centers of the disks, apply above steps and redimention the intersection points to Rd.
    //
    //output cech scale, point of intersection and wheter cech and vieti rips scale coincide
}

