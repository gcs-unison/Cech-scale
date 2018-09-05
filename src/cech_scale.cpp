
double rho_nonnegative(const std::vector< std::vector<double> >& disk_system,
                       double lambda_val)
{
    for(std::vector<double> disk1 : disk_system){
        for(std::vector<double> disk2 : disk_system){
            if(disk1 == disk2)
                continue;

            for(std::vector<double> disk3: disk_system){
                if(lambda(disk1, disk2, disk3, lambda) < 0)
                    return false;
            }

        }
    }

    return true;
}


bool cech_scale(std::string input_file = "", std::string output_file = "")
{
    //steps for the algorithm
    //read and validate disk system from text file

    std::cout << "Reading from file: " << input_file << std::endl;

    std::vector< std::vector<double> > disk_system;
    if(!read_file(disk_system, input_file)){
        return false;
    }

    //calculate max vietori rips scale of the disk system
    double number_disks = disk_system.size();
    double dimentions = disk_system[0].size() - 1;
    double max_vietori_rips = std::numeric_limits<double>::lowest();
    for(unsigned i = 0; i < number_disks - 1; ++i){
        for(unsigned j = i + 1; j < number_disks; ++j){
            vietori_rips_distance = std::max(vietori_rips_distance,
                                             vietori_rips(disk_system[i][0],
                                                          disk_system[i][1],
                                                          disk_system[i][2],
                                                          disk_system[j][0],
                                                          disk_system[j][1],
                                                          disk_system[j][2]));
        }
    }


    //verify if rho is positive (if the minimum output is positive)
    if(rho_nonnegative(disk_system, vietori_rips_distance)){
        write_file(vietori_rips_distance, vietori_rips_distance, 0, 0);
    }else{

    }

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

