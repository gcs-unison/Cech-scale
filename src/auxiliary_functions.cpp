#include "auxiliary_functions.h"

//############################################################################
double vectorial_distance(double x0, double y0, double x1, double y1)
{
    return std::sqrt( std::pow(x0-x1, 2.0) + std::pow(y0-y1, 2.0) );
}

//############################################################################

bool is_left(double x0, double y0,
             double x1, double y1,
             double ax, double ay)
{
    return ( (x1 - ax)*(ay - y0) - (y1 - y0)*(ax - x0) ) < 0;
}

//############################################################################

void left_intersecting_point(double x0, double y0, double r0,
                             double x1, double y1, double r1,
                             double& ax, double& ay)
{
    double x_intersect1;
    double y_intersect1;
    double x_intersect2;
    double y_intersect2;

    circle_circle_intersection(x0, y0, r0,
                               x1, y1, r1,
                               x_intersect1, y_intersect1,
                               x_intersect2, y_intersect2);

    if(is_left(x0, y0, x1, y1, x_intersect1, y_intersect1)){
        ax = x_intersect1;
        ay = y_intersect1;
    }else{
        ax = x_intersect2;
        ay = y_intersect2;
    }
}

void left_intersecting_point(const std::vector<double>& disk1,
                             const std::vector<double>& disk2,
                             double& ax, double& ay)
{
    left_intersecting_point(disk1[0], disk1[1], disk1[2],
                            disk2[0], disk2[1], disk2[2],
                            ax, ay);
}

//############################################################################

double rho(const std::vector< std::vector<double> >& disk_system,
                    double lambda_val)
{
    double rho_value = std::numeric_limits<double>::lowest();

    for(std::vector<double> disk1 : disk_system){
        for(std::vector<double> disk2 : disk_system){
            if(disk1 == disk2)
                continue;

            double rho_max = true;
            double min_rho = std::numeric_limits<double>::max();
            for(std::vector<double> disk3 : disk_system){
                if(disk3 == disk1 || disk3 == disk2)
                    continue;

                double lambda_k = lambda(disk1, disk2, disk3, lambda_val);
                min_rho = std::min(min_rho, lambda_k);
                if(lambda_k <= rho_value){
                    rho_max = false;
                    break;
                }
            }

            if(rho_max)
                rho_value = min_rho;
        }
    }

    return rho_value;
}

//############################################################################

double lambda(const std::vector<double>& disk1,
              const std::vector<double>& disk2,
              const std::vector<double>& disk3,
              double lambda_val)
{
    double x_intersect;
    double y_intersect;
    double scaled_radious1 = disk1.back()*lambda_val;
    double scaled_radious2 = disk2.back()*lambda_val;
    double scaled_radious3 = disk3.back()*lambda_val;

    left_intersecting_point(disk1[0], disk1[1], scaled_radious1,
                            disk2[0], disk2[1], scaled_radious2,
                            x_intersect, y_intersect);

    double intersect_distance_to_disk3 = vectorial_distance(x_intersect, y_intersect,
                                                            disk3[0], disk3[1]);

    return scaled_radious3 - intersect_distance_to_disk3;
}

//#############################################################################

double bisection(const std::vector< std::vector<double> >& disk_system,
                 double a, double b,
                 int dig_prec){
    double mp = 0;
    double rho_b = rho(disk_system, b);
    double TOLERANCE = 1e-12;

    int num_iter = ceil((log10(a*(sqrt(4.0/3.0)-1)) + dig_prec)/log10(2.0));
    for(int n=1; n <= num_iter; n++){
        mp = 0.5*(a + b);
        double rho_mp = rho(disk_system, mp);

        if(std::abs(rho_mp) <= TOLERANCE){
            break;
        }else{
            if(rho_mp * rho_b < 0){
                a = mp;
            }else{
                b = mp;
                rho_b = rho_mp;
            }
        }
    }

    return mp;
}

//#############################################################################

double vietori_rips(double x0, double y0, double r0,
                    double x1, double y1, double r1)
{
    return vectorial_distance(x0, y0, x1, y1) / (r0 + r1);
}

double vietori_rips(const std::vector<double>& disk1,
                    const std::vector<double>& disk2)
{
    return vietori_rips(disk1[0], disk1[1], disk1[2],
                        disk2[0], disk2[1], disk2[2]);
}

//#############################################################################

VietoriRipsScale max_vietori_rips(const std::vector< std::vector<double> >& disk_system)
{
    double number_disks = disk_system.size();
    double vietori_rips_scale = std::numeric_limits<double>::lowest();
    unsigned d1_idx = 0;
    unsigned d2_idx = 1;
    for(unsigned i = 0; i < number_disks - 1; ++i){
        for(unsigned j = i + 1; j < number_disks; ++j){


            double vietori_rips_scale_ij = vietori_rips(disk_system[i],
                                                        disk_system[j]);

            if(vietori_rips_scale_ij > vietori_rips_scale){
                vietori_rips_scale = vietori_rips_scale_ij;
                d1_idx = i;
                d2_idx = j;
            }
        }
    }

    double intersection_x;
    double intersection_y;
    const std::vector<double>& disk1 = disk_system[d1_idx];
    const std::vector<double>& disk2 = disk_system[d2_idx];
    left_intersecting_point(disk1[0], disk1[1], disk1[2]*vietori_rips_scale,
                            disk2[0], disk2[1], disk2[2]*vietori_rips_scale,
                            intersection_x, intersection_y);

    return VietoriRipsScale(vietori_rips_scale, intersection_x, intersection_y);
}

//#############################################################################

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
                    break;
                }
            }
            if(nonnegative)
                return true;
        }
    }

    return false;
}

//#############################################################################

bool read_file(std::vector< std::vector<double> >& disk_system,
                        std::string filename /*= "textfiles/Disks-system.txt"*/)
{

    std::ifstream file(filename);
    if(!file.is_open()){
        std::cout << "Error. Couldn't find file with path: " << filename << std::endl;
        return false;
    }

    int number_disks;
    int dimentions;
    file >> number_disks;
    file >> dimentions;

    if(number_disks <= 0 || dimentions <= 0 ||
       (dimentions > 2 && number_disks > 3)){
        file.close();
        return false;
    }

    //initialize N's size knowing d
    disk_system = std::vector< std::vector<double> >(number_disks,
                    std::vector<double>(dimentions+1));

    for(int i = 0; i < number_disks; ++i)
        for(int j = 0; j < dimentions+1; ++j)
            file >> disk_system[i][j];

    file.close(); //close the input file

    return true;
}

//#############################################################################

bool write_file(double cech_scale, double vietori_rips, std::vector<double> intersection,
                std::string filename /*= "textfiles/Cech-Scale.txt"*/)
{
    const double TOLERANCE = 1e-12;

    std::ofstream file(filename);
    if(!file.is_open()){
        std::cout << "Error. Couldn't create file: " << filename << std::endl;
        return false;
    }

    //checks check_scale == vietori_rips
    if(std::abs(cech_scale - vietori_rips) < TOLERANCE){
        file << "The cech scale and vietori rips coincide." << std::endl;
    }else{
        file << "The cech scale and vietori rips DO NOT coincide." << std::endl;
    }

    file << "Cech scale: " << cech_scale << std::endl;
    file << "The intersection point:" << std::endl;

    file << "(" << intersection[0];
    for(unsigned i = 1; i < intersection.size(); ++i)
        file << ", " << intersection[i];
    file << ")" << std::endl;

    file.close(); //close the input file

    return true;
}

//#############################################################################

