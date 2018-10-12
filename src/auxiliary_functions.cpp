#include "auxiliary_functions.h"

//############################################################################
double vectorial_distance(double x0, double y0, double x1, double y1)
{
    return std::sqrt( std::pow(x0-x1, 2.0) + std::pow(y0-y1, 2.0) );
}

//############################################################################
bool approx_equal(double d1, double d2, double tol /*= 1e-12*/)
{
    return std::abs(d1 -d2) < tol;
}

//############################################################################

bool is_left(double x0, double y0,
             double x1, double y1,
             double ax, double ay)
{
    return ( (x1 - ax)*(ay - y0) - (y1 - y0)*(ax - x0) ) < 0;
}

//############################################################################

std::vector<double> left_intersecting_point(double x0, double y0, double r0,
                                            double x1, double y1, double r1)
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
        return std::vector<double>({x_intersect1, y_intersect1});
    }else{
        return std::vector<double>({x_intersect2, y_intersect2});
    }
}

std::vector<double> left_intersection_scaled(const std::vector<double>& disk1,
                                             const std::vector<double>& disk2,
                                             double scale)
{
    return left_intersecting_point(disk1[0], disk1[1], disk1[2]*scale,
                                   disk2[0], disk2[1], disk2[2]*scale);
}

std::vector<double> right_intersection_scaled(const std::vector<double>& disk1,
                                              const std::vector<double>& disk2,
                                              double scale)
{
    double x0 = disk1[0];
    double y0 = disk1[1];
    double r0 = disk1[2]*scale;
    double x1 = disk2[0];
    double y1 = disk2[1];
    double r1 = disk2[2]*scale;
    double x_intersect1;
    double y_intersect1;
    double x_intersect2;
    double y_intersect2;

    circle_circle_intersection(x0, y0, r0,
                               x1, y1, r1,
                               x_intersect1, y_intersect1,
                               x_intersect2, y_intersect2);

    if(is_left(x0, y0, x1, y1, x_intersect1, y_intersect1)){ //reverse of left intersecting point
        return std::vector<double>({x_intersect2, y_intersect2});
    }else{
        return std::vector<double>({x_intersect1, y_intersect1});
    }
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

            bool rho_max = true;
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
    std::vector<double> intersection = left_intersection_scaled(disk1, disk2, lambda_val);

    double intersect_distance_to_disk3 = vectorial_distance(intersection[0],
                                                            intersection[1],
                                                            disk3[0],
                                                            disk3[1]);

    return disk3[2]*lambda_val - intersect_distance_to_disk3;
}

//#############################################################################

double bisection(const std::vector< std::vector<double> >& disk_system,
                 double a, double b,
                 int dig_prec)
{
    double mp = 0;
    double rho_b = rho(disk_system, b);
    double TOLERANCE = 1e-12;

    int num_iter = ceil((log10(a*(sqrt(4.0/3.0)-1)) + dig_prec)/log10(2.0));
    for(int n = 1; n <= num_iter; n++){
        mp = 0.5*(a + b);
        double rho_mp = rho(disk_system, mp);

        if(std::abs(rho_mp) <= TOLERANCE){
            break;
        }else if(rho_mp * rho_b < 0){
            a = mp;
        }else{
            b = mp;
            rho_b = rho_mp;
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

double max_vietori_rips(const std::vector< std::vector<double> >& disk_system)
{
    double vietori_rips_scale = std::numeric_limits<double>::lowest();
    for(std::vector<double> disk1 : disk_system)
        for(std::vector<double> disk2 : disk_system)
            vietori_rips_scale = std::max(vietori_rips_scale,
                                          vietori_rips(disk1, disk2));

    return vietori_rips_scale;
}

std::tuple<double, std::vector<double>>
max_vietori_rips_intersection(const std::vector< std::vector<double> >& disk_system)
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

    std::vector<double> intersection = left_intersection_scaled(disk_system[d1_idx],
                                                                disk_system[d2_idx],
                                                                vietori_rips_scale);

    return std::make_tuple(vietori_rips_scale, intersection);
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
        std::cerr << "Error. Couldn't find file with path: " << filename << std::endl;
        return false;
    }

    int number_disks;
    int dimensions;
    file >> number_disks;
    file >> dimensions;

    if(number_disks <= 0 || dimensions <= 0 ||
       (dimensions > 2 && number_disks > 3)){
        file.close();
        return false;
    }

    disk_system = std::vector< std::vector<double> >(number_disks,
                    std::vector<double>(dimensions+1));

    for(int i = 0; i < number_disks; ++i)
        for(int j = 0; j < dimensions+1; ++j)
            file >> disk_system[i][j];

    file.close(); //close the input file

    return true;
}

//#############################################################################

bool write_file(double cech_scale, double vietori_rips, std::vector<double> intersection,
                std::string filename /*= "textfiles/Cech-Scale.txt"*/)
{
    std::ofstream file(filename);
    if(!file.is_open()){
        std::cerr << "Error. Couldn't create file: " << filename << std::endl;
        return false;
    }

    if(approx_equal(cech_scale, vietori_rips)){
        file << "The Cech scale agrees with the Vietori-Rips scale." << std::endl;
    }else{
        file << "The Cech scale is greater than the Vietori-Rips scale." << std::endl;
        file << "Vietori-Rips scale: " << vietori_rips << std::endl;
    }

    file << "Cech scale: " << cech_scale << std::endl;
    file << std::endl;
    file << "The intersection point:" << std::endl;

    file << "(" << intersection[0];
    for(unsigned i = 1; i < intersection.size(); ++i)
        file << ", " << intersection[i];
    file << ")" << std::endl;

    file.close(); //close the input file

    return true;
}

//#############################################################################

std::vector< std::vector<double> > transform_disk_system(std::vector< std::vector<double> > disk_system)
{
    //create a new 3x3 disk_system in 2d
    std::vector< std::vector<double> > new_disk_system(3, std::vector<double>(3));
    double sum_diff_sq01 = 0;
    double sum_diff_sq02 = 0;
    double sum_diff_sq12 = 0;
    for(unsigned i = 0; i < disk_system[0].size() - 1; ++i){
        sum_diff_sq01 += std::pow(disk_system[0][i] - disk_system[1][i], 2);
        sum_diff_sq02 += std::pow(disk_system[0][i] - disk_system[2][i], 2);
        sum_diff_sq12 += std::pow(disk_system[1][i] - disk_system[2][i], 2);
    }

    //new disk 1
    new_disk_system[0][0] = 0;
    new_disk_system[0][1] = 0;
    new_disk_system[0][2] = disk_system[0].back();

    //new disk 2
    new_disk_system[1][0] = std::sqrt(sum_diff_sq01);
    new_disk_system[1][1] = 0;
    new_disk_system[1][2] = disk_system[1].back();

    //new disk 3
    new_disk_system[2][0] = (sum_diff_sq02 + sum_diff_sq01 - sum_diff_sq12) / (2 * std::sqrt(sum_diff_sq01));
    new_disk_system[2][1] = std::sqrt(sum_diff_sq02 - std::pow(new_disk_system[2][0], 2));
    new_disk_system[2][2] = disk_system[2].back();

    return new_disk_system;
}

//#############################################################################

std::vector<double> transform_intersection(std::vector< std::vector<double> > disk_system, std::vector< std::vector<double> > read_system, std::vector<double> c_star)
{
    int dimensions = read_system[0].size() - 1;
    double bc_numerator = disk_system[2][0]*disk_system[2][1] -
                          (disk_system[2][0] - disk_system[1][0])*disk_system[2][1];
    double bc_1 = (-disk_system[2][1]*(c_star[0] - disk_system[2][0]) +
                   (disk_system[2][0]-disk_system[1][0])*(c_star[1]-disk_system[2][1]))/bc_numerator;
    double bc_2 = (disk_system[2][1]*(c_star[0] - disk_system[2][0]) -
                   disk_system[2][0]*(c_star[1]-disk_system[2][1])) / bc_numerator;
    double bc_3 = 1 - bc_1 - bc_2;

    std::vector<double> intersection(dimensions);
    for(unsigned i = 0; i < read_system[0].size() - 1; ++i)
        intersection[i] = bc_1*read_system[0][i] + bc_2*read_system[1][i] + bc_3*read_system[2][i];

    return intersection;
}

//#############################################################################

