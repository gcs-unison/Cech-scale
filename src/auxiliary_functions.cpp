#include "auxiliary_functions.h"

//############################################################################
double vectorial_distance(double x0, double y0, double x1, double y1)
{
    return sqrt( pow(x0-x1, 2.0) + pow(y0-y1, 2.0) );
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

//############################################################################

double rho_improved(const std::vector< std::vector<double> >& disk_system,
                    double lambda_val,
                    double tolerance = 1e-12)
{
    unsigned number_disks = disk_system.size();
    unsigned dimentions = disk_system[0].size() - 1;
    double rho_value = std::numeric_limits<double>::lowest();
    //substract 2 from the size for the iterations where k==i and k==j
    std::vector<double> lambdas(number_disks - 2);
    lambdas.clear(); //set size to 0 but mantain capacity

    for(unsigned i = 0; i < number_disks; ++i){
        for(unsigned j = 0; j < dimentions; ++j){
            if(i == j)
                continue;

            bool update_rho = true;

            for(unsigned k = 0; k < number_disks; ++k){
                if(k == i || k == j)
                    continue;

                lambdas.push_back(lambda(disk_system[i],
                                         disk_system[j],
                                         disk_system[k],
                                         lambda_val,
                                         tolerance));

                if(lambdas.back() <= rho_value){
                    update_rho = false;
                    break;
                }
            }

            if(update_rho){
                //store the smaller lambda calculated in rho_value
                rho_value = *std::min_element(std::begin(lambdas), std::end(lambdas));
            }

            lambdas.clear();
        }
    }

    return rho_value;
}

//############################################################################

double lambda(const std::vector<double>& disk1,
              const std::vector<double>& disk2,
              const std::vector<double>& disk3,
              double lambda_val,
              double tolerance = 1e-12)
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

//############################################################################

double rho(double M[3][3], int m, double lambda, double prec){
    // M is 2-disks system.
    // m is the number is disks in M. In our script we will fix m=3.
    // lambda is the scale at which the function rho_M is evaluated.
    int j, D, k, k_index;
    int rho_ij_is_higher, rho_ji_is_higher, double_intersection;
    int k_values[m-2];

    double rho_value = std::numeric_limits<double>::lowest();
    double i_radious_reescaled, j_radious_reescaled, dist_ij, r_ij, a, b;
    double sum, distdkcenter; //auxiliar variables
    double i_center[2], j_center[2];
    double c_ij[2];
    double d_ij[2], d_ji[2];
    double k_center[2];
    double n_ij[2];
    double rho_ijk[m-2], rho_jik[m-2];

    for(int i=0;i<m-1;i++){
        for(j=i+1;j<m;j++){
            // Calculate the intersection points d_ij:
            for(D=0; D < 2; D++)
                i_center[D]=M[i][D];

            //should this line be inside for??
            i_radious_reescaled = lambda*M[i][2];
            for(D=0; D < 2; D++)
                j_center[D]=M[j][D];

            //should this line be inside for??
            j_radious_reescaled = lambda*M[j][2];

            for(D=0; D < 2; D++)
                c_ij[D]=j_center[D]-i_center[D];
            sum=0;
            for(D=0; D < 2; D++){
                sum+=pow(i_center[D]-j_center[D],2);
            }
            dist_ij=sqrt(sum);
            double_intersection=0;

            // Radii close enough, are equal_
            if(fabs(i_radious_reescaled-j_radious_reescaled) <= prec){
                j_radious_reescaled = i_radious_reescaled;
            }

            // Concentric disks
            if(dist_ij <= prec){
                for(D=0; D < 2; D++)
                    d_ij[D] = i_center[D];
            }else{
                // External tangency
                if(fabs(dist_ij - (i_radious_reescaled+j_radious_reescaled)) < prec){
                    for(D=0; D < 2; D++){
                        double dividend = (j_radious_reescaled*i_center[D] + i_radious_reescaled*j_center[D]);
                        double divisor = (i_radious_reescaled+j_radious_reescaled);
                        d_ij[D] = dividend / divisor;
                    }
                }else{
                    // Complete contention
                    if((dist_ij > prec) && ( dist_ij <= fabs(i_radious_reescaled-j_radious_reescaled)+prec)){
                        r_ij = fabs(i_radious_reescaled-j_radious_reescaled);
                        if(j_radious_reescaled > i_radious_reescaled){
                            for(D=0; D < 2; D++)
                                i_center[D] = j_center[D];

                            i_radious_reescaled = j_radious_reescaled;

                            for(D=0; D < 2; D++)
                                c_ij[D] = -c_ij[D];
                        }
                        for(D=0; D < 2; D++)
                            d_ij[D] = i_center[D] + i_radious_reescaled*c_ij[D]/r_ij;
                    }else{
                        // Two intersection points
                        n_ij[0]=-c_ij[1]; n_ij[1]=c_ij[0];

                        // Calculate the coefficients a and b such that
                        // d_ij = i_center + a * c_ij + b * n_ij and d_ji = i_center + a * c_ij - b * n_ij
                        a = 0.5*(1 + (i_radious_reescaled-j_radious_reescaled)*(i_radious_reescaled+j_radious_reescaled)/pow(dist_ij,2));
                        b = sqrt(pow(i_radious_reescaled,2) - pow(a*dist_ij,2))/dist_ij;
                        for(D=0; D < 2; D++)
                            d_ij[D] = i_center[D] + a * c_ij[D] + b * n_ij[D];
                        for(D=0; D < 2; D++)
                            d_ji[D] = i_center[D] + a * c_ij[D] - b * n_ij[D];
                        double_intersection = 1;
                    }
                }
            }

            // Take k such that k != i,j.
            int count = 0;
            for(D=0; D<m; D++)
                if(D!=i && D!=j)
                    k_values[count++] = D;

            rho_ij_is_higher = 1;
            k_index = 0;
            while(rho_ij_is_higher==1 && k_index != m-2){
                k = k_values[k_index];
                for(D=0; D < 2; D++)
                    k_center[D]=M[k][D];

                distdkcenter=0;
                for(D=0; D < 2; D++){
                    distdkcenter+=pow(d_ij[D]-k_center[D],2);
                }
                distdkcenter=sqrt(distdkcenter);
                rho_ijk[k_index] = lambda*M[k][2]-distdkcenter;

                if(rho_ijk[k_index] < rho_value){
                    rho_ij_is_higher = 0;
                }
                k_index++;
            }

            if(rho_ij_is_higher==1){
                rho_value=rho_ijk[0];
                for(D=1;D<m-2;D++)
                    if(rho_ijk[D]<rho_value)
                        rho_value=rho_ijk[D];
            }

            if(double_intersection==1){
                k_index = 0;
                rho_ji_is_higher = 1;

                while(rho_ji_is_higher==1 && k_index != m-2){
                    k = k_values[k_index];
                    for(D=0; D < 2; D++)
                        k_center[D]=M[k][D];

                    distdkcenter=0;
                    for(D=0; D < 2; D++){
                        distdkcenter+=pow(d_ji[D]-k_center[D],2);
                    }
                    distdkcenter=sqrt(distdkcenter);
                    rho_jik[k_index] = lambda*M[k][2]-distdkcenter;

                    if(rho_jik[k_index] < rho_value){
                        rho_ji_is_higher = 0;
                    }
                    k_index++;
                }

                if(rho_ji_is_higher==1){
                    rho_value=rho_jik[0];
                    for(D=1;D<m-2;D++)
                        if(rho_jik[D]<rho_value)
                            rho_value=rho_jik[D];
                }
            }
        }
    }
    return rho_value;
}

//#############################################################################

double bisection(double M[3][3], int m, double a, double b, int dig_prec, double prec){
    // dig_prec = number of digits of precision.
    double mp;
    double rho_M_mp;
    double rho_M_a;
    double rho_M_b;
    int n;
    int num_iter;

    rho_M_a = rho(M, m, a, prec);
    rho_M_b = rho(M, m, b, prec);
    num_iter = ceil((log10(a*(sqrt(4.0/3)-1)) + dig_prec)/log10(2));
    for(n=1; n <= num_iter; n++){
        mp = 0.5*(a + b);
        rho_M_mp = rho(M, m, mp, prec);

        if(fabs(rho_M_mp) <= prec){
            break;
        }else{
            if(rho_M_mp * rho_M_b < 0){
                a = mp;
            }else{
                b = mp;
                rho_M_b = rho_M_mp;
            }
        }
    }
    return mp;
}

//#############################################################################

double bisection_improved(const std::vector< std::vector<double> >& disk_system,
                          double a, double b,
                          int dig_prec, double prec){
    double mp = 0;
    double rho_b = rho_improved(disk_system, b, prec);

    int num_iter = ceil((log10(a*(sqrt(4.0/3.0)-1)) + dig_prec)/log10(2.0));
    for(int n=1; n <= num_iter; n++){
        mp = 0.5*(a + b);
        double rho_mp = rho_improved(disk_system, mp, prec);

        if(fabs(rho_mp) <= prec){
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

bool read_file(int& m, int& d, double M[3][3],
              std::vector<std::vector<double>>& N,
              std::string filename = "textfiles/Disk-system.txt")
{

    std::ifstream file(filename);
    if(!file.is_open()){
        std::cout << "Error. Couldn't find file with path: " << filename << std::endl;
        return false;
    }

    file >> m;
    file >> d;

    //initialize N's size knowing d
    N = std::vector< std::vector<double> >(3, std::vector<double>(d+1));

    if(d==2)
        for(int i=0;i<m;i++)
            for(int j=0;j<3;j++)
                file >> M[i][j];
    else if(m==3)
        for(int i=0; i<3; i++)
            for(int j=0; j<d+1; j++)
                file >> N[i][j];

    file.close(); //close the input file

    return true;
}

//#############################################################################

bool read_file_improved(std::vector< std::vector<double> >& disk_system,
                        std::string filename = "textfiles/Disk-system.txt")
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

    //initialize N's size knowing d
    disk_system = std::vector< std::vector<double> >(number_disks,
                    std::vector<double>(dimentions+1));

    for(int i = 0; i < number_disks; ++i)
        for(int j = 0; j < dimentions+1; ++j)
            file >> disk_system[i][j];

    file.close(); //close the input file

    return true;
}

