#ifndef AUXILIARY_FUNCTIONS_H_
#define AUXILIARY_FUNCTIONS_H_

#include <tuple> //for std::tuple
#include <cmath> //for std::sqrt, std::pow
#include <iostream> //for std::cout, std::endl
#include <fstream> // for std::ofstream, std::ifstream
#include <limits> //for std::numeric_limits<double>::lowest
#include <vector> //for std::vector
#include <string> //for std::string
#include <algorithm> //for std::min, std::max

#include "point.h"
#include "circle_circle_intersection.h"

/**
 * Calculates the euclidean distance between the points (x0, y0) (x1, y1)
 *
 * @param x0 x Coordinate of point1.
 * @param y0 y Coordinate of point1.
 * @param x1 x Coordinate of point2.
 * @param y1 y Coordinate of point2.
 *
 * @return Euclidean distance between both points.
 */
double vectorial_distance(double x0, double y0, double x1, double y1);

/**
 * Tells wheter the point (ax, ay) is to the left of the vector (x1-x0, y1-y0) or not.
 *
 * @param x0 x coordinate of starting point of vector.
 * @param y0 y coordinate of starting point of vector.
 * @param x1 x coordinate of ending point of vector.
 * @param y1 y coordinate of ending point of vector.
 * @param ax x coordinate of point to test against.
 * @param ay y coordinate of point to test against.
 *
 * @return True if the point (ax, ay) is to the left, false otherwise.
 */
bool is_left(double x0, double y0,
             double x1, double y1,
             double ax, double ay);

/**
 * Calculates the intersecting points of the 2d disks
 * d0:(x0, y0), r0 and d1:(x1, y1), r1
 * And selects the one to the left of the vector going from the center of d0 to
 * the center of d1.
 *
 * @param x0 x coordinate of disk 1.
 * @param y0 y coordinate of disk 1.
 * @param x1 x coordinate of disk 2.
 * @param y1 y coordinate of disk 2.
 *
 * @return The left intersection point of the two disks.
 */
std::vector<double> left_intersecting_point(double x0, double y0, double r0,
                              double x1, double y1, double r1);

/**
 * Calculates the intersecting points of the 2d disks
 * disk2:(x1, y1, r1) and disk2:(x2, y2, r2)
 * And selects the one to the left of the vector going from the center of disk1
 * to the center of disk2
 *
 * @param disk1 Vector with form (x1, y1, r1). x1,y1 are the coordinates of the
 *              center and r1 is its radius.
 * @param disk2 Vector with form (x2, y2, r2). x2,y2 are the coordinates of the
 *              center and r2 is its radius.
 *
 * @return The left intersection point of the two disks.
 */
std::vector<double> left_intersection_scaled(const std::vector<double>& disk1,
                               const std::vector<double>& disk2,
                               double scale);


/**
 * Calculates the vietori rips distance of two disks.
 *
 * @param x0 x coordinate of disk 1.
 * @param y0 y coordinate of disk 1.
 * @param x1 x coordinate of disk 2.
 * @param y1 y coordinate of disk 2.
 * @param r0 Radius of disk 1.
 * @param r1 Radius of disk 2.
 *
 * @return The vietori rips distance of the two disks
 */
double vietori_rips(double x0, double y0, double r0,
                    double x1, double y1, double r1);

/**
 * Calculates the vietori rips distance of two disks.
 *
 * @param disk1 Vector with form (x1, y1, r1). x1,y1 are the coordinates of the
 *              center and r1 is its radius.
 * @param disk2 Vector with form (x2, y2, r2). x2,y2 are the coordinates of the
 *              center and r2 is its radius.
 *
 * @return The vietori rips distance of the two disks
 */
double vietori_rips(const std::vector<double>& disk1,
                    const std::vector<double>& disk2);

/**
 * Calculates the maximum vietori rips distance of a system of disks.
 *
 * @param disk_system System of disks.
 *
 * @return The maximum vietori rips distance between two disks in the system.
 */
double max_vietori_rips(const std::vector< std::vector<double> >& disk_system);

/**
 * Calculates the maximum vietori rips distance of a system of disks.
 *
 * @param disk_system System of disks.
 *
 * @return The maximum vietori rips distance between two disks in the system
 *          and the intersection of the disks with this distance.
 */
std::tuple<double, std::vector<double>> max_vietori_rips_intersection(const std::vector< std::vector<double> >& disk_system);

/**
 * rho function. Evaluates the map rho_m(lambda)
 *
 * @param disk_system The system of disk in the space. The last column is the
 *                    radius of the disk and all previous are the coordinates
 *                    in the space.
 * @param lambda_val Scale to adjust the disk system.
 * @return rho
 */
double rho(const std::vector< std::vector<double> >& disk_system,
           double lambda_val);

/*
 * Es un escalar. Se calcula para una terna de discos. Por ejemplo, si tienes
 * los tres discos D1, D2 yD3, y tomas lambda=1.5, entonces Lambda_{1,2}^{3}
 * se calcula reescalando los radios de D1 y D2 por un factor de 1.5, luego
 * calculas el punto de intersección que está a la izquierda del vector que
 * va del centro de D1 al centro de D2 (esto es d_ij(lambda) = d_{1,2}(1.5)),
 * después calculas la distancia al centro de D3, esto es: ||d_12(1.5) - c3||.
 * Esa cantidad se la restas al radio de D3 reescalado por lambda=1.5.
 *
 * Given three disks in a space, lambda is the distance from the third disk
 * scaled by a number to the left intersection of the disk 1 and 2 scaled by
 * the same number.
 *
 * @param disk1 Vector with the positions of disk1 and its radius as last entry.
 * @param disk2 Vector with the positions of disk2 and its radius as last entry.
 * @param disk3 Vector with the positions of disk3 and its radius as last entry.
 * @param lambda_val The radius are scaled by this number.
 *
 * @return The distance from the third disk scaled by a number to the left
 * intersection of the disk 1 and 2 scaled by the same number.
 */

double lambda(const std::vector<double>& disk1,
              const std::vector<double>& disk2,
              const std::vector<double>& disk3,
              double lambda_val);

//#############################################################################

/**
 * Applies the bisection numeric method to the rho function (defined above)
 * using the specified disk system.
 *
 * @param disk_system Matrix which rows are disks in space and columns the
 *                    coordinates of those disks. The last column is the
 *                    radius of the disks.
 * @param a One of two starting points for bisection.
 * @param b One of two starting points for bisection.
 * @param dig_prec Digit precision. Used to calculate the number of iterations.
 *
 * @return Value such as rho(disk_system, value) is as close to zero as the
 *         method could give.
 */
double bisection(const std::vector< std::vector<double> >& disk_system,
                 double a, double b, int dig_prec);

//#############################################################################

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

//#############################################################################

/**
 * Reads the file and stores the contents in the disk system.
 *
 * @param disk_system Matrix of coordinates of the disks. Each row has the
 *                    coordinates of the disks. And the last column has the
 *                    radius of the disk.
 * @param filename The file's name where the disk system is defined.
 *                 The first number must be the number of disks, the second the
 *                 number of dimentions and then rows of coordinates of each
 *                 disk and the last column must be the radius of the disk.
 *                 "Disk-system.txt" by default
 *
 * @return True if the file could be opened and everything went all right, false otherwise
 */
bool read_file(std::vector< std::vector<double> >& disk_system,
                        std::string filename = "textfiles/Disks-system.txt");

//#############################################################################

/**
 * Reads the file and stores the contents in the disk system.
 * Writes to a file the results of the program.
 * Writes whether the cech scale is equal to the vietori rips scale and writes
 * the cech scale of the disk system.
 * Also outputs the intersection point.
 *
 * @param cech_scale The cech scale calculated by the program.
 * @param vietori_rips The vietori rips calculated by the program.
 * @param intersection Vector with the coordinates of the intersecion
 *                     calculated by the program.
 * @param filename The file's name to output the result to.
 *                 By default is "textfiles/Cech-Scale.txt"
 *
 * @return True if the file could be opened and everything went all right, false otherwise
 */
bool write_file(double cech_scale, double vietori_rips, std::vector<double> intersection,
                std::string filename = "textfiles/Cech-Scale.txt");

//#############################################################################

/**
 * Given a disk system of 3 disks and a dimention greater than 2, this function
 * applies a projection to the disk system to transform it into a system of 3
 * disks and dimention 2.
 *
 * @param disk_system The disk system of 3 disks and dimention greater than 2
 *
 * @return The projection of the disk system into a 2d plane: a disk system of
 *         3 disks and dimention 2.
 */
std::vector< std::vector<double> > transform_disk_system(std::vector< std::vector<double> > disk_system);

//#############################################################################

/**
 * Given a disk system of 3 disks, dimention greater than 2, its projection,
 * to a 2d plane, and the intersection of the projected system, calculates the
 * corresponding intersection of the original disk system.
 *
 * @param disk_system The projection of the disk system.
 * @param read_system The original disk system.
 * @param c_star The intersection of the projected disk system.
 *
 * @return The corresponding intersection of the original disk system.
 */
std::vector<double> transform_intersection(std::vector< std::vector<double> > disk_system, std::vector< std::vector<double> > read_system, std::vector<double> c_star);

//#############################################################################

#endif //AUXILIARY_FUNCTIONS_H_

