#ifndef AUXILIARY_FUNCTIONS_H_
#define AUXILIARY_FUNCTIONS_H_

#include <stdlib.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <math.h>
#include <fstream>
#include <cstring>
#include <limits> //for std::numeric_limits<double>::lowest
#include <vector> //for std::vector
#include <string> //for std::string
#include <algorithm> //for std::min_element

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
 * @param ax Reference where the x coordinate of left intersection will be stored.
 * @param ay Reference where the y coordinate of left intersection will be stored.
 *
 */
void left_intersecting_point(double x0, double y0, double r0,
                             double x1, double y1, double r1,
                             double& ax, double& ay);

/**
 * Improved rho function. Evaluates the map rho_m(lambda)
 *
 * @param disk_system The system of disk in the space. The last column is the
 *                    radius of the disk and all previous are the coordinates
 *                    in the space.
 * @param lambda Scale to adjust the disk system.
 * @return rho
 */
double rho_improved(const std::vector< std::vector<double> >& disk_system,
                    double lambda_val,
                    double tolerance/* = 1e-12*/);

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
 * @param tolerance The tolerance to declare that two disks are close enough to be
 *                  the same. Used by the intersection calculating function.
 *
 * @return The distance from the third disk scaled by a number to the left
 * intersection of the disk 1 and 2 scaled by the same number.
 */

double lambda(const std::vector<double>& disk1,
              const std::vector<double>& disk2,
              const std::vector<double>& disk3,
              double lambda_val,
              double tolerance /*= 1e-12*/);

/*      Function RHO     */
// For a definition, see the paper "A numerical approach for the filtered generalized Cech complex".

double rho(double M[3][3], int m, double lambda, double prec);

//#############################################################################

/*      Bisection's numerical method      */
double bisection(double M[3][3], int m, double a, double b, int dig_prec, double prec);

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
 * @param prec Precision used to decide whether a solution is good enough.
 *
 * @return Value such as rho(disk_system, value) is as close to zero as the
 *         method could give.
 */
double bisection_improved(const std::vector< std::vector<double> >& disk_system,
                          double a, double b,
                          int dig_prec, double prec);
//#############################################################################

/**
 * Reads the file and stores the contents in m and d first and the rest in M.
 *
 * @param m Number of disks.
 * @param d Number of dimensions of the disks.
 * @param M Matrix of coordinates of the disks. Each row has the coordinates of the disks.
 * @param N Matrix of coordinates of the disks. Each row has the coordinates of the disks.?????????
 * @param file The file's name. "Disk-system.txt" by default
 * @return True if the file could be opened and everything went all right, false otherwise
 */
bool read_file(int& m, int& d, double M[3][3],
              std::vector<std::vector<double>>& N,
              std::string filename /*= "textfiles/Disk-system.txt"*/);


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
bool read_file_improved(std::vector< std::vector<double> >& disk_system,
                        std::string filename /*= "textfiles/Disk-system.txt"*/);


#endif //AUXILIARY_FUNCTIONS_H_
