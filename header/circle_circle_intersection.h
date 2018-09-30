#ifndef CIRCLE_CIRCLE_INTERSECTION_H_
#define CIRCLE_CIRCLE_INTERSECTION_H_

/* circle_circle_intersection() *
 * Determine the points where 2 circles in a common plane intersect.
 *
 * bool circle_circle_intersection(
 *                                // center and radius of 1st circle
 *                                double x0, double y0, double r0,
 *                                // center and radius of 2nd circle
 *                                double x1, double y1, double r1,
 *                                // 1st intersection point
 *                                double *xi, double *yi,
 *                                // 2nd intersection point
 *                                double *xi_prime, double *yi_prime)
 *
 * This is a public domain work. 3/26/2005 Tim Voght
 *
 */
#include <cmath> //std::sqrt, std::abs, std::hypot

bool circle_circle_intersection(double x0, double y0, double r0,
                                double x1, double y1, double r1,
                                double& xi, double &yi,
                                double& xi_prime, double& yi_prime);

#endif //CIRCLE_CIRCLE_INTERSECTION_H_
