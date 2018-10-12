#include "circle_circle_intersection.h"

bool circle_circle_intersection(double x0, double y0, double r0,
                                double x1, double y1, double r1,
                                double& xi, double &yi,
                                double& xi_prime, double& yi_prime)
{
    double a, dx, dy, d, h, rx, ry;
    double x2, y2;

    /* dx and dy are the vertical and horizontal distances between
    * the circle centers.
    */
    dx = x1 - x0;
    dy = y1 - y0;

    /* Determine the straight-line distance between the centers. */
    //d = std::sqrt((dy*dy) + (dx*dx));
    d = std::hypot(dx,dy); // Suggested by Keith Briggs

    /* Check for solvability. */
    if (d > (r0 + r1) && approx_equal(d, r0+r1))
    /* no solution. circles do not intersect. */
        return false;

    if (d < std::abs(r0 - r1) && approx_equal(d, std::abs(r0-r1))){
        //no solution.
        //One circle is contained in the other so intersection is
        //defined as follows
        double dr = std::abs(r0 - r1);
        if(r0 > r1){
            xi = x0 + r0*(x1-x0)/dr;
            yi = y0 + r0*(y1-y0)/dr;
        }else{
            xi = x1 + r1*(x0-x1)/dr;
            yi = y1 + r1*(y0-y1)/dr;
        }
        xi_prime = xi;
        yi_prime = yi;

        return true;
    }

    /* 'point 2' is the point where the line through the circle
    * intersection points crosses the line between the circle
    * centers.
    */

    /* Determine the distance from point 0 to point 2. */
    a = ((r0*r0) - (r1*r1) + (d*d)) / (2.0 * d) ;

    /* Determine the coordinates of point 2. */
    x2 = x0 + (dx * a/d);
    y2 = y0 + (dy * a/d);

    /* Determine the distance from point 2 to either of the
    * intersection points.
    */
    h = std::sqrt( std::abs(r0*r0 - a*a) );

    /* Now determine the offsets of the intersection points from
    * point 2.
    */
    rx = -dy * (h/d);
    ry = dx * (h/d);

    /* Determine the absolute intersection points. */
    xi = x2 + rx;
    xi_prime = x2 - rx;
    yi = y2 + ry;
    yi_prime = y2 - ry;

    return true;
}

