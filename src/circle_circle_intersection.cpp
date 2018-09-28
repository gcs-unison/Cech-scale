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
  //d = sqrt((dy*dy) + (dx*dx));
  d = hypot(dx,dy); // Suggested by Keith Briggs

  /* Check for solvability. */
  if (d > (r0 + r1) && fabs(d - (r0+r1)) > 1e-12)
  {
    /* no solution. circles do not intersect. */
    return false;
  }
  if (d < fabs(r0 - r1) && fabs(d - fabs(r0-r1)) > 1e-12)
  {
    /* no solution. one circle is contained in the other */
    return false;
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
  h = std::sqrt( std::min(0.0, std::abs(r0*r0 - a*a)) );

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

