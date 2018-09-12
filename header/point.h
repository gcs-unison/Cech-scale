#ifndef POINT_H_
#define POINT_H_

//simple 2d point structure.
struct Point{
    double x;
    double y;

    Point(){}
    Point(double x, double y): x(x), y(y) {}
};


#endif //POINT_H_

