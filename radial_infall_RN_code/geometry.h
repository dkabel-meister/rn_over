// geometry.h
// Header file describing the functions performing various geometrical
// operations (intersections, area calculation)
//
// Author: Roland Haas and Peter Zimmerman
// File created Jun 18, 2006
//


#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

class Geometry
{
    public:
        double area_triangle(
                const double x1, const double t1, 
                const double x2, const double t2, 
                const double x3, const double t3
                );

        double area_quadrangle(
                const double x1, const double t1, 
                const double x2, const double t2, 
                const double x3, const double t3, 
                const double x4, const double t4
                );

        double area_pentangle(
                const double x1, const double t1, 
                const double x2, const double t2, 
                const double x3, const double t3,
                const double x4, const double t4, 
                const double x5, const double t5 
                ); 

        double area_hexagon(
                const double x1, const double t1,
                const double x2, const double t2,
                const double x3, const double t3, 
                const double x4, const double t4,
                const double x5, const double t5, 
                const double x6, const double t6
                );

        static void intersect_line(
                const double u1, const double v1,
                const double u2, const double v2,
                const double X1, const double Y1,
                const double X2, const double Y2,
                double *a, double *b
                );
};
#endif //_GEOMETRY_H_
