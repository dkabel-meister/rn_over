// geometry.cc
// contains functions to calculate areas of triangles and pentangles
// formulae are taken from Bronstein (3.278) on page 189
// contains a function to intersect two straight lines specified 
// by four points
//
// Author: Roland Haas and Peter Zimmerman
// File created on Feb 23, 2006

#include <cmath>
#include "geometry.h"

//
// area_triangle
// calculates the area of a triangle given its vertices
// 
double Geometry::area_triangle(const double x1,const double t1, const double x2, const double  t2, const double x3, const double  t3)
{
	return 0.5*fabs((x1-x2)*(t1+t2)+(x2-x3)*(t2+t3)+(x3-x1)*(t3+t1));
}

//
// area_quadrangle
// calculates the area of a quadrangle given its vertices
// 
double Geometry::area_quadrangle(const double x1,const double t1, const double x2, const double  t2, const double x3, const double  t3, const double x4, const double  t4)
{
	return 0.5*fabs((x1-x2)*(t1+t2)+(x2-x3)*(t2+t3)+(x3-x4)*(t3+t4)+(x4-x1)*(t4+t1));
}

//
// area_pentangle
// calculates the area of a pentangle given its vertices
// 
double Geometry::area_pentangle(const double x1,const double t1, const double x2, const double  t2, const double x3, const double  t3, const double x4, const double  t4, const double x5, const double  t5)
{
	return 0.5*fabs((x1-x2)*(t1+t2)+(x2-x3)*(t2+t3)+(x3-x4)*(t3+t4)+(x4-x5)*(t4+t5)+(x5-x1)*(t5+t1));
}

//
// area_hexagon
// calculates the area of a hexagon given its vertices
// 
double Geometry::area_hexagon(const double x1,const double t1, const double x2, const double  t2, const double x3, const double  t3, const double x4, const double  t4, const double x5, const double  t5, const double x6, const double  t6)
{
	return 0.5*fabs((x1-x2)*(t1+t2)+(x2-x3)*(t2+t3)+(x3-x4)*(t3+t4)+(x4-x5)*(t4+t5)+(x5-x6)*(t5+t6)+(x6-x1)*(t6+t1));
}

void Geometry::intersect_line(const double u1, const double v1, const double u2, const double v2, 
        const double X1, const double Y1, const double X2, const double Y2, double *a, double *b)
{
    // code produced by MAPLE
    // [> restart;
    // [> sols:=solve({u1+lambda*(u2-u1)=X1+mu*(X2-X1), v1+lambda*(v2-v1)= Y1+mu*(Y2-Y1)}, {mu, lambda});
    // [> ab:=subs(op(sols), [u1+lambda*(u2-u1), v1+lambda*(v2-v1)]);
    // [> CodeGeneration['C'](ab[1]);
    // [> CodeGeneration['C'](ab[2]);

    *a = u1 + (-X1 * v1 - X2 * Y1 + X2 * v1 + Y2 * X1 - Y2 * u1 + Y1 * u1) / (X1 * v2 - X1 * v1 - X2 * v2 + X2 * v1 + Y2 * u2 - Y2 * u1 - Y1 * u2 + Y1 * u1) * (u2 - u1);
    *b = v1 + (-X1 * v1 - X2 * Y1 + X2 * v1 + Y2 * X1 - Y2 * u1 + Y1 * u1) / (X1 * v2 - X1 * v1 - X2 * v2 + X2 * v1 + Y2 * u2 - Y2 * u1 - Y1 * u2 + Y1 * u1) * (v2 - v1);
}
