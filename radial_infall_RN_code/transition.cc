#include <cmath>
#include "field.h"

#if 0
double Field::transitionSource(const double t)
{
    if (t<t_switch)
    {
        const double tan_term = tan(M_PI*t/2./t_switch);
        return 0.5 + 0.5*tanh(1/M_PI*(tan_term-1./tan_term));
    }
    else 
        return 1.0;
}
#endif

double Field::transitionSource(const double t)
{
    //double k = 0.1;

    if (t<t_switch)
        return t/t_switch;
    else 
        return 1.0;
}
#if 0
double Field::transitionSource(const double t)
{
    double k = 10;

    //double t_switch = 100;

    if (t<t_switch)
        return tanh(k*t/t_switch);
    else 
        return 1.0;
}
#endif
