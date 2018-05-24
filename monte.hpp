/* type definitions for monte carlo simulation */

/*
 * HISTORY
 * 16-Oct-88 Jeffrey Trull (jt1j) at Carnegie-Mellon University
 *       Created
 */

/* note: in actual fact I am typing this in, 30 years later, from a paper copy,
 * and being a bit inexact about it
 */

#ifndef MONTE_HPP
#define MONTE_HPP

// constants

/*
 * evaluate the constant of proportionality between the acoustic
 * scattering rate and the square root of the electron energy
 */
double const scatter_const = (0.449E18)*(0.015813)*(300.0)*(49)/((5.37)*(270.4E9));

/*
 * evaluate the accelerations due to the constant electric field
 * between collisions
 */
double const accel_const = 2.0*3.14159*(1.6E-19)*(10000.0)/(6.63E-34);

/*
 * evaluate k-to-velocity conversion constant
 */
double const vel_const = (10000.0)*(6.63E-34)/((2.0*3.14159)*(0.063)*(9.11E-30));

// types

template<typename Float>
struct vector_str {
    vector_str() : x(0.0), y(0.0), z(0.0) {}
    Float x, y, z;
};

template<typename Float>
using vector_ptr = vector_str<Float> *;

template<typename Float>
using vector_cptr = vector_str<Float> const *;

// utility functions

/* k to energy conversion */
template<typename Float>
Float energy_from_k(vector_cptr<Float> kvecptr)
/*
 * determine squared velocity, then multiply by effective mass
 * and divide by two, converting to electron volts
 */
{
    Float    vx, vy, vz;      /* velocity components */
    Float    temp;

    vx = vel_const * kvecptr->x;
    vy = vel_const * kvecptr->y;
    vz = vel_const * kvecptr->z;

    temp = vx*vx + vy*vy + vz*vz;

    temp = (0.5)*(0.063)*(9.11E-30)*temp*(0.0001)/(1.6E-19);

    return temp;
}

template<typename Float>
void    coord_convert(vector_cptr<Float> firstkptr,
                      Float theta_r,
                      Float phi_r,
                      vector_ptr<Float> lastkptr)
{
    Float               theta, phi;
    vector_str<Float>   k2prime;
    Float               totalk;

    phi   = atan(firstkptr->y/firstkptr->x);
    theta = atan(sqrt(firstkptr->x*firstkptr->x +
                      firstkptr->y*firstkptr->y)/firstkptr->z);

    totalk = sqrt(firstkptr->x*firstkptr->x +
                  firstkptr->y*firstkptr->y +
                  firstkptr->z*firstkptr->z);

    k2prime.x = totalk * sin(theta_r) * cos(phi_r);
    k2prime.y = totalk * sin(theta_r) * sin(phi_r);
    k2prime.z = totalk * cos(theta_r);

    lastkptr->x = cos(phi) * cos(theta) * k2prime.x
        + cos(phi) * sin(theta) * k2prime.z
        - sin(phi) * k2prime.y;

    lastkptr->y = sin(phi) * cos(theta) * k2prime.x
        + sin(phi) * sin(theta) * k2prime.z
        + cos(phi) * k2prime.y;

    lastkptr->z = cos(theta) * k2prime.z - sin(theta) * k2prime.z;
}

#endif // MONTE_HPP
