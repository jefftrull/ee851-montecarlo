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

#include <boost/math/constants/constants.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>

// constants
constexpr double two_pi  = boost::math::constants::two_pi<double>();
constexpr double echarge = boost::units::si::constants::codata::e / boost::units::si::coulomb;
constexpr double dirac   = boost::units::si::constants::codata::hbar /
                           (boost::units::si::joule * boost::units::si::second);
constexpr double T       = 300.0;  // room temperature in Kelvin

/*
 * evaluate the constant of proportionality between the acoustic
 * scattering rate and the square root of the electron energy
 */
// components of the scattering rate constant:
constexpr double E1 = 7;       // Acoustic deformation potential in eV
constexpr double rho = 5.37;   // Crystal density in g/cm^3
constexpr double u = 5.2E5;    // speed of sound in crystal, cm/s
constexpr double meff = 0.063; // effective electron mass in gamma (000) valley (unitless)

constexpr double scatter_const =
    (0.449E18 * pow(meff, 1.5) * T * (E1*E1)) /
    (rho * (u*u)) ;

/*
 * evaluate the accelerations due to the constant electric field
 * between collisions
 */
constexpr double Efield = 10000.0;   // Applied field (V/cm)
constexpr double accel_const = echarge * Efield / dirac ;

/*
 * evaluate k-to-velocity conversion constant
 */
constexpr double m_e       = boost::units::si::constants::codata::m_e /   // electron mass
                             boost::units::si::kilogram ;
constexpr double vel_const = (Efield * dirac)/(meff * m_e) ;


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

    temp = 0.5 * meff * m_e * temp/(echarge * Efield);

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
