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

#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <boost/units/systems/si.hpp>
#include <boost/units/systems/si/codata/universal_constants.hpp>
#include <boost/units/systems/si/codata/electron_constants.hpp>
#include <boost/units/systems/si/codata/electromagnetic_constants.hpp>

namespace monte {

using namespace boost::units::si;

// constants
constexpr double two_pi  = boost::math::constants::two_pi<double>();
constexpr auto   echarge = constants::codata::e;
constexpr auto   dirac   = constants::codata::hbar;
constexpr auto   T       = 300.0 * kelvin;  // room temperature
constexpr auto   cm      = 0.01 * meter;

/*
 * evaluate the constant of proportionality between the acoustic
 * scattering rate and the square root of the electron energy
 */
// components of the scattering rate constant:
constexpr auto   eV = echarge * volt;
constexpr auto   E1 = 7.0 * eV;            // Acoustic deformation potential
constexpr auto   g  = 10e-3 * kilogram;
constexpr auto   cm3 = cm*cm*cm;
constexpr auto   rho = 5.37 * g / cm3;     // Crystal density
constexpr auto   u = 5.2E5 * cm / second;  // speed of sound in GaAs
constexpr double meff = 0.063; // effective electron mass in gamma (000) valley (unitless)

// Shur 2-3-12 as modified by assignment (terms after gamma^1/2 removed)
constexpr auto scatter_const =
    // units for this initial constant not given in Shur :-/
    0.449E18 * ((g / cm3) * (cm / second) * (cm / second) /
                (kelvin * eV * eV * sqrt(eV) * second)) *
    std::pow(meff, 1.5) * T * (E1*E1) /
    (rho * (u*u)) ;

/*
 * evaluate the accelerations due to the constant electric field
 * between collisions
 */
constexpr auto Efield = 10000.0 * volt / cm;
constexpr auto accel_const = echarge * Efield / dirac ;

/*
 * evaluate k-to-velocity conversion constant
 */

// The instantaneous velocity (in real space) is the gradient of E with respect to k
// The Shur book seems to suggest calculating initial and final energies from a given timestep
// then assuming a straight line between. Here I am just using the instantaneous value,
// which is wrong.
constexpr auto   m_e       = constants::codata::m_e;   // electron mass
// If E(k) = (dirac^2 * k^2)/(2 * m_eff), and v = (1/dirac) * dE(k)/dt, then
// v(k) = dirac * k / m_eff
constexpr auto   vel_const = dirac/(meff * m_e) ;

// types

template<typename Float>
struct vector_str {
    vector_str() = default;
    boost::units::quantity<wavenumber, Float> x, y, z;
};

template<typename Float>
using vector_ptr = vector_str<Float> *;

template<typename Float>
using vector_cptr = vector_str<Float> const *;

// utility functions

/* k to energy conversion */
template<typename Float>
boost::units::quantity<energy, Float>
energy_from_k(vector_cptr<Float> kvecptr)
/*
 * This appears to be very complicated and depends on which region
 * of reciprocal lattice space the electron currently occupies.
 * To make this tractable I'm assuming a parabolic relationship
 * given by E(k) = (dirac^2 * k^2)/(2 * m_eff)
 * See Shur p14.
 */
{
    constexpr auto C = (dirac * dirac) / (2.0 * meff * m_e) ;

    boost::units::quantity<energy> nrg =
        C * ((kvecptr->x * kvecptr->x) + (kvecptr->y * kvecptr->y) + (kvecptr->z * kvecptr->z));

    // returning above calculation directly gets a compile error
    // but we can explicitly convert it to another (possibly lower) precision *shrug*
    return boost::units::quantity<energy, Float>{nrg};

}

template<typename Float>
void    coord_convert(vector_cptr<Float> firstkptr,
                      Float theta_r,
                      Float phi_r,
                      vector_ptr<Float> lastkptr)
{
    using namespace boost::units;

    quantity<si::plane_angle, Float>  theta, phi;
    vector_str<Float>                 k2prime;
    quantity<si::wavenumber, Float>   totalk;

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

}

#endif // MONTE_HPP
