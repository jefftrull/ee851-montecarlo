/*
 * monte.cpp - monte carlo simulation of acoustic phonon scattering
 */

/*
 * HISTORY
 * 16-Oct-88 Jeffrey Trull (jt1j) at Carnegie-Mellon University
 *       Created
 * 23-May-18 Jeffrey Trull (edaskel@att.net)
 *       Converted to Modern C++, fixed bugs
 */

// Note to interested readers:
// This program is an updated version of a programming assignment from a semiconductor
// physics class in 1988... The textbook was Michael Shur's "GaAs Devices and Circuits".
// The assignment was to estimate the electron mobility in the material by simulating
// one particular scattering mechanism: acoustic phonons. This was done not because they
// are the dominant mechanism, but because it makes for a tractable problem :)
// The general approach followed is given in Chapter 2 of Shur's textbook; another helpful
// source I (in retrospect) found is:
// Jacoboni and Reggiani: Monte Carlo Method in Transport (Rev Mod Phys, July 1983)
// which gives a more detailed explanation/motivation of the process


#include <vector>
#include <iostream>
#include <random>

#include <cmath>

#include <boost/units/systems/si.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>

#include "monte.hpp"

using namespace boost::units;
using namespace boost::units::si;
using namespace monte;

using vector = vector_str<float>;

vector                kinit;         /* k vector prior to collision */
float                 capgamma_f;    /* total scattering rate */
int                   numtrials;     /* number of scattering events to perform */
quantity<frequency>   maxlambda;     /* highest lambda ever calculated */
int                   num_real_events;

int main(int argc, char **argv)
    /* initialize, then perform a number of scattering events specified by the
     * user.  Print relevant statistics after the run
     */
{
    (void)argc;          /* unused */
    (void)argv;


    float    phi, theta;
    float    seed;

    int      cur_trial;

    num_real_events = 0;     /* initialize some statistics */
    maxlambda = quantity<frequency>{0};

    /* get number of trials from user */
    numtrials = 0;
    int scanf_err = 1;
    while ((numtrials <= 0) || (scanf_err != 1)) {
        printf("Number of scattering events to perform: ");
        scanf_err = scanf("%d", &numtrials);
    }
    std::vector<quantity<si::time, float>> scat_times(numtrials);
    std::vector<quantity<si::velocity, float>> vel_mags(numtrials);

    /* get the total scattering rate from the user */
    capgamma_f = 0.0;
    scanf_err = 1;
    while ((capgamma_f <= 0.0) || (scanf_err != 1)) {
        printf("Total scattering rate: ");
        scanf_err = scanf("%f", &capgamma_f);
    }
    quantity<frequency> capgamma = capgamma_f * hertz;

    /* get a seed for the random number generator */
    seed = 0.0;
    scanf_err = 1;
    while ((seed <= 0.0) || (scanf_err != 1)) {
        printf("Random number generator seed : ");
        scanf_err = scanf("%f", &seed);
    }
    std::mt19937 randeng(seed);
    std::uniform_real_distribution<double> drand_dist(0, 1.0);
    auto drand = [&]() { return drand_dist(randeng); };

    /* initialize and begin scattering */
    cur_trial = 0;

    quantity<wavenumber>    lastkx;

    // Use Boost.Accumulators to produce an average velocity
    using namespace boost::accumulators;
    accumulator_set<quantity<velocity>,      // what we are storing
                    stats<tag::mean>,        // what we want to calculate
                    quantity<si::time>>      // the weight for each data point
        vel_acc;

    while (cur_trial < numtrials) {
        /* determine the time until the scattering event */
        quantity<si::time> ts = -(1.0/capgamma)*log(drand());
        lastkx = kinit.x;

        /* accelerate the particle accordingly */
        kinit.x = kinit.x - accel_const * ts;

        /* record the current value of ts and the velocity */
        scat_times[cur_trial] = ts;
        vel_mags[cur_trial] = vel_const * kinit.x;
        cur_trial++;

        /* do average velocity calculations */
        quantity<velocity> cur_avg = vel_const * (lastkx + kinit.x)/2.0;
        vel_acc(cur_avg, weight = ts);

        /* determine the new acoustic scattering rate lambda */
        quantity<energy, float> energy = kinit.get_energy();
        quantity<frequency> lambda = scatter_const * sqrt(energy);
        if (lambda > maxlambda) {
            maxlambda = lambda;              /* possibly useful statistic */
        }

        /* determine if acoustic scattering event occurred */
        if (drand() >= (lambda/capgamma)) {
            /* no, keep going */
            // this is a "self-scattering" event, a fake event that makes inverting the simulation
            // (choosing a flight time at random instead of using fixed intervals and testing
            // each one) work out right
            continue;
        }
        num_real_events++;

        /* determine the angles of the new vector */

        // If I understand this correctly we are calculating the angle theta between
        // the original k (taken facing along the Z axis) and k'
        // theta is the polar angle (declination from Z ("up"))
        // Taking dP(theta)d(theta) ~ sin(theta)d(theta), we need integration constants
        // that give the right total from 0 to pi and the right values at the extremes
        // (0 at 0, 1 at pi). P(theta) = 0.5*(1-cos(theta)) fits.
        theta = acos(1 - 2.0*drand());
        // phi, the "azimuthal angle" (from x in the direction of y) is uniformly distributed
        phi = two_pi * drand();
        // This is consistent with the approach in "Bulk Monte Carlo Method Described",
        // https://nanohub.org/resources/4844/download/montecarlocodedescribed.pdf
        // Should I have used the more complex procedure described in Shur p23?

        /* determine the resultant vector and replace */
        kinit = kinit.collision_result(theta, phi);
    }

    /* now print out results - good idea to pipe this through more */
    std::cout << num_real_events << " real events out of " << numtrials << ", maximum lambda " << maxlambda << "\n" ;
    std::cout << "average x velocity " << no_prefix << mean(vel_acc) << "\n";
    std::cout << "estimated mobility: " << ((mean(vel_acc) / Efield) / ((cm*cm)/(volt*second))).value() << " cm^2/Vs\n";
    std::cout << "event\t\tscattering time\t\tvelocity\n";
    for (cur_trial = 0; cur_trial < numtrials; cur_trial++) {
        std::cout << cur_trial << "\t\t" << engineering_prefix << scat_times[cur_trial] << "\t\t" << no_prefix << vel_mags[cur_trial] << "\n";
    }
}
