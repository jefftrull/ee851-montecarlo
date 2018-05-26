/*
 * monte.cpp - monte carlo simulation of acoustic phonon scattering
 */

#include <vector>
#include <iostream>

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <boost/units/systems/si.hpp>
#include <boost/units/cmath.hpp>
#include <boost/units/io.hpp>
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
    srand48(seed);

    /* initialize and begin scattering */
    cur_trial = 0;

    quantity<wavenumber>    lastkx;
    quantity<velocity>  xvel_avg{0};

    while (cur_trial < numtrials) {
        /* determine the time until the scattering event */
        quantity<si::time> ts = -(1.0/capgamma)*log(drand48());
        lastkx = kinit.x;

        /* accelerate the particle accordingly */
        kinit.x = kinit.x - accel_const * ts;

        /* record the current value of ts and the velocity */
        scat_times[cur_trial] = ts;
        vel_mags[cur_trial] = vel_const * kinit.x;
        cur_trial++;

        /* do average velocity calculations */
        quantity<velocity> cur_avg = vel_const * (lastkx + kinit.x)/2.0;
        xvel_avg = (xvel_avg*(double)(cur_trial - 1) + cur_avg)/(double)cur_trial;

        /* determine the new acoustic scattering rate lambda */
        quantity<energy, float> energy = kinit.get_energy();
        quantity<frequency> lambda = scatter_const * sqrt(energy);
        if (lambda > maxlambda) {
            maxlambda = lambda;              /* possibly useful statistic */
        }

        /* determine if acoustic scattering event occurred */
        if (drand48() >= (lambda/capgamma)) {
            /* no, keep going */
            continue;
        }
        num_real_events++;

        /* determine the angles of the new vector */
        theta =acos(2.0*drand48()-1);
        phi = two_pi * drand48();

        /* determine the resultant vector and replace */
        kinit = kinit.collision_result(theta, phi);
    }

    /* now print out results - good idea to pipe this through more */
    std::cout << num_real_events << " real events out of " << numtrials << ", maximum lambda " << maxlambda << "\n" ;
    std::cout << "average x velocity " << xvel_avg << "\n";
    std::cout << "event\t\tscattering time\t\tvelocity\n";
    for (cur_trial = 0; cur_trial < numtrials; cur_trial++) {
        std::cout << cur_trial << "\t\t" << scat_times[cur_trial] << "\t\t" << vel_mags[cur_trial] << "\n";
    }
}
