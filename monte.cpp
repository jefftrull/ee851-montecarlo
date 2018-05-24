/*
 * monte.cpp - monte carlo simulation of acoustic phonon scattering
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "monte.hpp"

using vector = vector_str<float>;

vector    kfinal;        /* k vector after collision */
vector    kinit;         /* k vector prior to collision */
float     capgamma;      /* total scattering rate */
int       numtrials;     /* number of scattering events to perform */
float     maxlambda;     /* highest lambda ever calculated */
float     xvel_avg;
int       num_real_events;

int main(int argc, char **argv)
    /* initialize, then perform a number of scattering events specified by the
     * user.  Print relevant statistics after the run
     */
{
    (void)argc;          /* unused */
    (void)argv;


    float    energy;     /* temporaries */
    float    ts;
    float    lambda;
    float    phi, theta;
    float    lastkx;
    float    cur_avg;
    float    seed;

    float    *scat_times;    /* output arrays */
    float    *vel_mags;
    int      cur_trial;

    num_real_events = 0;     /* initialize some statistics */
    maxlambda = 0.0;
    xvel_avg = 0.0;

    /* get number of trials from user */
    numtrials = 0;
    int scanf_err = 1;
    while ((numtrials <= 0) || (scanf_err != 1)) {
        printf("Number of scattering events to perform: ");
        scanf_err = scanf("%d", &numtrials);
    }
    scat_times = (float *)malloc(sizeof(float)*numtrials);
    vel_mags   = (float *)malloc(sizeof(float)*numtrials);

    /* get the total scattering rate from the user */
    capgamma = 0.0;
    scanf_err = 1;
    while ((capgamma <= 0.0) || (scanf_err != 1)) {
        printf("Total scattering rate: ");
        scanf_err = scanf("%f", &capgamma);
    }

    /* get a seed for the random number generator */
    seed = 0.0;
    scanf_err = 1;
    while ((seed <= 0.0) || (scanf_err != 1)) {
        printf("Random number generator seed : ");
        scanf_err = scanf("%f", &seed);
    }
    srand48(seed);

/*
    ts = -(1.0/capgamma)*log(drand48());
    long newseed = 1;
    while (fabs(ts - 1.155194e-11) > 1e-18) {
        srand48(newseed++);
        ts = -(1.0/capgamma)*log(drand48());
    }
    printf("timestep %g is within %g of desired value %g\n", ts, 1e-13, 1.155194e-11);
    printf("because their difference is %g\n", ts - 1.155194e-11);

    printf("succeeded with seed = %ld\n", newseed-1);
*/

    /* initialize and begin scattering */
    energy = 0.0;
    cur_trial = 0;

    while (cur_trial < numtrials) {
        /* determine the time until the scattering event */
        ts = -(1.0/capgamma)*log(drand48());
        lastkx = kinit.x;

        /* accelerate the particle accordingly */
        kinit.x = kinit.x - accel_const * ts;

        /* record the current value of ts and the velocity */
        scat_times[cur_trial] = ts;
        vel_mags[cur_trial] = vel_const * kinit.x;
        cur_trial++;

        /* do average velocity calculations */
        cur_avg = vel_const * (lastkx + kinit.x)/(2.0);
        xvel_avg = (xvel_avg*((float)(cur_trial - 1)) + cur_avg)/((float)cur_trial);

        /* determine the new acoustic scattering rate lambda */
        energy = energy_from_k(&kinit);
        lambda = scatter_const * sqrt(energy);
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
        phi = 2.0 * 3.14159*drand48();

        /* determine the resultant vector, kfinal */
        coord_convert(&kinit, theta, phi, &kfinal);
        bcopy(&kfinal, &kinit, sizeof(vector));
    }

    /* now print out results - good idea to pipe this through more */
    printf("%d real events out of %d, maximum lambda %e\n",
           num_real_events, numtrials, maxlambda);
    printf("average x velocity %e\n", xvel_avg);
    printf("event\t\tscattering time\t\tvelocity\n");
    for (cur_trial = 0; cur_trial < numtrials; cur_trial++) {
        printf("%d", cur_trial);
        printf("\t\t%e\t\t%e\n", scat_times[cur_trial], vel_mags[cur_trial]);
    }
}
