/*
 * monte.cpp - monte carlo simulation of acoustic phonon scattering
 */

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "monte.h"

void    init_vector(vector_ptr vecptr);
float   energy_from_k(vector_ptr kvecptr);
void    coord_convert(vector_ptr firstkptr,
                      float theta_r,
                      float phi_r,
                      vector_ptr lastkptr);

vector    kfinal;        /* k vector after collision */
vector    kinit;         /* k vector prior to collision */
float     capgamma;      /* total scattering rate */
int       numtrials;     /* number of scattering events to perform */
float     scatter_const; /* constants of proportionality */
float     accel_const;
float     vel_const;
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
    /*
     * evaluate the constant of proportionality between the acoustic
     * scattering rate and the square root of the electron energy
     */
    scatter_const = (0.449E18)*(0.015813)*(300.0)*(49)/((5.37)*(270.4E9));

    /*
     * evaluate the accelerations due to the constant electric field
     * between collisions
     */
    accel_const = 2.0*3.14159*(1.6E-19)*(10000.0)/(6.63E-34);

    /*
     * evaluate k-to-velocity conversion constant
     */
    vel_const = (10000.0)*(6.63E-34)/((2.0*3.14159)*(0.063)*(9.11E-30));

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

    /* initialize and begin scattering */
    energy = 0.0;
    init_vector(&kinit);
    init_vector(&kfinal);
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

/* k to energy conversion */
float energy_from_k(vector_ptr kvecptr)
/*
 * determine squared velocity, then multiply by effective mass
 * and divide by two, converting to electron volts
 */
{
    float    vx, vy, vz;      /* velocity components */
    float    temp;

    vx = vel_const * kvecptr->x;
    vy = vel_const * kvecptr->y;
    vz = vel_const * kvecptr->z;

    temp = vx*vx + vy*vy + vz*vz;

    temp = (0.5)*(0.063)*(9.11E-30)*temp*(0.0001)/(1.6E-19);

    return temp;
}

void    init_vector(vector_ptr vecptr)
{
    vecptr->x = vecptr->y = vecptr->z = 0.0;
}

void    coord_convert(vector_ptr firstkptr,
                      float theta_r,
                      float phi_r,
                      vector_ptr lastkptr)
{
    float    theta, phi;
    vector   k2prime;
    float    totalk;

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
