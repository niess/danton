/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C99 executable dedicated to the sampling of decaying
 * taus from ultra high energy neutrinos.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

/* Standard library includes. */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Linux / POSIX includes. */
#include <getopt.h>

/* The Physics APIs. */
#include "alouette.h"
#include "ent.h"
#include "pumas.h"

/* The spherical Earth radius, in m. */
#define EARTH_RADIUS 6371.E+03

/* The radius of the geostationary orbit, in m. */
#define GEO_ORBIT 42164E+03

/* Avogadro's number. */
#define PHYS_NA 6.022E+23

/* Biasing factor for backward tau decays. */
#define DECAY_BIAS 6.

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Handles for the transport engines. */
static struct ent_physics * physics = NULL;
static struct pumas_context * ctx_pumas = NULL;

/* The tau lepton mass, in GeV / c^2. */
static double tau_mass;

/* The tau lifetime, in m / c. */
static double tau_ctau0;

/* Finalise and exit to the OS. */
static void gracefully_exit(int rc)
{
        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        pumas_context_destroy(&ctx_pumas);
        pumas_finalise();
        alouette_finalise();
        exit(rc);
}

/* Show help and exit. */
static void exit_with_help(int code)
{
        // clang-format off
        fprintf(stderr,
"Usage: danton [OPTION]... [PID]\n"
"Simulate taus decaying in the Earth atmosphere, originating from neutrinos\n"
"with the given flavour (PID).\n"
"\n"
"Configuration options:\n"
"      --altitude=Z           switch to a point estimate at an altitude of Z\n"
"      --altitude-max=Z       set the maximum decay altitude to Z [1E+05]\n"
"      --altitude-min=Z       set the minimum decay altitude to Z [1E+00]\n"
"  -c, --cos-theta=C          switch to a point estimate with cos(theta)=C\n"
"      --cos-theta-max=C      set the maximum value of cos(theta) to C [0.25]\n"
"      --cos-theta-min=C      set the minimum value of cos(theta) to C [0.15]\n"
"      --elevation=A          switch to a point estimate at an elevation\n"
"                               angle of A\n"
"      --elevation-max=A      set the maximum elevation angle to A [10.]\n"
"      --elevation-min=A      set the minimum elevation angle to A [-10.]\n"
"  -e, --energy=E             switch to a monokinetic beam of primaries with\n"
"                               energy E\n"
"      --energy-analog        switch to an analog sampling of the primary\n"
"                               energy\n"
"      --energy-cut=E         set to E the energy under which all particles\n"
"                               are killed [1E+03]\n"
"      --energy-max=E         set to E the maximum primary energy [1E+12]\n"
"      --energy-min=E         set to E the minimum primary energy [1E+07]\n"
"      --forward              switch to forward Monte-Carlo\n"
"      --longitudinal         disable the transverse transport\n"
"  -n                         set the number of incoming neutrinos [10000] or\n"
"                               the number of bins for the grammage sampling\n"
"                               [1001]\n"
"      --pem-no-sea           Replace seas with rocks in the Preliminary\n"
"                               Earth Model (PEM)\n"
"  -t, --taus=M               request at most M > 0 tau decays\n"
"\n"
"Control options:\n"
"      --append               append to the output file if it already exists\n"
"      --grammage             disable neutrino interactions but compute the\n"
"                               total grammage along the path.\n"
"      --flux=Z               switch to the sampling of the tau flux at an\n"
"                               altitude of Z\n"
"  -h, --help                 show this help and exit\n"
"  -o, --output-file=FILE     set the output FILE for the computed flux\n"
"      --pdf-file=FILE        specify a pdf FILE for partons distributions\n"
"                               [(builtin)/CT14nnlo_0000.dat]. The format\n"
"                               must be `lhagrid1` complient.\n"
"\n"
"Energies (E) must be given in GeV, altitudes (Z) in m and angles (A) in deg.\n"
"PID must be one of 15 (tau) or -15 (tau_bar) in backward mode and -12 \n"
"(nu_e_bar), 16 (nu_tau) or -16 (nu_tau_bar) in forward mode.\n"
"\n"
"The default behaviour is to randomise the primary neutrino energy over a\n"
"1/E^2 spectrum with a log bias. Use the --energy-analog option in order to\n"
"to an analog simulation. The primary direction is randomised uniformly over\n"
"a solid angle specified by the min an max value of the cosine of the angle\n"
"theta with the local vertical at the atmosphere entrance. Use the -c,\n"
"--cos-theta option in order to perform a point estimate instead.\n"
"\n"
"Note that altitudes (Z) must be strictly positive.\n");
        exit(code);
        // clang-format on
}

/* Prototype for long options setters. */
typedef void opt_setter_t(char * optarg, void * argument, void * variable);

/* Container for setting a long option. */
struct opt_setter {
        void * variable;
        opt_setter_t * set;
        void * argument;
};

/* String copy setter. */
static void opt_copy(char * optarg, void * argument, void * variable)
{
        char ** s = variable;
        *s = optarg;
}

/* Setter for a double variable. */
static void opt_strtod(char * optarg, void * argument, void * variable)
{
        double * d = variable;
        *d = strtod(optarg, argument);
}

/* Error handler for ENT. */
static void handle_ent(enum ent_return rc, ent_function_t * caller)
{
        /* Dump an error message. */
        ent_error_print(stderr, rc, caller, "\t", "\n");
        fprintf(stderr, "\n");

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_FAILURE);
}

/* Error handler for PUMAS. */
static void handle_pumas(
    enum pumas_return rc, pumas_function_t * caller, struct pumas_error * error)
{
        /* Dump the error summary. */
        fprintf(stderr, "error : ");
        pumas_error_print(stderr, rc, caller, error);
        fprintf(stderr, "\n");

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_FAILURE);
}

/* Density according to the Preliminary Earth Model (PEM). */
static double pem_model0(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a2 = -8.8381E+03;
        *density = 13.0885E+03 + a2 * x * x;
        const double xg = (x <= 5E-02) ? 5E-02 : x;
        return 0.01 * EARTH_RADIUS / fabs(2. * a2 * xg);
}

static double pem_model1(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 1.2638E+03;
        *density = 12.58155E+03 + x * (-a + x * (-3.6426E+03 - x * 5.5281E+03));
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model2(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 6.4761E+03;
        *density = 7.9565E+03 + x * (-a + x * (2.5283E+03 - x * 3.0807E+03));
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model3(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 1.4836E+03;
        *density = 5.3197E+03 - a * x;
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model4(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 8.0298E+03;
        *density = 11.2494E+03 - a * x;
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model5(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 3.8045E+03;
        *density = 7.1089E+03 - a * x;
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model6(double r, double * density)
{
        const double x = r / EARTH_RADIUS;
        const double a = 0.6924E+03;
        *density = 2.691E+03 + a * x;
        return 0.01 * EARTH_RADIUS / a;
}

static double pem_model7(double r, double * density)
{
        *density = 2.9E+03;
        return 0.;
}

static double pem_model8(double r, double * density)
{
        *density = 2.6E+03;
        return 0.;
}

static double pem_model9(double r, double * density)
{
        *density = 1.02E+03;
        return 0.;
}

/* The U.S. standard atmosphere model. */
#define USS_MODEL(INDEX, B, C)                                                 \
        static double uss_model##INDEX(double r, double * density)             \
        {                                                                      \
                *density = B / C * exp(-(r - EARTH_RADIUS) / C);               \
                return 0.01 * C;                                               \
        }

USS_MODEL(0, 12226.562, 9941.8638)
USS_MODEL(1, 11449.069, 8781.5355)
USS_MODEL(2, 13055.948, 6361.4304)
USS_MODEL(3, 5401.778, 7721.7016)

/* Outer space density model. */
static double space_model0(double r, double * density)
{
        *density = 1.E-21; /* ~10^6 H per m^-3. */
        return 0.;
}

/* Generic Monte-Carlo state with stepping data. */
struct generic_state {
        union {
                struct ent_state ent;
                struct pumas_state pumas;
        } base;
        double r;
        int is_tau;
        int is_inside;
        int has_crossed;
};

/* Density callbacks for ENT. */
#define DENSITY(MODEL, INDEX)                                                  \
        static double density_##MODEL##INDEX(struct ent_medium * medium,       \
            struct ent_state * state, double * density)                        \
        {                                                                      \
                struct generic_state * s = (struct generic_state *)state;      \
                return MODEL##_model##INDEX(s->r, density);                    \
        }

DENSITY(pem, 0)
DENSITY(pem, 1)
DENSITY(pem, 2)
DENSITY(pem, 3)
DENSITY(pem, 4)
DENSITY(pem, 5)
DENSITY(pem, 6)
DENSITY(pem, 7)
DENSITY(pem, 8)
DENSITY(pem, 9)
DENSITY(uss, 0)
DENSITY(uss, 1)
DENSITY(uss, 2)
DENSITY(uss, 3)
DENSITY(space, 0)

/* Local callbacks for PUMAS. */
#define LOCALS(MODEL, INDEX)                                                   \
        static double locals_##MODEL##INDEX(struct pumas_medium * medium,      \
            struct pumas_state * state, struct pumas_locals * locals)          \
        {                                                                      \
                struct generic_state * s = (struct generic_state *)state;      \
                const double step =                                            \
                    MODEL##_model##INDEX(s->r, &locals->density);              \
                memset(locals->magnet, 0x0, sizeof(locals->magnet));           \
                return step;                                                   \
        }

LOCALS(pem, 0)
LOCALS(pem, 1)
LOCALS(pem, 2)
LOCALS(pem, 3)
LOCALS(pem, 4)
LOCALS(pem, 5)
LOCALS(pem, 6)
LOCALS(pem, 7)
LOCALS(pem, 8)
LOCALS(pem, 9)
LOCALS(uss, 0)
LOCALS(uss, 1)
LOCALS(uss, 2)
LOCALS(uss, 3)
LOCALS(space, 0)

/* Flag for switching to a tau flux computation. */
static int flux_mode = 0;

/* Altitude at which the flux is to be computed. */
static double flux_altitude = 0.;

/* Generic medium callback. */
static double medium(const double * position, const double * direction,
    int * index, struct generic_state * state)
{
        *index = -1;
        double step = 0.;

        const double r2 = position[0] * position[0] +
            position[1] * position[1] + position[2] * position[2];
        if (r2 > GEO_ORBIT * GEO_ORBIT) return step;
        const double r = sqrt(r2);
        state->r = r;

        if (flux_mode && state->is_tau && state->has_crossed >= 0) {
                /* Check the flux boundary in forward MC. */
                const double rf = EARTH_RADIUS + flux_altitude;
                if (state->is_inside < 0)
                        state->is_inside = (r < rf) ? 1 : 0;
                else if ((state->is_inside && (r >= rf)) ||
                    (!state->is_inside && (r <= rf))) {
                        state->has_crossed = 1;
                        return step;
                }
        };

        const double ri[16] = { 1221.5E+03, 3480.E+03, 5701.E+03, 5771.E+03,
                5971.E+03, 6151.E+03, 6346.6E+03, 6356.E+03, 6368.E+03,
                EARTH_RADIUS, EARTH_RADIUS + 4.E+03, EARTH_RADIUS + 1.E+04,
                EARTH_RADIUS + 4.E+04, EARTH_RADIUS + 1.E+05, GEO_ORBIT,
                2 * GEO_ORBIT };

        /* Kill neutrinos that exit the atmosphere.  */
        if ((!state->is_tau) && (state->r > ri[13])) return step;

        int i;
        for (i = 0; i < sizeof(ri) / sizeof(*ri) - 1; i++) {
                if (r <= ri[i]) {
                        *index = i;

                        /* Outgoing intersection. */
                        const double b = position[0] * direction[0] +
                            position[1] * direction[1] +
                            position[2] * direction[2];
                        const double r1 = ri[i + 1];
                        const double d2 = b * b + r1 * r1 - r * r;
                        const double d = (d2 <= 0.) ? 0. : sqrt(d2);
                        step = d - b - 1.;

                        if ((i > 0) && (b < 0.)) {
                                /* This is a downgoing trajectory. First, let
                                 * us compute the intersection with the lower
                                 * radius.
                                 */
                                const double r1 = ri[i - 1];
                                const double d2 = b * b + r1 * r1 - r * r;
                                if (d2 > 0.) {
                                        const double d = sqrt(d2);
                                        const double s = d - b - 1.;
                                        if (s < step) step = s;
                                }

                                if (i > 1) {
                                        /* Let us check for an intersection with
                                         * the below lower radius.
                                         */
                                        const double r1 = ri[i - 2];
                                        const double d2 =
                                            b * b + r1 * r1 - r * r;
                                        if (d2 > 0.) {
                                                const double d = sqrt(d2);
                                                const double s = d - b - 1.;
                                                if (s < step) step = s;
                                        }
                                }
                        }
                        if (step < 1.) step = 1.;
                        break;
                }
        }
        return step;
}

/* Media table for ENT. */
#define ZR 13.
#define AR 26.
#define ZW 10.
#define AW 18.
#define ZA 7.32
#define AA 14.72

static struct ent_medium media_ent[] = { { ZR, AR, &density_pem0 },
        { ZR, AR, &density_pem1 }, { ZR, AR, &density_pem2 },
        { ZR, AR, &density_pem3 }, { ZR, AR, &density_pem4 },
        { ZR, AR, &density_pem5 }, { ZR, AR, &density_pem6 },
        { ZR, AR, &density_pem7 }, { ZR, AR, &density_pem8 },
        { ZW, AW, &density_pem9 }, { ZA, AA, &density_uss0 },
        { ZA, AA, &density_uss1 }, { ZA, AA, &density_uss2 },
        { ZA, AA, &density_uss3 }, { ZA, AA, &density_space0 } };

/* Medium callback encapsulation for ENT. */
static double medium_ent(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium_ptr)
{
        int index;
        const double step = medium(state->position, state->direction, &index,
            (struct generic_state *)state);
        if (index >= 0)
                *medium_ptr = media_ent + index;
        else
                *medium_ptr = NULL;
        return step;
}

#undef ZR
#undef AR
#undef ZA
#undef AA

/* Media table for PUMAS. */
static struct pumas_medium media_pumas[] = { { 0, &locals_pem0 },
        { 0, &locals_pem1 }, { 0, &locals_pem2 }, { 0, &locals_pem3 },
        { 0, &locals_pem4 }, { 0, &locals_pem5 }, { 0, &locals_pem6 },
        { 0, &locals_pem7 }, { 0, &locals_pem8 }, { 1, &locals_pem9 },
        { 2, &locals_uss0 }, { 2, &locals_uss1 }, { 2, &locals_uss2 },
        { 2, &locals_uss3 }, { 2, &locals_space0 } };

/* Medium callback encapsulation for PUMAS. */
double medium_pumas(struct pumas_context * context, struct pumas_state * state,
    struct pumas_medium ** medium_ptr)
{
        int index;
        const double step = medium(state->position, state->direction, &index,
            (struct generic_state *)state);
        if (index >= 0)
                *medium_ptr = media_pumas + index;
        else
                *medium_ptr = NULL;
        return step;
}

/* Data for the Mersenne Twister PRNG. */
static struct {
#define MT_PERIOD 624
        int index;
        unsigned long data[MT_PERIOD];
} random_mt;

/* Initialise the PRNG random seed. */
static void random_initialise()
{
        /* Get a seed from /dev/urandom*/
        unsigned long seed;
        static const char * urandom = "/dev/urandom";
        FILE * fp = fopen(urandom, "rb");
        if (fp == NULL) goto error;
        if (fread(&seed, sizeof(long), 1, fp) <= 0) {
                fclose(fp);
                goto error;
        }
        fclose(fp);

        /* Set the Mersenne Twister initial state. */
        random_mt.data[0] = seed & 0xffffffffUL;
        int j;
        for (j = 1; j < MT_PERIOD; j++) {
                random_mt.data[j] =
                    (1812433253UL * (random_mt.data[j - 1] ^
                                        (random_mt.data[j - 1] >> 30)) +
                        j);
                random_mt.data[j] &= 0xffffffffUL;
        }
        random_mt.index = MT_PERIOD;

        return;
error:
        fprintf(
            stderr, "danton: could not Initialise PRNG from %s.\n", urandom);
        gracefully_exit(EXIT_FAILURE);
}

/* Uniform pseudo random distribution over [0,1] from a Mersenne Twister. */
static double random_uniform01(void * context)
{
        /* Check the buffer. */
        if (random_mt.index < MT_PERIOD - 1) {
                random_mt.index++;
        } else {
                /* Update the MT state. */
                const int M = 397;
                const unsigned long UPPER_MASK = 0x80000000UL;
                const unsigned long LOWER_MASK = 0x7fffffffUL;
                static unsigned long mag01[2] = { 0x0UL, 0x9908b0dfUL };
                unsigned long y;
                int kk;
                for (kk = 0; kk < MT_PERIOD - M; kk++) {
                        y = (random_mt.data[kk] & UPPER_MASK) |
                            (random_mt.data[kk + 1] & LOWER_MASK);
                        random_mt.data[kk] = random_mt.data[kk + M] ^ (y >> 1) ^
                            mag01[y & 0x1UL];
                }
                for (; kk < MT_PERIOD - 1; kk++) {
                        y = (random_mt.data[kk] & UPPER_MASK) |
                            (random_mt.data[kk + 1] & LOWER_MASK);
                        random_mt.data[kk] =
                            random_mt.data[kk + (M - MT_PERIOD)] ^ (y >> 1) ^
                            mag01[y & 0x1UL];
                }
                y = (random_mt.data[MT_PERIOD - 1] & UPPER_MASK) |
                    (random_mt.data[0] & LOWER_MASK);
                random_mt.data[MT_PERIOD - 1] =
                    random_mt.data[M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];
                random_mt.index = 0;
        }

        /* Tempering. */
        unsigned long y = random_mt.data[random_mt.index];
        y ^= (y >> 11);
        y ^= (y << 7) & 0x9d2c5680UL;
        y ^= (y << 15) & 0xefc60000UL;
        y ^= (y >> 18);

        /* Convert to a floating point and return. */
        return y * (1.0 / 4294967295.0);
}

/* Loader for PUMAS. */
void load_pumas()
{
        const enum pumas_particle particle = PUMAS_PARTICLE_TAU;
        const char * dump = "materials.b";

        /* First, attempt to load any binary dump. */
        FILE * stream = fopen(dump, "rb");
        if (stream != NULL) {
                pumas_error_catch(1); /* catch any error. */
                pumas_load(stream);
                fclose(stream);
                pumas_error_raise();
                pumas_particle(NULL, &tau_ctau0, &tau_mass);
                return;
        }

        /* If no binary dump, initialise from the MDF and dump. */
        pumas_initialise(particle, NULL, NULL, NULL);
        pumas_particle(NULL, NULL, &tau_mass);

        /* Dump the library configuration. */
        stream = fopen(dump, "wb+");
        if (stream == NULL) handle_pumas(PUMAS_RETURN_IO_ERROR, NULL, NULL);
        pumas_dump(stream);
        fclose(stream);
}

/* Flag to disable the transverse transport. */
static int longitudinal = 0;

/* Set a neutrino state from a tau decay product. */
static void copy_neutrino(const struct pumas_state * tau, int pid,
    const double * momentum, struct ent_state * neutrino, double * direction)
{
        neutrino->pid = pid;
        neutrino->energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2]);
        memcpy(neutrino->position, tau->position, sizeof(neutrino->position));
        if (longitudinal) {
                neutrino->direction[0] = direction[0];
                neutrino->direction[1] = direction[1];
                neutrino->direction[2] = direction[2];
        } else {
                neutrino->direction[0] = momentum[0] / neutrino->energy;
                neutrino->direction[1] = momentum[1] / neutrino->energy;
                neutrino->direction[2] = momentum[2] / neutrino->energy;
        }
        neutrino->distance = tau->distance;
        neutrino->grammage = tau->grammage;
        neutrino->weight = tau->weight;
}

/* Output file name for the results. */
static char * output_file = NULL;

/* Open the output stream. */
static FILE * output_open(void)
{
        if (output_file == NULL) return stdout;
        return fopen(output_file, "a");
}

/* Close the output stream. */
static void output_close(FILE * stream)
{
        if (output_file != NULL) fclose(stream);
}

/* Utility functions for formating the results to the ouytput stream. */
static void format_ancester(long eventid, const struct ent_state * ancester)
{
        FILE * stream = output_open();
        fprintf(stream, "%10ld %4d %12.5lE %12.5lE %12.5lE %12.5lE %13.3lf "
                        "%13.3lf %13.3lf %12.5lE\n",
            eventid + 1, ancester->pid, ancester->energy,
            ancester->direction[0], ancester->direction[1],
            ancester->direction[2], ancester->position[0],
            ancester->position[1], ancester->position[2], ancester->weight);
        output_close(stream);
}

static void format_tau(int generation, int pid,
    const struct pumas_state * production, const struct pumas_state * decay)
{
        FILE * stream = output_open();
        fprintf(stream, "%10d %4d %12.5lE %12.5lE %12.5lE %12.5lE %13.3lf "
                        "%13.3lf %13.3lf\n",
            generation, pid, production->kinetic, production->direction[0],
            production->direction[1], production->direction[2],
            production->position[0], production->position[1],
            production->position[2]);
        fprintf(stream, "%10c %4c %12.5lE %12.5lE %12.5lE %12.5lE %13.3lf "
                        "%13.3lf %13.3lf\n",
            ' ', ' ', decay->kinetic, decay->direction[0], decay->direction[1],
            decay->direction[2], decay->position[0], decay->position[1],
            decay->position[2]);
        output_close(stream);
}

static void format_decay_product(int pid, const double momentum[3])
{
        FILE * stream = output_open();
        fprintf(stream, "%10c %4d %12c %12.5lE %12.5lE %12.5lE\n", ' ', pid,
            ' ', momentum[0], momentum[1], momentum[2]);
        output_close(stream);
}

static void format_grammage(double cos_theta, double grammage)
{
        FILE * stream = output_open();
        fprintf(stream, "%12.5lE %12.5lE\n", cos_theta, grammage);
        output_close(stream);
}

/* Print the header of the output file. */
static void print_header_decay(FILE * stream)
{
        fprintf(stream,
            "    Event   PID    Energy             Direction or "
            "Momentum                         Position               "
            "      Weight\n                    (GeV)                "
            " (1 or GeV/c)                                 (m)\n     "
            "                           ux or Px     uy or Py    "
            "uz or Pz         X             Y             Z\n");
}

static void print_header_grammage(FILE * stream)
{
        fprintf(stream, "  cos(theta)    Grammage\n                (kg/m^2)\n");
}

/* The requested number of tau decays. */
static int n_taus = 0;

/* The lower (upper) energy bound under (above) which all particles are
 * killed.
 */
static double energy_cut = -1.;

/* Forward transport routine, recursive. */
static void transport_forward(struct ent_context * ctx_ent,
    struct ent_state * neutrino, long eventid, int generation,
    int primary_dumped, const struct ent_state * ancester, int * done)
{
        if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
            (abs(neutrino->pid) != ENT_PID_NU_TAU))
                return;
        if ((n_taus > 0) && (*done >= n_taus)) return;

        /* Backup the initial direction if the transverse transport is
         * disabled.
         */
        double direction[3];
        if (longitudinal)
                memcpy(direction, neutrino->direction, sizeof(direction));

        struct ent_state product;
        enum ent_event event;
        for (;;) {
                /* Neutrino transport with ENT. */
                ent_transport(physics, ctx_ent, neutrino, &product, &event);
                if ((event == ENT_EVENT_EXIT) ||
                    (neutrino->energy <= energy_cut + FLT_EPSILON))
                        break;
                if (longitudinal) {
                        memcpy(neutrino->direction, direction,
                            sizeof(neutrino->direction));
                        memcpy(product.direction, direction,
                            sizeof(product.direction));
                }
                if (abs(neutrino->pid) == ENT_PID_TAU) {
                        /* Exchange the lepton and product state. */
                        struct ent_state tmp;
                        memcpy(&tmp, neutrino, sizeof(tmp));
                        memcpy(neutrino, &product, sizeof(*neutrino));
                        memcpy(&product, &tmp, sizeof(product));
                }

                if (abs(product.pid) == ENT_PID_TAU) {
                        /* Tau transport with PUMAS. */
                        const double charge = (product.pid > 0) ? -1. : 1.;
                        const double kinetic = product.energy - tau_mass;
                        struct generic_state tau_data = {
                                .base.pumas = { charge, kinetic,
                                    product.distance, product.grammage, 0.,
                                    product.weight, .decayed = 0 },
                                .r = 0.,
                                .is_tau = 1,
                                .is_inside = -1
                        };
                        struct pumas_state * tau = &tau_data.base.pumas;
                        memcpy(&tau->position, &product.position,
                            sizeof(tau->position));
                        memcpy(&tau->direction, &product.direction,
                            sizeof(tau->direction));
                        struct pumas_state tau_prod;
                        memcpy(&tau_prod, tau, sizeof(tau_prod));
                        pumas_transport(ctx_pumas, tau);
                        if (tau->decayed) {
                                /* Tau decay with ALOUETTE/TAUOLA. */
                                const double p = sqrt(tau->kinetic *
                                    (tau->kinetic + 2. * tau_mass));
                                double momentum[3] = { p * tau->direction[0],
                                        p * tau->direction[1],
                                        p * tau->direction[2] };
                                int trials;
                                for (trials = 0; trials < 20; trials++) {
                                        if (alouette_decay(product.pid,
                                                momentum, tau->direction) ==
                                            ALOUETTE_RETURN_SUCCESS)
                                                break;
                                }
                                int pid, nprod = 0;
                                struct generic_state nu_e_data, nu_t_data;
                                memset(&nu_e_data, 0x0, sizeof(nu_e_data));
                                memset(&nu_t_data, 0x0, sizeof(nu_t_data));
                                struct ent_state *nu_e = NULL, *nu_t = NULL;
                                while (alouette_product(&pid, momentum) ==
                                    ALOUETTE_RETURN_SUCCESS) {
                                        if (abs(pid) == 16) {
                                                /* Update the neutrino state
                                                 * with the nu_tau daughter.
                                                 */
                                                nu_t = &nu_t_data.base.ent;
                                                copy_neutrino(tau, pid,
                                                    momentum, nu_t, direction);
                                                continue;
                                        } else if (pid == -12) {
                                                nu_e = &nu_e_data.base.ent;
                                                copy_neutrino(tau, pid,
                                                    momentum, nu_e, direction);
                                                continue;
                                        } else if (flux_mode || (pid == 12) ||
                                            (abs(pid) == 13) ||
                                            (abs(pid) == 14))
                                                continue;

                                        /* Log the decay if in air. */
                                        int medium_index;
                                        medium(tau->position, tau->direction,
                                            &medium_index, &tau_data);
                                        if (medium_index < 10) continue;
                                        if (nprod == 0) {
                                                if (!primary_dumped) {
                                                        format_ancester(
                                                            eventid, ancester);
                                                        primary_dumped = 1;
                                                }
                                                format_tau(generation,
                                                    product.pid, &tau_prod,
                                                    tau);
                                        }
                                        format_decay_product(pid, momentum);
                                        nprod++;
                                }
                                if (nprod > 0) (*done)++;
                                if ((n_taus > 0) && (*done >= n_taus)) return;
                                generation++;

                                /* Process any additional nu_e~ or nu_tau. */
                                if (nu_e != NULL) {
                                        transport_forward(ctx_ent, nu_e,
                                            eventid, generation, primary_dumped,
                                            ancester, done);
                                }
                                if (nu_t != NULL)
                                        transport_forward(ctx_ent, nu_t,
                                            eventid, generation, primary_dumped,
                                            ancester, done);
                        } else if (tau_data.has_crossed) {
                                if (!primary_dumped) {
                                        format_ancester(eventid, ancester);
                                        primary_dumped = 1;
                                }
                                format_tau(
                                    generation, product.pid, &tau_prod, tau);
                        }
                }
                if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
                    (abs(neutrino->pid) != ENT_PID_NU_TAU))
                        break;
        }
}

/* Ancester callback for ENT. */
static double ancester_cb(struct ent_context * context, enum ent_pid ancester,
    struct ent_state * daughter)
{
        if (daughter->pid == ENT_PID_NU_BAR_E) {
                if (ancester == ENT_PID_NU_BAR_E) return 1.;
                return 0.;
        } else if (daughter->pid == ENT_PID_NU_TAU) {
                if (ancester == ENT_PID_NU_TAU)
                        return 1.;
                else if (ancester == ENT_PID_TAU)
                        return 1E-03;
                return 0.;
        } else if (daughter->pid == ENT_PID_NU_BAR_TAU) {
                if (ancester == ENT_PID_NU_BAR_TAU)
                        return 1.;
                else if (ancester == ENT_PID_TAU_BAR)
                        return 1E-03;
                return 0.;
        } else if (daughter->pid == ENT_PID_TAU) {
                if (ancester == ENT_PID_NU_TAU) return 1.;
                return 0.;
        } else if (daughter->pid == ENT_PID_TAU_BAR) {
                if (ancester == ENT_PID_NU_BAR_TAU) return 1.;
                return 0.;
        } else
                return 0.;
}

/* Polarisation callback for ALOUETTE. A 100% longitudinal polarisation is
 * assumed.
 */
void polarisation_cb(int pid, const double momentum[3], double * polarisation)
{
        double nrm = momentum[0] * momentum[0] + momentum[1] * momentum[1] +
            momentum[2] * momentum[2];
        if (nrm <= 0.) {
                polarisation[0] = 0.;
                polarisation[1] = 0.;
                polarisation[2] = 0.;
                return;
        }
        nrm = 1. / sqrt(nrm);
        polarisation[0] = momentum[0] * nrm;
        polarisation[1] = momentum[1] * nrm;
        polarisation[2] = momentum[2] * nrm;
}

/* Backward transport routine. */
static void transport_backward(struct ent_context * ctx_ent,
    struct pumas_state * tau, long eventid, int generation,
    struct pumas_state * tau_at_decay, struct pumas_state * tau_at_production,
    int * done)
{
        /* Backup the final tau state at decay. */
        if (generation == 1) memcpy(tau_at_decay, tau, sizeof(*tau_at_decay));

        /* Backward propagate the tau state. */
        const double Kf = tau->kinetic;
        pumas_transport(ctx_pumas, tau);
        if (!tau->decayed ||
            (tau->kinetic + tau_mass >= energy_cut - FLT_EPSILON))
                return;

        /* Apply the BMC weight for the tau decay. */
        const double Pi = sqrt(tau->kinetic * (tau->kinetic + 2. * tau_mass));
        tau->weight *= tau_ctau0 * Pi / tau_mass;
        if (!flux_mode || (generation > 1)) {
                const double Pf = sqrt(Kf * (Kf + 2. * tau_mass));
                tau->weight /= tau_ctau0 * Pf / tau_mass;
        }

        /* Backup the tau state at production. */
        if (generation == 1)
                memcpy(tau_at_production, tau, sizeof(*tau_at_production));

        /* Backward generate the production vertex. */
        enum ent_pid pid = (tau->charge < 0.) ? ENT_PID_TAU : ENT_PID_TAU_BAR;
        const double * const r = tau->position;
        const double * const u = tau->direction;
        struct generic_state g_state = {
                .base.ent = { pid, tau->kinetic + tau_mass, tau->distance,
                    tau->grammage, tau->weight, { r[0], r[1], r[2] },
                    { u[0], u[1], u[2] } },
                .r = 0.,
                .is_tau = 0,
                .is_inside = -1,
                .has_crossed = -1
        };
        struct ent_state * state = &g_state.base.ent;
        struct ent_medium * medium;
        medium_ent(ctx_ent, state, &medium);
        if (medium == NULL) return;
        ent_vertex(physics, ctx_ent, state, medium, ENT_PROCESS_NONE, NULL);

        /* Append the effective BMC weight for the transport, in order to
         * recover a flux convention.
         */
        double cs;
        ent_physics_cross_section(physics, pid, state->energy, medium->Z,
            medium->A, ENT_PROCESS_NONE, &cs);
        double density;
        medium->density(medium, state, &density);
        state->weight *= 1E+03 * cs * PHYS_NA * density / medium->A;

        /* Reset the initial direction if transvserse transport is disabled. */
        if (longitudinal)
                memcpy(
                    state->direction, tau->direction, sizeof(state->direction));

        /* Backward propagate the neutrino. */
        enum ent_event event = ENT_EVENT_NONE;
        while ((event != ENT_EVENT_EXIT) &&
            (state->energy < energy_cut - FLT_EPSILON)) {
                ent_transport(physics, ctx_ent, state, NULL, &event);
                if (longitudinal)
                        memcpy(state->direction, tau->direction,
                            sizeof(state->direction));

                if (event == ENT_EVENT_DECAY_TAU) {
                        /* Backward randomise the tau decay. */
                        double momentum[3] = { state->energy *
                                    state->direction[0],
                                state->energy * state->direction[1],
                                state->energy * state->direction[2] };
                        double weight;
                        int trials;
                        for (trials = 0; trials < 20; trials++) {
                                if (alouette_undecay(state->pid, momentum,
                                        &polarisation_cb, DECAY_BIAS,
                                        &weight) == ALOUETTE_RETURN_SUCCESS)
                                        break;
                        }

                        int pid1;
                        if ((alouette_product(&pid1, momentum) !=
                                ALOUETTE_RETURN_SUCCESS) ||
                            (abs(pid1) != ENT_PID_TAU))
                                return;
                        const double p12 = momentum[0] * momentum[0] +
                            momentum[1] * momentum[1] +
                            momentum[2] * momentum[2];
                        const double E1 = sqrt(p12 + tau_mass * tau_mass);
                        if (E1 >= energy_cut - FLT_EPSILON) return;

                        /* Update the tau state and start again the BMC
                         * transport.
                         */
                        tau->charge = (pid1 > 0) ? -1. : 1.;
                        tau->kinetic = E1 - tau_mass;
                        tau->distance = state->distance;
                        tau->grammage = state->grammage;
                        tau->time = 0.;
                        tau->weight = state->weight * weight * state->energy *
                            state->energy / p12;
                        memcpy(tau->position, state->position,
                            sizeof(tau->position));
                        if (!longitudinal)
                                memcpy(tau->direction, state->direction,
                                    sizeof(tau->direction));
                        tau->decayed = 0;
                        generation++;
                        transport_backward(ctx_ent, tau, eventid, generation,
                            tau_at_decay, tau_at_production, done);
                        return;
                }
        }
        if (event != ENT_EVENT_EXIT) return;
        enum ent_pid pid0 =
            (tau_at_decay->charge < 0.) ? ENT_PID_NU_TAU : ENT_PID_NU_BAR_TAU;
        if (state->pid != pid0) return;
        (*done)++;

        /* This is a valid event. In flux mode let us dump the tau state and
         * then return.
         */
        if (flux_mode) {
                /* Log the primary neutrino state. */
                format_ancester(eventid, state);

                /* Log the end points of the tau state. */
                format_tau(generation, pid, tau_at_production, tau_at_decay);
                return;
        }

        /* In full mode let us perform the tau decay with ALOUETTE/TAUOLA. */
        const double p = sqrt(
            tau_at_decay->kinetic * (tau_at_decay->kinetic + 2. * tau_mass));
        double momentum[3] = { p * tau_at_decay->direction[0],
                p * tau_at_decay->direction[1],
                p * tau_at_decay->direction[2] };
        int trials;
        for (trials = 0; trials < 20; trials++) {
                if (alouette_decay(pid, momentum, tau_at_decay->direction) ==
                    ALOUETTE_RETURN_SUCCESS)
                        break;
        }

        int pid1, nprod = 0;
        while (alouette_product(&pid1, momentum) == ALOUETTE_RETURN_SUCCESS) {
                if (abs(pid1 == 12) || (abs(pid1) == 13) || (abs(pid1) == 14) ||
                    (abs(pid1) == 16))
                        continue;
                if (nprod == 0) {
                        /* Log the primary neutrino state. */
                        format_ancester(eventid, state);

                        /* Log the end points of the tau state. */
                        format_tau(
                            generation, pid, tau_at_production, tau_at_decay);
                }
                format_decay_product(pid1, momentum);
                nprod++;
        }
}

int main(int argc, char * argv[])
{
        /* Clear the error status. */
        errno = 0;

        /* Set the input arguments. */
        int n_events = 0;
        int use_append = 0, do_interaction = 1, theta_interval = 1;
        int z_interval = 1, elevation_interval = 1, mode_forward = 0;
        double cos_theta_min = 0.15, cos_theta_max = 0.25;
        double elevation_min = -10., elevation_max = 10.;
        double z_min = 1E+00, z_max = 1E+05;
        double energy_min = 1E+07, energy_max = 1E+12;
        int energy_spectrum = 1, energy_analog = 0;
        int pem_sea = 1;
        char * pdf_file = DANTON_DIR "/ent/data/pdf/CT14nnlo_0000.dat";

        /* Parse the optional arguments. */
        for (;;) {
                /* Short options. */
                const char * short_options = "c:e:hn:t:o:";

                /* Long options. */
                struct option long_options[] = { /* Configuration options. */
                        { "altitude", required_argument, NULL, 0 },
                        { "altitude-max", required_argument, NULL, 0 },
                        { "altitude-min", required_argument, NULL, 0 },
                        { "cos-theta", required_argument, NULL, 'c' },
                        { "cos-theta-max", required_argument, NULL, 0 },
                        { "cos-theta-min", required_argument, NULL, 0 },
                        { "elevation", required_argument, NULL, 0 },
                        { "elevation-max", required_argument, NULL, 0 },
                        { "elevation-min", required_argument, NULL, 0 },
                        { "energy", required_argument, NULL, 'e' },
                        { "energy-analog", no_argument, &energy_analog, 1 },
                        { "energy-cut", required_argument, NULL, 0 },
                        { "energy-max", required_argument, NULL, 0 },
                        { "energy-min", required_argument, NULL, 0 },
                        { "forward", no_argument, &mode_forward, 1 },
                        { "longitudinal", no_argument, &longitudinal, 1 },
                        { "pem-no-sea", no_argument, &pem_sea, 0 },
                        { "taus", required_argument, NULL, 't' },

                        /* Control options. */
                        { "append", no_argument, &use_append, 1 },
                        { "flux", required_argument, NULL, 0 },
                        { "grammage", no_argument, &do_interaction, 0 },
                        { "help", no_argument, NULL, 'h' },
                        { "output-file", required_argument, NULL, 'o' },
                        { "pdf-file", required_argument, NULL, 0 },

                        { 0, 0, 0, 0 }
                };

                /* Process the next argument. */
                char * endptr = NULL; /* For parsing with str* functions. */
                int option_index = 0;
                int c = getopt_long(
                    argc, argv, short_options, long_options, &option_index);
                if (c == -1)
                        break; /* No more options to parse. */
                else if (c > 0) {
                        if (c == 'c') {
                                theta_interval = 0;
                                cos_theta_min = strtod(optarg, &endptr);
                        } else if (c == 'e') {
                                energy_spectrum = 0;
                                energy_min = strtod(optarg, &endptr);
                        } else if (c == 'h')
                                exit_with_help(EXIT_SUCCESS);
                        else if (c == 'n')
                                n_events = strtol(optarg, &endptr, 0);
                        else if (c == 't')
                                n_taus = strtol(optarg, &endptr, 0);
                        else if (c == 'o')
                                output_file = optarg;
                        else {
                                /*
                                 * getopt should already have reported
                                 * an error.
                                 */
                                exit(EXIT_FAILURE);
                        }
                } else {
                        /* Setters for long options. */
                        struct opt_setter setters[] = {
                                /* Configuration options. */
                                { &z_min, &opt_strtod, &endptr },
                                { &z_max, &opt_strtod, &endptr },
                                { &z_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* cos-theta */
                                { &cos_theta_max, &opt_strtod, &endptr },
                                { &cos_theta_min, &opt_strtod, &endptr },
                                { &elevation_min, &opt_strtod, &endptr },
                                { &elevation_max, &opt_strtod, &endptr },
                                { &elevation_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* energy */
                                { NULL, NULL, NULL }, /* energy-analog */
                                { &energy_cut, &opt_strtod, &endptr },
                                { &energy_max, &opt_strtod, &endptr },
                                { &energy_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* forward */
                                { NULL, NULL, NULL }, /* longitudinal */
                                { NULL, NULL, NULL }, /* pem-no-sea */
                                { NULL, NULL, NULL }, /* taus */

                                /* Control options. */
                                { NULL, NULL, NULL }, /* append */
                                { &flux_altitude, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* grammage */
                                { NULL, NULL, NULL }, /* help */
                                { NULL, NULL, NULL }, /* output-file */
                                { &pdf_file, &opt_copy, NULL },
                        };

                        /* Set the long option. */
                        struct opt_setter * s = setters + option_index;
                        if (s->variable == NULL) continue;
                        s->set(optarg, s->argument, s->variable);

                        /* Check for extra flags to set. */
                        if (strcmp(long_options[option_index].name,
                                "altitude") == 0)
                                z_interval = 0;
                        else if (strcmp(long_options[option_index].name,
                                     "elevation") == 0)
                                elevation_interval = 0;
                        else if (strcmp(long_options[option_index].name,
                                     "flux") == 0)
                                flux_mode = 1;
                }

                /* Check the parsing. */
                if (endptr == optarg) errno = EINVAL;
                if (errno != 0) {
                        perror("danton");
                        exit(EXIT_FAILURE);
                }
        }

        /* Check the consistency of the options. */
        const double cmin = (mode_forward) ? 0. : -1.;
        if (((cos_theta_min < cmin) || (cos_theta_min > 1.)) ||
            ((theta_interval &&
                ((cos_theta_min >= cos_theta_max) || (cos_theta_max > 1.) ||
                    (cos_theta_max < cmin))))) {
                fprintf(stderr, "danton: inconsistent cos(theta) value(s)."
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }
        if ((z_min <= 0.) || (z_interval && (z_min >= z_max))) {
                fprintf(stderr, "danton: inconsistent altitude value(s)."
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }
        if ((elevation_min < -90.) || (elevation_min > 90.) ||
            (elevation_max > 90.) ||
            (elevation_interval && (elevation_min >= elevation_max))) {
                fprintf(stderr, "danton: inconsistent elevation value(s)."
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }
        if (n_events < 1) n_events = do_interaction ? 10000 : 1001;
        if (!do_interaction && !theta_interval) n_events = 1;
        if (!do_interaction && theta_interval && n_events < 2) {
                fprintf(stderr, "danton: numbers of bins must be 2 or more. "
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }
        if (energy_cut <= 0.) energy_cut = (mode_forward) ? 1E+03 : 1E+12;
        if ((energy_cut < 1E+02) || (energy_min < 1E+02)) {
                fprintf(stderr, "danton: energies must be at least 100 GeV. "
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }
        if (energy_spectrum && (energy_max <= energy_min)) {
                fprintf(stderr, "danton: inconsistent energy range. "
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }

        /* Set the primary. */
        int projectile;
        if (do_interaction) {
                /* Parse and check the mandatory arguments. */
                if (argc - optind != 1) {
                        fprintf(stderr, "danton: wrong number of arguments. "
                                        "Call with -h, --help for usage.\n");
                        exit(EXIT_FAILURE);
                }

                if (sscanf(argv[optind++], "%d", &projectile) != 1)
                        exit_with_help(EXIT_FAILURE);

                if (mode_forward) {
                        if ((projectile != ENT_PID_NU_TAU) &&
                            (projectile != ENT_PID_NU_BAR_TAU) &&
                            (projectile != ENT_PID_NU_BAR_E)) {
                                fprintf(stderr, "danton: invalid neutrino PID."
                                                "Call with -h, --help for "
                                                "usage.\n");
                                exit(EXIT_FAILURE);
                        }
                } else {
                        if (abs(projectile) != ENT_PID_TAU) {
                                fprintf(stderr, "danton: invalid daughter PID."
                                                "Call with -h, --help for "
                                                "usage.\n");
                                exit(EXIT_FAILURE);
                        }
                }

        } else {
                /* This is a grammage scan. Let's use some arbitrary
                 * primary. */
                mode_forward = 1;
                energy_spectrum = 0;
                energy_min = 1E+09;
                projectile = ENT_PID_NU_TAU;
        }

        /* Configure the geometry. */
        if (!pem_sea) {
                /* Replace seas with rock. */
                memcpy(media_ent + 9, media_ent + 8, sizeof(*media_ent));
                memcpy(media_pumas + 9, media_pumas + 8, sizeof(*media_pumas));
        }

        /* Configure the output stream. */
        int print_header = !use_append;
        FILE * stream;
        if (output_file == NULL)
                stream = stdout;
        else if (use_append) {
                stream = fopen(output_file, "a");
                if (stream == NULL) {
                        stream = fopen(output_file, "w+");
                        print_header = 1;
                }
        } else
                stream = fopen(output_file, "w+");

        if (stream == NULL) {
                fprintf(stderr, "danton: could not open the output file. "
                                "Call with -h, --help for usage.\n");
                exit(EXIT_FAILURE);
        }

        if (print_header) {
                if (do_interaction)
                        print_header_decay(stream);
                else
                        print_header_grammage(stream);
                if (output_file != NULL) fclose(stream);
        }

        /* Register the error handlers. */
        ent_error_handler_set(&handle_ent);
        pumas_error_handler_set(&handle_pumas);

        /* Create a new neutrino Physics environment. */
        if (do_interaction) ent_physics_create(&physics, pdf_file);

        /* Initialise the PUMAS transport engine. */
        load_pumas();

        /* Initialise ALOUETTE/TAUOLA. */
        enum alouette_return a_rc;
        if ((a_rc = alouette_initialise(1, NULL)) != ALOUETTE_RETURN_SUCCESS) {
                fprintf(stderr, "alouette_initialise: %s\n",
                    alouette_strerror(a_rc));
                gracefully_exit(EXIT_FAILURE);
        };

        /* Initialise the random engine. */
        random_initialise();

        /* Initialise the Monte-Carlo contexts. */
        struct ent_context ctx_ent = { &medium_ent,
                (ent_random_cb *)&random_uniform01, NULL };
        pumas_context_create(0, &ctx_pumas);
        ctx_pumas->medium = &medium_pumas;
        ctx_pumas->random = (pumas_random_cb *)&random_uniform01;
        ctx_pumas->kinetic_limit = energy_cut - tau_mass;
        ctx_pumas->longitudinal = longitudinal;

        int done = 0;
        if (mode_forward) {
                /* Run a bunch of forward Monte-Carlo events. */
                long i;
                for (i = 0; i < n_events; i++) {
                        double ct;
                        if (theta_interval) {
                                const double u = do_interaction ?
                                    random_uniform01(NULL) :
                                    i / (n_events - 1.);
                                ct = (cos_theta_max - cos_theta_min) * u +
                                    cos_theta_min;
                        } else
                                ct = cos_theta_min;
                        const double st = sqrt(1. - ct * ct);
                        double energy, weight = 1.;
                        if (energy_spectrum) {
                                if (energy_analog) {
                                        const double ei0 = 1. / energy_min;
                                        energy =
                                            1. /
                                            (ei0 +
                                                random_uniform01(NULL) *
                                                    (ei0 - 1. / energy_max));
                                } else {
                                        const double r =
                                            log(energy_max / energy_min);
                                        energy = energy_min *
                                            exp(r * random_uniform01(NULL));
                                        weight = r * energy_max * energy_min /
                                            ((energy_max - energy_min) *
                                                     energy);
                                }
                        } else
                                energy = energy_min;
                        struct generic_state state = {
                                .base.ent = { projectile, energy, 0., 0.,
                                    weight, { 0., 0., -EARTH_RADIUS - 1E+05 },
                                    { st, 0., ct } },
                                .r = 0.,
                                .is_tau = 0,
                                .is_inside = -1,
                                .has_crossed = 0
                        };
                        struct ent_state ancester;
                        memcpy(&ancester, &state.base.ent, sizeof(ancester));
                        transport_forward(&ctx_ent, (struct ent_state *)&state,
                            i, 1, 0, &ancester, &done);
                        if ((n_taus > 0) && (done >= n_taus)) break;
                        if (!do_interaction)
                                format_grammage(ct, state.base.ent.grammage);
                }
        } else {
                /* Run a bunch of backward Monte-Carlo events. */
                ctx_ent.ancester = &ancester_cb;
                ctx_pumas->forward = 0;

                cos_theta_min = cos((90. - elevation_min) * M_PI / 180.);
                if (elevation_interval)
                        cos_theta_max =
                            cos((90. - elevation_max) * M_PI / 180.);

                long i;
                for (i = 0; i < n_events; i++) {
                        double ct, weight = 1.;
                        if (elevation_interval) {
                                const double u = do_interaction ?
                                    random_uniform01(NULL) :
                                    i / (n_events - 1.);
                                const double dc = cos_theta_max - cos_theta_min;
                                ct = dc * u + cos_theta_min;
                                weight *= dc;
                        } else
                                ct = cos_theta_min;
                        const double st = sqrt(1. - ct * ct);
                        double energy;
                        if (energy_spectrum) {
                                const double r = log(energy_max / energy_min);
                                energy = energy_min *
                                    exp(r * random_uniform01(NULL));
                                weight *= r * energy;
                        } else
                                energy = energy_min;
                        double z0;
                        if (flux_mode)
                                z0 = flux_altitude;
                        else if (z_interval) {
                                const double r = log(z_max / z_min);
                                z0 = z_min * exp(r * random_uniform01(NULL));
                                weight *= r * z0;
                        } else
                                z0 = z_min;
                        const double charge = (projectile > 0) ? -1. : 1.;
                        struct generic_state state = {
                                .base.pumas = { charge, energy - tau_mass, 0.,
                                    0., 0., weight,
                                    { 0., 0., EARTH_RADIUS + z0 },
                                    { st, 0., ct }, 0 },
                                .r = 0.,
                                .is_tau = 1,
                                .is_inside = -1,
                                .has_crossed = -1
                        };
                        struct pumas_state tau_at_decay, tau_at_production;
                        transport_backward(&ctx_ent,
                            (struct pumas_state *)&state, i, 1, &tau_at_decay,
                            &tau_at_production, &done);
                        if ((n_taus > 0) && (done > n_taus)) break;
                }
        }

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_SUCCESS);
}
