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

/* The DANTON API. */
#include "danton.h"

/* Finalise and exit to the OS. */
static void gracefully_exit(int rc)
{
        /* Finalise and exit to the OS. */
        danton_finalise();
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
"      --altitude-min=Z       set the minimum decay altitude to Z [1E-03]\n"
"  -c, --cos-theta=C          switch to a point estimate with cos(theta)=C\n"
"      --cos-theta-max=C      set the maximum value of cos(theta) to C [0.25]\n"
"      --cos-theta-min=C      set the minimum value of cos(theta) to C [0.15]\n"
"      --decay-disable        disable the final tau decays and switch to a\n"
"                               flux sampling\n"
"      --elevation=A          switch to a point estimate at an elevation\n"
"                               angle of A\n"
"      --elevation-max=A      set the maximum elevation angle to A [10.]\n"
"      --elevation-min=A      set the minimum elevation angle to A [-10.]\n"
"  -e, --energy=E             switch to a monokinetic beam of primaries with\n"
"                               energy E\n"
"      --energy-max=E         set to E the maximum primary energy [1E+12]\n"
"      --energy-min=E         set to E the minimum primary energy [1E+07]\n"
"      --forward              switch to forward Monte-Carlo\n"
"      --longitudinal         disable the transverse transport\n"
"  -n                         set the number of incoming neutrinos [10000] or\n"
"                               the number of bins for the grammage sampling\n"
"                               [1001]\n"
"      --pem-no-sea           Replace seas with rocks in the Preliminary\n"
"                               Earth Model (PEM)\n"
"\n"
"Control options:\n"
"      --grammage             disable neutrino interactions but compute the\n"
"                               total grammage along the path.\n"
"  -h, --help                 show this help and exit\n"
"  -o, --output-file=FILE     set the output FILE for the computed flux\n"
"      --pdf-file=FILE        specify a pdf FILE for partons distributions\n"
"                               [(builtin)/CT14nnlo_0000.dat]. The format\n"
"                               must be `lhagrid1` compliant.\n"
"\n"
"Energies (E) must be given in GeV, altitudes (Z) in m and angles (A) in deg.\n"
"PID must be one of 15 (tau) or -15 (tau_bar) in backward mode and -12 \n"
"(nu_e_bar), 16 (nu_tau) or -16 (nu_tau_bar) in forward mode.\n"
"\n"
"The default behaviour is to randomise the primary neutrino energy over a\n"
"1/E^2 spectrum with a log bias. The primary direction is randomised uniformly\n"
"over a solid angle specified by the min and max value of the cosine of the\n"
"angle theta with the local vertical at the atmosphere entrance. Use the -c,\n"
"--cos-theta option in order to perform a point estimate instead.\n"
"\n");
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

int main(int argc, char * argv[])
{
        /* Set the input arguments. */
        int n_events = 1000;
        int longitudinal = 0, grammage = 0, decay = 1, forward = 0;
        int cos_theta_range = 1, elevation_range = 1, altitude_range = 1,
            energy_range = 1;
        double cos_theta_min = 0.15, cos_theta_max = 0.25;
        double elevation_min = -10., elevation_max = 10.;
        double altitude_min = 1E-03, altitude_max = 1E+05;
        double energy_min = 1E+07, energy_max = 1E+12;
        int pem_sea = 1;
        char * pdf_file = NULL;
        char * output_file = NULL;

        /* Parse the optional arguments. */
        for (;;) {
                /* Short options. */
                const char * short_options = "c:e:hn:o:";

                /* Long options. */
                struct option long_options[] = { /* Configuration options. */
                        { "altitude", required_argument, NULL, 0 },
                        { "altitude-max", required_argument, NULL, 0 },
                        { "altitude-min", required_argument, NULL, 0 },
                        { "cos-theta", required_argument, NULL, 'c' },
                        { "cos-theta-max", required_argument, NULL, 0 },
                        { "cos-theta-min", required_argument, NULL, 0 },
                        { "decay-disable", no_argument, &decay, 0 },
                        { "elevation", required_argument, NULL, 0 },
                        { "elevation-max", required_argument, NULL, 0 },
                        { "elevation-min", required_argument, NULL, 0 },
                        { "energy", required_argument, NULL, 'e' },
                        { "energy-max", required_argument, NULL, 0 },
                        { "energy-min", required_argument, NULL, 0 },
                        { "forward", no_argument, &forward, 1 },
                        { "longitudinal", no_argument, &longitudinal, 1 },
                        { "pem-no-sea", no_argument, &pem_sea, 0 },

                        /* Control options. */
                        { "grammage", no_argument, &grammage, 1 },
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
                                cos_theta_range = 0;
                                cos_theta_min = strtod(optarg, &endptr);
                        } else if (c == 'e') {
                                energy_range = 0;
                                energy_min = strtod(optarg, &endptr);
                        } else if (c == 'h')
                                exit_with_help(EXIT_SUCCESS);
                        else if (c == 'n')
                                n_events = strtol(optarg, &endptr, 0);
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
                                { &altitude_min, &opt_strtod, &endptr },
                                { &altitude_max, &opt_strtod, &endptr },
                                { &altitude_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* cos-theta */
                                { &cos_theta_max, &opt_strtod, &endptr },
                                { &cos_theta_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* decay-disable */
                                { &elevation_min, &opt_strtod, &endptr },
                                { &elevation_max, &opt_strtod, &endptr },
                                { &elevation_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* energy */
                                { &energy_max, &opt_strtod, &endptr },
                                { &energy_min, &opt_strtod, &endptr },
                                { NULL, NULL, NULL }, /* forward */
                                { NULL, NULL, NULL }, /* longitudinal */
                                { NULL, NULL, NULL }, /* pem-no-sea */

                                /* Control options. */
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
                                altitude_range = 0;
                        else if (strcmp(long_options[option_index].name,
                                     "elevation") == 0)
                                elevation_range = 0;
                }

                /* Check the parsing. */
                if (endptr == optarg) errno = EINVAL;
                if (errno != 0) {
                        perror("danton");
                        exit(EXIT_FAILURE);
                }
        }

        /* Set the primary. */
        int target;
        if (!grammage) {
                /* Parse and check the mandatory arguments. */
                if (argc - optind != 1) {
                        fprintf(stderr, "danton: wrong number of arguments. "
                                        "Call with -h, --help for usage.\n");
                        exit(EXIT_FAILURE);
                }
                if (sscanf(argv[optind++], "%d", &target) != 1)
                        exit_with_help(EXIT_FAILURE);
        } else {
                /* Set the target to a nu tau for a grammage scan. */
                target = 16;
        }

        /* Initialise DANTON. */
        if (danton_initialise(pdf_file, NULL, NULL) != EXIT_SUCCESS) {
                fprintf(stderr, "danton: couldn't initialise the library.\n");
                exit(EXIT_FAILURE);
        }
        if (!pem_sea) danton_pem_dry();

        /* Create a sampler. */
        struct danton_sampler * sampler = danton_sampler_create();
        if (sampler == NULL) {
                fprintf(stderr, "danton: could not create the sampler.\n");
                exit(EXIT_FAILURE);
        }

        /* Configure the sampler. */
        sampler->altitude[0] = altitude_min;
        sampler->altitude[1] = altitude_range ? altitude_max : altitude_min;
        if (forward) {
                sampler->cos_theta[0] = cos_theta_min;
                sampler->cos_theta[1] =
                    cos_theta_range ? cos_theta_max : cos_theta_min;
        } else {
                sampler->elevation[0] = elevation_min;
                sampler->elevation[1] =
                    elevation_range ? elevation_max : elevation_min;
        }
        sampler->energy[0] = energy_min;
        sampler->energy[1] = energy_range ? energy_max : energy_min;
        enum danton_particle index = danton_particle_index(target);
        sampler->weight[index] = 1.;
        if (danton_sampler_update(sampler) != EXIT_SUCCESS) exit(EXIT_FAILURE);

        /* Create a new simulation context. */
        struct danton_context * context = danton_context_create();
        if (context == NULL) {
                fprintf(stderr,
                    "danton: could not create the simulation context.\n");
                exit(EXIT_FAILURE);
        }

        /* Configure the simulation context. */
        context->forward = forward;
        context->longitudinal = longitudinal;
        context->decay = decay;
        context->grammage = grammage;
        context->sampler = sampler;
        context->output = output_file;

        /* Run the simulation. */
        danton_run(context, n_events);

        /* Finalise and exit to the OS. */
        danton_context_destroy(&context);
        danton_sampler_destroy(&sampler);
        gracefully_exit(EXIT_SUCCESS);
}
