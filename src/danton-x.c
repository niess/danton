/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C99 executable simulating the coupled transport of ultra
 * high energy taus and neutrinos through the Earth, by Monte-Carlo.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>
 */

/* Standard library includes. */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The DANTON API. */
#include "danton.h"
#include "danton/primary/discrete.h"
#include "danton/primary/powerlaw.h"
#include "danton/recorder/text.h"

/* Jasmine with some tea, for parsing the data card in JSON format. */
#include "jsmn-tea.h"

/* Handle for the simulation context. */
static struct danton_context * context = NULL;

/* Handle for parsing JSON card(s). */
static struct jsmn_tea * tea = NULL;

/* Handle for handling errors. */
static struct roar_handler handler = { NULL, NULL, NULL, NULL };

/* Path of the current card. */
static const char * card_path = NULL;

/* Global options for the stepping. */
static struct {
        char * path;
        int append;
        int verbosity;
} stepping_options = { NULL, 0, 0 };

/* Finalise and exit to the OS. */
static int gracefully_exit(int rc)
{
        /* Finalise and exit to the OS. */
        jsmn_tea_destroy(&tea);
        int i;
        for (i = 0; i < DANTON_PARTICLE_N_NU; i++)
                danton_destroy((void **)&context->primary[i]);
        danton_destroy((void **)&context->recorder);
        danton_destroy((void **)&context->sampler);
        danton_context_destroy(&context);
        danton_finalise();
        free(stepping_options.path);
        exit(rc);
}

/* Post-error callback. */
int handle_post_error(
    struct roar_handler * handler, roar_function_t * referent, int code)
{
        /* Finalise and exit to the OS. */
        return gracefully_exit(EXIT_FAILURE);
}

/* Callback for catching errors. */
int catch_error(
    struct roar_handler * handler, roar_function_t * referent, int code)
{
        return EXIT_SUCCESS;
}

/* Show help and exit. */
static void exit_with_help(int code)
{
        // clang-format off
        fprintf(stderr,
"Usage: danton [DATACARD.JSON]\n"
"Simulate the coupled transport of ultra high energy taus and neutrinos\n"
"through the Earth, by Monte-Carlo.\n"
"\n"
"Data card:\n"
"Syntax and examples available from https://github.com/niess/danton.\n"
"\n"
"Exit status:\n"
" %d  if OK,\n"
" %d  if an error occurred.\n"
"\n"
"License: GNU LGPLv3\n"
"Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC.\n"
"Author: Valentin NIESS (niess@in2p3.fr)\n"
"\n", EXIT_SUCCESS, EXIT_FAILURE);
        exit(code);
        // clang-format on
}

/* Get the next token in the JSON card as a double or as an array of
 * two doubles.
 */
static void card_get_range(const char * field, double * range)
{
        handler.pre = &catch_error;
        int size;
        int rc = jsmn_tea_next_array(tea, &size);
        handler.pre = NULL;
        if (rc < 0) {
                jsmn_tea_next_number(tea, JSMN_TEA_TYPE_DOUBLE, range);
                range[1] = range[0];
        } else if (size != 2) {
                ROAR_ERRNO_FORMAT(&handler, &card_get_range, EINVAL,
                    "[%s #%d] invalid array size for field `%s`", card_path,
                    tea->index, field);
        } else {
                jsmn_tea_next_number(tea, JSMN_TEA_TYPE_DOUBLE, range);
                jsmn_tea_next_number(tea, JSMN_TEA_TYPE_DOUBLE, range + 1);
        }
}

/* Update the event recorder according to the data card. */
static void card_update_recorder(void)
{
        char * output_file;
        jsmn_tea_next_string(tea, 0, &output_file);
        danton_destroy((void **)&context->recorder);
        context->recorder =
            (struct danton_recorder *)danton_text_create(output_file);
        if (context->recorder == NULL) gracefully_exit(EXIT_FAILURE);
}

/* Update DANTON's run mode according to the data card. */
static void card_update_mode(void)
{
        char * s;
        jsmn_tea_next_string(tea, 0, &s);
        if (strcmp(s, "backward") == 0)
                context->mode = DANTON_MODE_BACKWARD;
        else if (strcmp(s, "forward") == 0)
                context->mode = DANTON_MODE_FORWARD;
        else if (strcmp(s, "grammage") == 0)
                context->mode = DANTON_MODE_GRAMMAGE;
        else {
                ROAR_ERRNO_FORMAT(&handler, &card_update_mode, EINVAL,
                    "[%s #%d] invalid mode `%s`", card_path, tea->index, s);
        }
}

/* Update DANTON's PRNG seed according to the data card. */
static void card_update_seed(int * seeded, unsigned long * seed)
{
        handler.pre = &catch_error;
        int rc = jsmn_tea_next_null(tea);
        handler.pre = NULL;
        if (rc < 0) {
                jsmn_tea_next_number(tea, JSMN_TEA_TYPE_UNSIGNED_LONG, seed);
                *seeded = 1;
        }
}

/* List of particle names, following DANTON's ordering. */
static const char * particle_name[DANTON_PARTICLE_N] = { "nu_tau~", "nu_mu~",
        "nu_e~", "nu_e", "nu_mu", "nu_tau", "tau~", "tau" };

/* Get a particle name from the card and convert it to a DANTON index. */
static enum danton_particle card_get_particle(enum danton_particle index_max)
{
        char * particle;
        jsmn_tea_next_string(tea, 1, &particle);
        enum danton_particle index;
        for (index = 0; index < index_max; index++) {
                if (strcmp(particle, particle_name[index]) == 0) return index;
        }
        return ROAR_ERRNO_FORMAT(&handler, &card_get_particle, EINVAL,
            "[%s #%d] invalid particle name `%s`", card_path, tea->index,
            particle);
}

/* Update the sampler according to the data card. */
static void card_update_sampler(void)
{
        /* Create and initialise the sampler if required. */
        struct danton_sampler * sampler;
        if (context->sampler == NULL) {
                sampler = danton_sampler_create();
                if (sampler == NULL) {
                        ROAR_ERRWP_MESSAGE(&handler, &card_update_sampler, -1,
                            "danton error", danton_error_pop(NULL));
                }
                sampler->latitude = 45.;
                sampler->longitude = 0.;
                sampler->altitude[0] = 1E-03;
                sampler->altitude[1] = 1E+04;
                sampler->azimuth[0] = -180.;
                sampler->azimuth[1] = 180.;
                sampler->elevation[0] = -10.;
                sampler->elevation[1] = 10.;
                sampler->energy[0] = 1E+06;
                sampler->energy[1] = 1E+12;
                context->sampler = sampler;
        } else
                sampler = context->sampler;

        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * tag;
                jsmn_tea_next_string(tea, 1, &tag);
                if (strcmp(tag, "latitude") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &sampler->latitude);
                else if (strcmp(tag, "longitude") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &sampler->longitude);
                else if (strcmp(tag, "altitude") == 0)
                        card_get_range(tag, sampler->altitude);
                else if (strcmp(tag, "azimuth") == 0)
                        card_get_range(tag, sampler->azimuth);
                else if (strcmp(tag, "elevation") == 0)
                        card_get_range(tag, sampler->elevation);
                else if (strcmp(tag, "energy") == 0)
                        card_get_range(tag, sampler->energy);
                else if (strcmp(tag, "weight") == 0) {
                        int j;
                        for (jsmn_tea_next_object(tea, &j); j; j--) {
                                const int k =
                                    card_get_particle(DANTON_PARTICLE_N);
                                jsmn_tea_next_number(tea, JSMN_TEA_TYPE_DOUBLE,
                                    sampler->weight + k);
                        }
                } else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_sampler,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, tag);
                }
        }
}

/* Update the primary flux with a power law model. */
static struct danton_primary * card_update_primary_powerlaw(void)
{
        /* Parse the model parameters. */
        double energy[2] = { 1E+06, 1E+12 };
        double weight = 1.;
        double exponent = -2.;
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * field;
                jsmn_tea_next_string(tea, 1, &field);
                if (strcmp(field, "energy") == 0)
                        card_get_range(field, energy);
                else if (strcmp(field, "exponent") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &exponent);
                else if (strcmp(field, "weight") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &weight);
                else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_sampler,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, field);
                }
        }

        /* Create the model. */
        struct danton_primary * primary =
            (struct danton_primary *)danton_powerlaw_create(
                energy[0], energy[1], exponent, weight);
        if (primary == NULL) {
                ROAR_ERRWP_MESSAGE(&handler, &card_update_primary_powerlaw, -1,
                    "danton error", danton_error_pop(NULL));
        }
        return primary;
}

/* Update the primary flux with a discrete model. */
static struct danton_primary * card_update_primary_discrete(void)
{
        /* Parse the model parameters. */
        double energy = 1E+12;
        double weight = 1.;
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * field;
                jsmn_tea_next_string(tea, 1, &field);
                if (strcmp(field, "energy") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &energy);
                else if (strcmp(field, "weight") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &weight);
                else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_sampler,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, field);
                }
        }

        /* Create the model. */
        struct danton_primary * primary =
            (struct danton_primary *)danton_discrete_create(energy, weight);
        if (primary == NULL) {
                ROAR_ERRWP_MESSAGE(&handler, &card_update_primary_discrete, -1,
                    "danton error", danton_error_pop(NULL));
        }
        return primary;
}

/* Update the primary flux according to the data card. */
static void card_update_primary(void)
{
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                /* Parse the particle index. */
                const int index = card_get_particle(DANTON_PARTICLE_N_NU);

                /* Check that one has a length 2 array. */
                int size;
                jsmn_tea_next_array(tea, &size);
                if (size != 2) {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_primary,
                            EINVAL,
                            "[%s #%d] invalid array size for field `%s`",
                            card_path, tea->index, particle_name[index]);
                }

                /* Clear any existing model. */
                danton_destroy((void **)&context->primary[index]);

                /* Parse the primary model.  */
                char * model;
                jsmn_tea_next_string(tea, 0, &model);
                if (strcmp(model, "power-law") == 0)
                        context->primary[index] =
                            card_update_primary_powerlaw();
                else if (strcmp(model, "discrete") == 0)
                        context->primary[index] =
                            card_update_primary_discrete();
                else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_primary,
                            EINVAL,
                            "[%s #%d] invalid primary model `%s` for particle "
                            "`%s`",
                            card_path, tea->index, model, particle_name[index]);
                }
        }
}

/* Update the Earth model according to the data card. */
static void card_update_earth_model(void)
{
        /* Set the default parameter values. */
        char * reference = NULL;
        char * topography = NULL;
        char * material = NULL;
        double zf = 0., density = 0.;
        int sea = -1, flat_topography = 0;

        /* Parse the data card. */
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * field;
                jsmn_tea_next_string(tea, 1, &field);
                if (strcmp(field, "reference") == 0) {
                        jsmn_tea_next_string(tea, 0, &reference);
                } else if (strcmp(field, "topography") == 0) {
                        handler.pre = &catch_error;
                        int rc = jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &zf);
                        handler.pre = NULL;
                        if (rc < 0) {
                                jsmn_tea_next_string(tea, 0, &topography);
                        } else {
                                flat_topography = 1;
                        }
                } else if (strcmp(field, "material") == 0) {
                        jsmn_tea_next_string(tea, 0, &material);
                } else if (strcmp(field, "density") == 0) {
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_DOUBLE, &density);
                } else if (strcmp(field, "sea") == 0) {
                        jsmn_tea_next_bool(tea, &sea);
                } else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_earth_model,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, field);
                }
        }

        /* Check the topography */
        char buffer[128];
        if (flat_topography) {
                snprintf(buffer, 127, "flat://%lf", zf);
                topography = buffer;
        }

        /* Configure the Earth model. */
        int * ptr_sea = (sea >= 0) ? &sea : NULL;
        if (danton_earth_model(reference, topography, material,
                density, ptr_sea) != EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(&handler, &card_update_earth_model, -1,
                    "danton error", danton_error_pop(NULL));
        }
}

/* Update the stepping options according to the data card. */
static void card_update_stepping(void)
{
        /* Parse the data card. */
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * field;
                jsmn_tea_next_string(tea, 1, &field);
                if (strcmp(field, "append") == 0) {
                        jsmn_tea_next_bool(tea, &stepping_options.append);
                } else if (strcmp(field, "path") == 0) {
                        free(stepping_options.path);
                        char * s;
                        jsmn_tea_next_string(tea, 0, &s);
                        if (s != NULL) {
                                const int n = strlen(s) + 1;
                                stepping_options.path = malloc(n);
                                memcpy(stepping_options.path, s, n);
                        }
                } else if (strcmp(field, "verbosity") == 0) {
                        jsmn_tea_next_number(tea, JSMN_TEA_TYPE_INT,
                            &stepping_options.verbosity);
                } else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_stepping,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, field);
                }
        }
}

/* Update DANTON's configuration according to the content of the data card. */
static void card_update(
    int * n_events, int * n_requested, int * seeded, unsigned long * seed)
{
        /* Loop over the fields. */
        int i;
        for (jsmn_tea_next_object(tea, &i); i; i--) {
                char * tag;
                jsmn_tea_next_string(tea, 1, &tag);
                if (strcmp(tag, "events") == 0)
                        jsmn_tea_next_number(tea, JSMN_TEA_TYPE_INT, n_events);
                else if (strcmp(tag, "requested") == 0)
                        jsmn_tea_next_number(
                            tea, JSMN_TEA_TYPE_INT, n_requested);
                else if (strcmp(tag, "output-file") == 0)
                        card_update_recorder();
                else if (strcmp(tag, "mode") == 0)
                        card_update_mode();
                else if (strcmp(tag, "seed") == 0)
                        card_update_seed(seeded, seed);
                else if (strcmp(tag, "decay") == 0)
                        jsmn_tea_next_bool(tea, &context->decay);
                else if (strcmp(tag, "longitudinal") == 0)
                        jsmn_tea_next_bool(tea, &context->longitudinal);
                else if (strcmp(tag, "particle-sampler") == 0)
                        card_update_sampler();
                else if (strcmp(tag, "primary-flux") == 0)
                        card_update_primary();
                else if (strcmp(tag, "earth-model") == 0)
                        card_update_earth_model();
                else if (strcmp(tag, "stepping") == 0)
                        card_update_stepping();
                else {
                        ROAR_ERRNO_FORMAT(&handler, &card_update_sampler,
                            EINVAL, "[%s #%d] invalid key `%s`", card_path,
                            tea->index, tag);
                }
        }
}

/* Dump the given step to the opened file. */
static void dump_step(
    FILE * fid, int index, int medium, struct danton_state * state)
{
        if (index > 0) fputs(", ", fid);
        fprintf(fid, "[%3d, %3d, %.5E, %.3f, %.3f, %.3f]", medium, state->pid,
            state->energy, state->position[0], state->position[1],
            state->position[2]);
}

/* Custom run action for dumping an event's steps. */
static int dump_steps(struct danton_context * context,
    enum danton_run_event event, int medium, struct danton_state * state)
{
        static FILE * fid = NULL;
        static int index = 0;

        if (event == DANTON_RUN_EVENT_START) {
                fid = fopen(stepping_options.path, "a");
                if (fid == NULL) {
                        ROAR_ERRNO_MESSAGE(
                            &handler, &dump_steps, 0, stepping_options.path);
                        return EXIT_FAILURE;
                }
                index = 0;
                fprintf(fid, "[");
                dump_step(fid, index++, medium, state);
        } else if (event == DANTON_RUN_EVENT_STEP) {
                dump_step(fid, index++, medium, state);
        } else {
                fprintf(fid, "]\n");
                fclose(fid);
                fid = NULL;
        }

        return EXIT_SUCCESS;
}

int main(int argc, char * argv[])
{
        /* Configure the error handler. */
        errno = 0;
        handler.stream = stderr;
        handler.post = &handle_post_error;

        /* If no data card was provided let us show the help message and
         * exit.
         */
        if (argc <= 1) exit_with_help(EXIT_SUCCESS);

        /* Initialise DANTON. */
        if (danton_initialise(NULL, NULL, NULL, NULL, NULL) != EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(&handler, &main, -1, "danton error",
                    danton_error_pop(NULL));
        }

        /* Create a simulation context. */
        context = danton_context_create();
        if (context == NULL) {
                ROAR_ERRWP_MESSAGE(&handler, &main, -1, "danton error",
                    danton_error_pop(NULL));
        }

        /* Set the input arguments from the JSON card(s). */
        int n_events = 10000, n_requested = 0, seeded = 0;
        unsigned long seed;
        for (argv++; *argv != NULL; argv++) {
                tea = jsmn_tea_create(*argv, JSMN_TEA_MODE_LOAD, &handler);
                card_path = *argv;
                card_update(&n_events, &n_requested, &seeded, &seed);
                jsmn_tea_destroy(&tea);
        }

        /* Update the PRNG, if needed. */
        if (seeded) {
                danton_context_random_set(context, &seed);
        }

        /* Initialise any stepping dump. */
        if (stepping_options.path != NULL) {
                context->run_action = &dump_steps;
                if (!stepping_options.append) {
                        FILE * fid = fopen(stepping_options.path, "w+");
                        fclose(fid);
                }
        } else {
                context->run_action = NULL;
        }

        /* Update the particle sampler. */
        if (danton_sampler_update(context->sampler) != EXIT_SUCCESS) {
                ROAR_ERRWP_MESSAGE(&handler, &main, -1, "danton error",
                    danton_error_pop(NULL));
        }

        /* Run the simulation. */
        if (danton_context_run(context, n_events, n_requested) != EXIT_SUCCESS)
                ROAR_ERRWP_MESSAGE(&handler, &main, -1, "danton error",
                    danton_error_pop(context));

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_SUCCESS);
}
