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

/* Jasmine, for parsing the data card in JSON format. */
#include "jsmn.h"

/* Handle for the simulation context. */
static struct danton_context * context = NULL;

/* Container for the parsing of the data card. */
static struct {
        jsmn_parser parser;
        int cursor;
        int n_tokens;
        char * path;
        jsmntok_t * tokens;
        char * buffer;
} card = { { 0, 0, 0 }, 0, 0, NULL, NULL, NULL };

/* Finalise and exit to the OS. */
static void gracefully_exit(int rc)
{
        /* Finalise and exit to the OS. */
        free(card.tokens);
        free(card.buffer);
        int i;
        for (i = 0; i < DANTON_PARTICLE_N_NU; i++)
                danton_destroy((void **)&context->primary[i]);
        danton_destroy((void **)&context->recorder);
        danton_destroy((void **)&context->sampler);
        danton_context_destroy(&context);
        danton_finalise();
        exit(rc);
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

/* Some recurrent error messages. */
static const char * string_json_error =
    "%s (%d): invalid JSON file `%s` (%d).\n";
static const char * string_key_error =
    "%s (%d): invalid key `%s` in JSON file `%s`.\n";
static const char * string_memory_error =
    "%s (%d): could not allocate memory.\n";
static const char * string_particle_error =
    "%s (%d): invalid particle name `%s` in file `%s`\n";
static const char * string_size_error =
    "%s (%d): invalid array size (%d != %d) for field `%s` in file `s`.\n";
static const char * string_type_error =
    "%s (%d): unexpected type (%d : %s) in JSON file `%s`.\n";
static const char * string_value_error = "%s (%d): value error in file `%s` "
                                         "while parsing field `%s`. Expected "
                                         "%s (%s).\n";

/* Helper macro for a memory error. */
#define MEMORY_ERROR                                                           \
        {                                                                      \
                fprintf(stderr, string_memory_error, __FILE__, __LINE__);      \
                exit(EXIT_FAILURE);                                            \
        }

/* Load a data card from a file and parse the JSON tokens. */
static void card_load(char * path)
{
        size_t card_size;
        FILE * stream = fopen(path, "rb");
        if (stream == NULL) {
                fprintf(stderr, "%s (%d): could not open file `%s`.\n",
                    __FILE__, __LINE__, path);
                gracefully_exit(EXIT_FAILURE);
        }
        fseek(stream, 0, SEEK_END);
        card_size = ftell(stream);
        rewind(stream);
        card.buffer = malloc((card_size + 1) * sizeof(*card.buffer));
        if (card.buffer == NULL) MEMORY_ERROR
        const int nread =
            fread(card.buffer, sizeof(*card.buffer), card_size, stream);
        fclose(stream);
        if (nread != card_size) {
                fprintf(stderr, "%s (%d): error while reading file `%s`.\n",
                    __FILE__, __LINE__, path);
                gracefully_exit(EXIT_FAILURE);
        }
        card.buffer[card_size] = 0;
        card.path = path;

        /* Parse the JSON data card. 1st let us check the number of tokens.*/
        jsmn_init(&card.parser);
        card.n_tokens =
            jsmn_parse(&card.parser, card.buffer, card_size, card.tokens, 0);
        if (card.n_tokens <= 0) {
                fprintf(stderr, string_json_error, __FILE__, __LINE__,
                    card.path, card.n_tokens);
                gracefully_exit(EXIT_FAILURE);
        }

        /* Then let us allocate the tokens array and parse the file again. */
        card.tokens = malloc(card.n_tokens * sizeof(*card.tokens));
        if (card.tokens == NULL) MEMORY_ERROR
        jsmn_init(&card.parser);
        const int jrc = jsmn_parse(
            &card.parser, card.buffer, card_size, card.tokens, card.n_tokens);
        if (jrc < 0) {
                fprintf(stderr, string_json_error, __FILE__, __LINE__,
                    card.path, jrc);
                gracefully_exit(EXIT_FAILURE);
        }
}

/* Get the next token in the JSON card as a dictionary header. */
static int json_get_object(void)
{
        jsmntok_t * token = card.tokens + card.cursor++;
        if (token->type != JSMN_OBJECT) {
                card.buffer[token->end] = 0;
                const char * s = card.buffer + token->start;
                fprintf(stderr, string_type_error, __FILE__, __LINE__,
                    token->type, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
        return token->size;
}

/* Get the next token in the JSON card as a string. */
static const char * json_get_string(void)
{
        jsmntok_t * token = card.tokens + card.cursor++;
        card.buffer[token->end] = 0;
        const char * s = card.buffer + token->start;
        if ((token->type == JSMN_PRIMITIVE) && (strcmp(s, "null") == 0))
                return NULL;
        if (token->type != JSMN_STRING) {
                fprintf(stderr, string_type_error, __FILE__, __LINE__,
                    token->type, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
        return s;
}

/* Get the next token in the JSON card as an integer. */
static int json_get_int(const char * field)
{
        jsmntok_t * token = card.tokens + card.cursor++;
        card.buffer[token->end] = 0;
        const char * s = card.buffer + token->start;
        if (token->type != JSMN_PRIMITIVE) {
                fprintf(stderr, string_type_error, __FILE__, __LINE__,
                    token->type, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
        char * endptr;
        long l = strtol(s, &endptr, 0);
        if (*endptr != 0) {
                fprintf(stderr, string_value_error, __FILE__, __LINE__,
                    card.path, field, "an integer", s);
                gracefully_exit(EXIT_FAILURE);
        }
        return (int)l;
}

/* Get the next token in the JSON card as a boolean. */
static int json_get_bool(const char * field)
{
        jsmntok_t * token = card.tokens + card.cursor++;
        card.buffer[token->end] = 0;
        const char * s = card.buffer + token->start;
        if (token->type != JSMN_PRIMITIVE) {
                fprintf(stderr, string_type_error, __FILE__, __LINE__,
                    token->type, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
        if (strcmp(s, "true") == 0)
                return 1;
        else if (strcmp(s, "false") == 0)
                return 0;
        else {
                fprintf(stderr, string_value_error, __FILE__, __LINE__,
                    card.path, field, "a boolean", s);
                gracefully_exit(EXIT_FAILURE);
        }
        return -1; /* For compiler warning ... */
}

/* Get the next token in the JSON card as a double. */
static double json_get_double(const char * field)
{
        jsmntok_t * token = card.tokens + card.cursor++;
        card.buffer[token->end] = 0;
        const char * s = card.buffer + token->start;
        if (token->type != JSMN_PRIMITIVE) {
                fprintf(stderr, string_type_error, __FILE__, __LINE__,
                    token->type, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
        char * endptr;
        double d = strtod(s, &endptr);
        if (*endptr != 0) {
                fprintf(stderr, string_value_error, __FILE__, __LINE__,
                    card.path, field, "a floating number", s);
                gracefully_exit(EXIT_FAILURE);
        }
        return d;
}

/* Get the next token in the JSON card as a double or as an array of
 * two doubles.
 */
static void json_get_range(const char * field, double * range)
{
        jsmntok_t * token = card.tokens + card.cursor;
        if (token->type == JSMN_ARRAY) {
                if (token->size != 2) {
                        fprintf(stderr, string_size_error, __FILE__, __LINE__,
                            token->size, 2, field, card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
                card.cursor++;
                range[0] = json_get_double(field);
                range[1] = json_get_double(field);
        } else
                range[0] = range[1] = json_get_double(field);
}

/* Update the event recorder according to the data card. */
static void card_update_recorder(void)
{
        const char * output_file = json_get_string();
        danton_destroy((void **)&context->recorder);
        context->recorder =
            (struct danton_recorder *)danton_text_create(output_file);
        if (context->recorder == NULL) gracefully_exit(EXIT_FAILURE);
}

/* Update DANTON's run mode according to the data card. */
static void card_update_mode(void)
{
        const char * s = json_get_string();
        if (strcmp(s, "backward") == 0)
                context->mode = DANTON_MODE_BACKWARD;
        else if (strcmp(s, "forward") == 0)
                context->mode = DANTON_MODE_FORWARD;
        else if (strcmp(s, "grammage") == 0)
                context->mode = DANTON_MODE_GRAMMAGE;
        else {
                fprintf(stderr, "%s (%d): invalid mode `%s` in "
                                "JSON file `%s`.\n",
                    __FILE__, __LINE__, s, card.path);
                gracefully_exit(EXIT_FAILURE);
        }
}

/* List of particle names, following DANTON's ordering. */
static const char * particle_name[DANTON_PARTICLE_N] = { "nu_tau~", "nu_mu~",
        "nu_e~", "nu_e", "nu_mu", "nu_tau", "tau~", "tau" };

/* Get a particle name from the card and convert it to a DANTON index. */
static enum danton_particle card_get_particle(enum danton_particle index_max)
{
        const char * particle = json_get_string();
        enum danton_particle index;
        for (index = 0; index < index_max; index++) {
                if (strcmp(particle, particle_name[index]) == 0) return index;
        }
        fprintf(stderr, string_particle_error, __FILE__, __LINE__, particle,
            card.path);
        gracefully_exit(EXIT_FAILURE);
        return DANTON_PARTICLE_UNKNOWN; /* For warnings ... */
}

/* Update the sampler according to the data card. */
static void card_update_sampler(void)
{
        /* Create and initialise the sampler if required. */
        struct danton_sampler * sampler;
        if (context->sampler == NULL) {
                sampler = danton_sampler_create();
                if (sampler == NULL) {
                        fprintf(stderr, "%s\n", danton_error_pop(NULL));
                        gracefully_exit(EXIT_FAILURE);
                }
                sampler->altitude[0] = 1E-03;
                sampler->altitude[1] = 1E+04;
                sampler->elevation[0] = -10.;
                sampler->elevation[1] = 10.;
                sampler->energy[0] = 1E+06;
                sampler->energy[1] = 1E+12;
                context->sampler = sampler;
        } else
                sampler = context->sampler;

        const int size = card.tokens[card.cursor++].size;
        int i;
        for (i = 0; i < size; i++) {
                const char * tag = json_get_string();
                if (strcmp(tag, "altitude") == 0)
                        json_get_range(tag, sampler->altitude);
                else if (strcmp(tag, "elevation") == 0)
                        json_get_range(tag, sampler->elevation);
                else if (strcmp(tag, "energy") == 0)
                        json_get_range(tag, sampler->energy);
                else if (strcmp(tag, "weight") == 0) {
                        const int size_j = json_get_object();
                        int j;
                        for (j = 0; j < size_j; j++) {
                                const int k =
                                    card_get_particle(DANTON_PARTICLE_N);
                                sampler->weight[k] = json_get_double(tag);
                        }
                } else {
                        fprintf(stderr, string_key_error, __FILE__, __LINE__,
                            tag, card.path);
                        gracefully_exit(EXIT_FAILURE);
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
        const int size = json_get_object();
        int i;
        for (i = 0; i < size; i++) {
                const char * field = json_get_string();
                if (strcmp(field, "energy") == 0)
                        json_get_range(field, energy);
                else if (strcmp(field, "exponent") == 0)
                        exponent = json_get_double(field);
                else if (strcmp(field, "weight") == 0)
                        weight = json_get_double(field);
                else {
                        fprintf(stderr, string_key_error, __FILE__, __LINE__,
                            "power-law", card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
        }

        /* Create the model. */
        struct danton_primary * primary =
            (struct danton_primary *)danton_powerlaw_create(
                energy[0], energy[1], exponent, weight);
        if (primary == NULL) {
                fprintf(stderr, "%s\n", danton_error_pop(NULL));
                gracefully_exit(EXIT_FAILURE);
        }
        return primary;
}

/* Update the primary flux with a discrete model. */
static struct danton_primary * card_update_primary_discrete(void)
{
        /* Parse the model parameters. */
        double energy = 1E+12;
        double weight = 1.;
        const int size = json_get_object();
        int i;
        for (i = 0; i < size; i++) {
                const char * field = json_get_string();
                if (strcmp(field, "energy") == 0)
                        energy = json_get_double(field);
                else if (strcmp(field, "weight") == 0)
                        weight = json_get_double(field);
                else {
                        fprintf(stderr, string_key_error, __FILE__, __LINE__,
                            "discrete", card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
        }

        /* Create the model. */
        struct danton_primary * primary =
            (struct danton_primary *)danton_discrete_create(energy, weight);
        if (primary == NULL) {
                fprintf(stderr, "%s\n", danton_error_pop(NULL));
                gracefully_exit(EXIT_FAILURE);
        }
        return primary;
}

/* Update the primary flux according to the data card. */
static void card_update_primary(void)
{
        const int size = card.tokens[card.cursor++].size;
        int i;
        for (i = 0; i < size; i++) {
                /* Parse the particle index. */
                const int index = card_get_particle(DANTON_PARTICLE_N_NU);

                /* Check that one has a length 2 array. */
                jsmntok_t * token = card.tokens + card.cursor;
                if (token->type != JSMN_ARRAY) {
                        fprintf(stderr, string_type_error, __FILE__, __LINE__,
                            token->type, particle_name[index], card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
                if (token->size != 2) {
                        fprintf(stderr, string_size_error, __FILE__, __LINE__,
                            token->size, 2, particle_name[index], card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
                card.cursor++;

                /* Clear any existing model. */
                danton_destroy((void **)&context->primary[index]);

                /* Parse the primary model.  */
                const char * model = json_get_string();
                if (strcmp(model, "power-law") == 0)
                        context->primary[index] =
                            card_update_primary_powerlaw();
                else if (strcmp(model, "discrete") == 0)
                        context->primary[index] =
                            card_update_primary_discrete();
                else {
                        fprintf(stderr, "%s (%d): invalid primary model `%s` "
                                        "for particle `%s` in file `%s`.\n",
                            __FILE__, __LINE__, model, particle_name[index],
                            card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
        }
}

/* Update the Earth model according to the data card. */
static int card_update_pem(void)
{
        /* Default Earth configuration. */
        int pem_sea = 1;

        /* Parse the data card. */
        const int size = json_get_object();
        int i;
        for (i = 0; i < size; i++) {
                const char * field = json_get_string();
                if (strcmp(field, "sea") == 0)
                        pem_sea = json_get_bool(field);
                else {
                        fprintf(stderr, string_key_error, __FILE__, __LINE__,
                            "earth-model", card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
        }
        return pem_sea;
}

/* Update DANTON's configuration according to the content of the data card. */
static void card_update(int * n_events, int * pem_sea)
{
        /* Check that the root is a dictionary. */
        json_get_object();

        /* Loop over the fields. */
        while (card.cursor < card.n_tokens) {
                const char * tag = json_get_string();
                if (strcmp(tag, "events") == 0)
                        *n_events = json_get_int(tag);
                else if (strcmp(tag, "output-file") == 0)
                        card_update_recorder();
                else if (strcmp(tag, "mode") == 0)
                        card_update_mode();
                else if (strcmp(tag, "decay") == 0)
                        context->decay = json_get_bool(tag);
                else if (strcmp(tag, "longitudinal") == 0)
                        context->longitudinal = json_get_bool(tag);
                else if (strcmp(tag, "particle-sampler") == 0)
                        card_update_sampler();
                else if (strcmp(tag, "primary-flux") == 0)
                        card_update_primary();
                else if (strcmp(tag, "earth-model") == 0)
                        *pem_sea = card_update_pem();
                else {
                        fprintf(stderr, string_key_error, __FILE__, __LINE__,
                            tag, card.path);
                        gracefully_exit(EXIT_FAILURE);
                }
        }
}

int main(int argc, char * argv[])
{
        /* If no data card was provided let us show the help message and
         * exit.
         */
        if (argc <= 1) exit_with_help(EXIT_SUCCESS);

        /* Initialise DANTON. */
        if (danton_initialise(NULL, NULL, NULL) != EXIT_SUCCESS) {
                fprintf(stderr, "%s\n", danton_error_pop(NULL));
                exit(EXIT_FAILURE);
        }

        /* Create a simulation context. */
        context = danton_context_create();
        if (context == NULL) {
                fprintf(stderr, "%s\n", danton_error_pop(NULL));
                exit(EXIT_FAILURE);
        }

        /* Load the card from the disk. */
        card_load(argv[1]);

        /* Set the input arguments from the JSON card. */
        int n_events = 10000, pem_sea = 1;
        card_update(&n_events, &pem_sea);

        /* Configure the Earth model. */
        if (!pem_sea) danton_pem_dry();

        /* Update the particle sampler. */
        if (danton_sampler_update(context->sampler) != EXIT_SUCCESS) {
                fprintf(stderr, "%s\n", danton_error_pop(NULL));
                exit(EXIT_FAILURE);
        }

        /* Run the simulation. */
        if (danton_run(context, n_events) != EXIT_SUCCESS)
                fprintf(stderr, "%s\n", danton_error_pop(context));

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_SUCCESS);
}
