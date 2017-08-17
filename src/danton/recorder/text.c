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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The DANTON API. */
#include "danton.h"
#include "danton/recorder/text.h"

/* Low level data structure for the text recorder. */
struct text_recorder {
        struct danton_text api;
        long last_id;
        char path[];
};

/* Create or append to the output stream. */
static FILE * output_open(
    struct danton_context * context, struct text_recorder * text)
{
        if (text->path[0] == 0x0) return stdout;
        const char * mode =
            (text->api.mode == DANTON_TEXT_MODE_CREATE) ? "w+" : "a+";
        FILE * stream = fopen(text->path, mode);
        if (stream == NULL) {
                danton_error_push(context,
                    "%s (%d): could not open file `%s`\n", __FILE__, __LINE__,
                    text->path);
                return NULL;
        }
        return stream;
}

/* Close the output stream. */
static void output_close(struct text_recorder * text, FILE * stream)
{
        if (text->path[0] != 0x0) fclose(stream);
}

/* Format the header for a Monte-Carlo event. */
static void format_header_event(FILE * stream)
{
        fprintf(stream,
            "    Event   PID    Energy             Direction or "
            "Momentum                         Position               "
            "      Weight\n                    (GeV)                "
            " (1 or GeV/c)                                 (m)\n     "
            "                           ux or Px     uy or Py    "
            "uz or Pz         X             Y             Z\n");
}

/* Utility function for formating a state to the output stream. */
static void format_state(FILE * stream, long index, const int * pid,
    const struct danton_state * state, const double * weight)
{
        if (index > 0)
                fprintf(stream, "%10ld ", index);
        else
                fprintf(stream, "%10c ", ' ');
        if (pid != NULL)
                fprintf(stream, "%4d ", *pid);
        else
                fprintf(stream, "%4c ", ' ');
        fprintf(stream, "%12.5lE %12.5lE %12.5lE %12.5lE %13.3lf "
                        "%13.3lf %13.3lf",
            state->energy, state->direction[0], state->direction[1],
            state->direction[2], state->position[0], state->position[1],
            state->position[2]);
        if (weight != NULL)
                fprintf(stream, " %12.5lE\n", *weight);
        else
                fprintf(stream, "\n");
}

/* Utility function for formating a decay product to the output stream. */
static void format_product(FILE * stream, const struct danton_product * product)
{
        fprintf(stream, "%10c %4d %12c %12.5lE %12.5lE %12.5lE\n", ' ',
            product->pid, ' ', product->momentum[0], product->momentum[1],
            product->momentum[2]);
}

/* Format the header for a grammage data. */
static void format_header_grammage(FILE * stream)
{
        fprintf(stream, "   elevation    Grammage\n     (deg)      (kg/m^2)\n");
}

/* Format a grammage data point. */
static void format_grammage(
    FILE * stream, const struct danton_grammage * grammage)
{
        fprintf(
            stream, "%12.5lf %12.5lE\n", grammage->elevation, grammage->value);
}

/* Callback for recording an event. */
static int record_event(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_event * event)
{
        /* Unpack the text recorder object. */
        struct text_recorder * text = (struct text_recorder *)recorder;

        /* Open the output stream. */
        FILE * stream = output_open(context, text);
        if (stream == NULL) return EXIT_FAILURE;

        /* Print the header if the file is new. */
        if (text->api.mode == DANTON_TEXT_MODE_CREATE) {
                format_header_event(stream);

                /* Let us go back to append mode for next calls. */
                text->api.mode = DANTON_TEXT_MODE_APPEND;
        }

        /* Dump the Monte-Carlo states. Make sure that the event id,
         * the event weight and the generation index are printed once and
         * only once.
         */
        long index[3] = { event->id + 1, event->generation, 0 };
        const double * weight[3] = { &event->weight, NULL, NULL };
        int i = 0;
        if ((event->primary != NULL) && (event->id != text->last_id)) {
                format_state(stream, index[i], &event->primary->pid,
                    event->primary, weight[i]);
                i++;
                text->last_id = event->id;
        }
        int dump_pid = 1;
        if (event->vertex != NULL) {
                format_state(stream, index[i], &event->vertex->pid,
                    event->vertex, weight[i]);
                dump_pid = 0;
                i++;
        }
        int * pid = dump_pid ? &event->final->pid : NULL;
        format_state(stream, index[i], pid, event->final, weight[i]);

        /* Dump the decay products. */
        struct danton_product * p;
        for (i = 0, p = event->product; i < event->n_products; i++, p++)
                format_product(stream, p);

        /* Close the output stream. */
        output_close(text, stream);

        return EXIT_SUCCESS;
}

/* Callback for recording a grammage point. */
static int record_grammage(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_grammage * grammage)
{
        /* Unpack the text recorder object. */
        struct text_recorder * text = (struct text_recorder *)recorder;

        /* Open the output stream. */
        FILE * stream = output_open(context, text);
        if (stream == NULL) return EXIT_FAILURE;

        /* Print the header if the file is new. */
        if (text->api.mode == DANTON_TEXT_MODE_CREATE) {
                format_header_grammage(stream);

                /* Let us go back to append mode for next calls. */
                text->api.mode = DANTON_TEXT_MODE_APPEND;
        }

        /* Format the grammage data point. */
        format_grammage(stream, grammage);

        /* Close the output stream. */
        output_close(text, stream);

        return EXIT_SUCCESS;
}

/* API function for creating a new text recorder. */
struct danton_text * danton_text_create(const char * path)
{
        /* Allocate the memory for the new text recorder. */
        const int n = (path == NULL) ? 1 : strlen(path) + 1;
        struct text_recorder * text;
        if ((text = malloc(sizeof(*text) + n)) == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory\n",
                    __FILE__, __LINE__);
                return NULL;
        }

        /* Initialise the text recorder and return. */
        text->api.base.record_event = &record_event;
        text->api.base.record_grammage = &record_grammage;
        text->api.mode = DANTON_TEXT_MODE_CREATE;
        text->last_id = -1;
        if (n > 1)
                memcpy(text->path, path, n);
        else
                text->path[0] = 0x0;

        return &text->api;
}

/* API function for checking for a text recorder type. */
int danton_text_check(struct danton_recorder * recorder)
{
        return ((recorder->record_event == &record_event) &&
            (recorder->record_grammage == &record_grammage));
}
