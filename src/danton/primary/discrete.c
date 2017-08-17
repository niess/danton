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

/* The DANTON API. */
#include "danton.h"
#include "danton/primary/discrete.h"

/* Data container for a discrete spectrum. */
struct danton_discrete {
        struct danton_primary base;
        double weight;
};

/* Flux callback for the discrete spectrum. */
static double flux(struct danton_primary * primary, double energy)
{
        if (energy > 0.)
                return 0.;
        else {
                struct danton_discrete * discrete =
                    (struct danton_discrete *)primary;
                return discrete->weight;
        }
}

/* API function for setting the parameters of a discrete spectrum. */
void danton_discrete_set(
    struct danton_discrete * discrete, double energy, double weight)
{
        /* TODO: check the parameter values. */
        discrete->base.energy[0] = discrete->base.energy[1] = energy;
        discrete->weight = weight;
}

/* API function for getting the parameters of a discrete spectrum. */
void danton_discrete_get(
    const struct danton_discrete * discrete, double * energy, double * weight)
{
        if (energy != NULL) *energy = discrete->base.energy[0];
        if (weight != NULL) *weight = discrete->weight;
}

/* API function for creating a new discrete spectrum. */
struct danton_discrete * danton_discrete_create(double energy, double weight)
{
        /* Allocate the memory for the new discrete spectrum. */
        struct danton_discrete * discrete;
        if ((discrete = malloc(sizeof(*discrete))) == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory\n",
                    __FILE__, __LINE__);
                return NULL;
        }

        /* Initialise the discrete spectrum and return. */
        discrete->base.flux = &flux;
        danton_discrete_set(discrete, energy, weight);
        return discrete;
}

/* API function for checking for a discrete type. */
int danton_discrete_check(struct danton_primary * primary)
{
        return (primary->flux == &flux);
}
