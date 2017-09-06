/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C99 executable dedicated to the sampling of decaying
 * taus from ultra high energy neutrinos.
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
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* The DANTON API. */
#include "danton.h"
#include "danton/primary/powerlaw.h"

/* Flux callback for the power law spectrum. */
static double flux(struct danton_primary * primary, double energy)
{
        struct danton_powerlaw * powerlaw = (struct danton_powerlaw *)primary;
        if (energy > 0.) {
                double pdf;
                if (powerlaw->base.energy[0] == powerlaw->base.energy[1])
                        pdf = 0.;
                else if (powerlaw->exponent == -1.) {
                        const double r = log(powerlaw->base.energy[1] /
                            powerlaw->base.energy[0]);
                        pdf = 1. / (r * energy);
                } else {
                        const double a1 = powerlaw->exponent + 1.;
                        const double r = pow(powerlaw->base.energy[1], a1) -
                            pow(powerlaw->base.energy[0], a1);
                        pdf = a1 * pow(energy, powerlaw->exponent) / r;
                }
                return powerlaw->weight * pdf;
        } else
                return powerlaw->weight;
}

/* API function for creating a new power law spectrum. */
struct danton_powerlaw * danton_powerlaw_create(
    double energy_min, double energy_max, double exponent, double weight)
{
        /* Check the arguments. */
        if ((energy_min <= 0.) || (energy_min >= energy_max) || (weight < 0.)) {
                danton_error_push(NULL,
                    "%s (%d): invalid argument(s) (%.5lE, %.5lE).", __FILE__,
                    __LINE__);
                return NULL;
        }

        /* Allocate the memory for the new power law spectrum. */
        struct danton_powerlaw * powerlaw;
        if ((powerlaw = malloc(sizeof(*powerlaw))) == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory.",
                    __FILE__, __LINE__);
                return NULL;
        }

        /* Initialise the power law spectrum and return. */
        powerlaw->base.flux = &flux;
        powerlaw->base.energy[0] = energy_min;
        powerlaw->base.energy[1] = energy_max;
        powerlaw->exponent = exponent;
        powerlaw->weight = weight;

        return powerlaw;
}

/* API function for checking for a power law type. */
int danton_powerlaw_check(struct danton_primary * primary)
{
        return (primary->flux == &flux);
}
