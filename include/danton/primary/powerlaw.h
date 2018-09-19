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

#ifndef danton_powerlaw_h
#define danton_powerlaw_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DANTON_API
#define DANTON_API
#endif

#include "danton.h"

/**
 * Data structure for a powerlaw danton_primary.
 *
 * This is an implementation of a powerlaw danton_primary. The exposed
 * data can be directly modified.
 */
struct danton_powerlaw {
        /** The base danton_primary object. */
        struct danton_primary base;
        /** The exponent of the power law. */
        double exponent;
        /** The weight of the primary source. */
        double weight;
};

/**
 * Create a powerlaw danton_primary.
 *
 * @param  energy_min   The minimum energy of the primary neutrino, in GeV.
 * @param  energy_max   The maximum energy of the primary neutrino, in GeV.
 * @param  exponent     The exponent of the powerlaw.
 * @param  weight       The intensity of the powerlaw source.
 * @return              The corresponding powerlaw danton_primary, or `ǸULL`.
 */
DANTON_API struct danton_powerlaw * danton_powerlaw_create(
    double energy_min, double energy_max, double exponent, double weight);

/**
 * Check if a danton_primary is of *powerlaw* type.
 *
 * @param primary  The danton_primary.
 * @return         `1` if the primary is a powerlaw one, `0` otherwise.
 */
DANTON_API int danton_powerlaw_check(struct danton_primary * primary);

#ifdef __cplusplus
}
#endif
#endif
