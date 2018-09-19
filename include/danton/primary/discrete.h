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

#ifndef danton_discrete_h
#define danton_discrete_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DANTON_API
#define DANTON_API
#endif

#include "danton.h"

/** Opaque structure for a discrete danton_primary. */
struct danton_discrete;

/**
 * Create a discrete danton_primary.
 *
 * @param  energy   The energy of the primary neutrino, in GeV.
 * @param  weight   The intensity of the discrete source.
 * @return          The corresponding discrete danton_primary, or `ǸULL`.
 */
DANTON_API struct danton_discrete * danton_discrete_create(
    double energy, double weight);

/**
 * Set the properties of a discrete danton_primary.
 *
 * @param discrete  The discrete danton_primary.
 * @param  energy   The energy of the primary neutrino, in GeV.
 * @param  weight   The intensity of the discrete source.
 * @return          `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
DANTON_API int danton_discrete_set(
    struct danton_discrete * discrete, double energy, double weight);

/**
 * Get the properties of a discrete danton_primary.
 *
 * @param discrete  The discrete danton_primary.
 * @param  energy   The energy of the primary neutrino, in GeV.
 * @param  weight   The intensity of the discrete source.
 */
DANTON_API void danton_discrete_get(
    const struct danton_discrete * discrete, double * energy, double * weight);

/**
 * Check if a danton_primary is of *discrete* type.
 *
 * @param primary  The danton_primary.
 * @return         `1` if the primary is a discrete one, `0` otherwise.
 */
DANTON_API int danton_discrete_check(struct danton_primary * primary);

#ifdef __cplusplus
}
#endif
#endif
