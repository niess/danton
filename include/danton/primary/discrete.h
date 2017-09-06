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

struct danton_discrete;

DANTON_API struct danton_discrete * danton_discrete_create(
    double energy, double weight);

DANTON_API int danton_discrete_set(
    struct danton_discrete * discrete, double energy, double weight);
DANTON_API void danton_discrete_get(
    const struct danton_discrete * discrete, double * energy, double * weight);

DANTON_API int danton_discrete_check(struct danton_primary * primary);

#ifdef __cplusplus
}
#endif
#endif
