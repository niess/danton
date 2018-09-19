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

#ifndef danton_text_h
#define danton_text_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DANTON_API
#define DANTON_API
#endif

#include "danton.h"

/** Operations mode for a text danton_recorder. */
enum danton_text_mode {
        /** Append to an existing text file, or create it. */
        DANTON_TEXT_MODE_APPEND = 0,
        /** Create a new text file, or override it. */
        DANTON_TEXT_MODE_CREATE,
        /** The number of operations modes.  */
        DANTON_TEXT_MODE_N
};

/**
 * Data structure for a text danton_recorder.
 *
 * This is an implementation of a danton_recorder to a *text* file. The exposed
 * data can be directly modified.
 */
struct danton_text {
        /** The base danton_recorder. */
        struct danton_recorder base;
        /** The operation mode. */
        enum danton_text_mode mode;
};

/**
 * Create a text danton_recorder.
 *
 * @param  path  The path to the text file
 * @return       The corresponding text danton_recorder, or `ǸULL`.
 */
DANTON_API struct danton_text * danton_text_create(const char * path);

/**
* Check if a danton_recorder is of *text* type.
*
* @param primary  The danton_recorder.
* @return         `1` if the recorder is a text one, `0` otherwise.
*/
DANTON_API int danton_text_check(struct danton_recorder * recorder);

#ifdef __cplusplus
}
#endif
#endif
