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

#ifndef danton_text_h
#define danton_text_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DANTON_API
#define DANTON_API
#endif

#include "danton.h"

enum danton_text_mode {
        DANTON_TEXT_MODE_APPEND = 0,
        DANTON_TEXT_MODE_CREATE,
        DANTON_TEXT_MODE_N
};

struct danton_text {
        struct danton_recorder base;
        enum danton_text_mode mode;
};

DANTON_API struct danton_text * danton_text_create(const char * path);
DANTON_API void danton_text_destroy(struct danton_text ** text);

#ifdef __cplusplus
}
#endif
#endif
