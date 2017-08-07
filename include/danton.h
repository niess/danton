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

struct danton_primary;

typedef double danton_primary_cb(
    struct danton_primary * primary, double energy);

struct danton_primary {
        danton_primary_cb * flux;
        int pid;
        double energy_min;
        double energy_max;
};

struct danton_state {
        int pid;
        double energy;
        double position[3];
        double direction[3];
};

struct danton_product {
        int pid;
        double momentum[3];
};

struct danton_decay {
        int generation;
        struct danton_state at_decay;
        struct danton_state at_production;
        int n_products;
        struct danton_product * product;
};

struct danton_event {
        long id;
        double weight;
        struct danton_state state;
        int n_decays;
        struct danton_decay * decay;
};

struct danton_context {
        int forward;
        int longitudinal;
        int decay;
        int grammage;
        const char * output;
};

typedef void danton_event_cb(
    struct danton_context * context, struct danton_event * event);

typedef int danton_lock_cb(void);

int danton_initialise(
    const char * pdf, danton_lock_cb * lock, danton_lock_cb * unlock);

void danton_finalise(void);

/* Replace the sea layer of the PEM with Standard Rock. */
void danton_pem_dry(void);

struct danton_context * danton_context_create(void);
void danton_context_destroy(struct danton_context * context);

void danton_altitude(struct danton_context * context, double altitude);
void danton_altitude_range(
    struct danton_context * context, double altitude_min, double altitude_max);

void danton_energy(struct danton_context * context, double energy);
void danton_energy_range(
    struct danton_context * context, double energy_min, double energy_max);

void danton_elevation(struct danton_context * context, double elevation);
void danton_elevation_range(struct danton_context * context,
    double elevation_min, double elevation_max);

void danton_cos_theta(struct danton_context * context, double cos_theta);
void danton_cos_theta_range(struct danton_context * context,
    double cos_theta_min, double cos_theta_max);

void danton_target_clear(struct danton_context * context);
void danton_target_set(struct danton_context * context, int pid, double weight);
void danton_target_unset(struct danton_context * context, int pid);
int danton_target_active(struct danton_context * context, int pid);
double danton_target_weight(struct danton_context * context, int pid);

int danton_run(struct danton_context * context, long events);

const char * danton_error(struct danton_context * context);
