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

/* Data structures for the primary neutrino flux. */
struct danton_primary;

typedef double danton_primary_cb(
    struct danton_primary * primary, double energy);

struct danton_primary {
        danton_primary_cb * flux;
        int pid;
        double energy_min;
        double energy_max;
};

/** Indices of DANTON particles. */
enum danton_particle {
        DANTON_PARTICLE_UNKNOWN = -1,
        DANTON_PARTICLE_NU_BAR_TAU = 0,
        DANTON_PARTICLE_NU_BAR_MU,
        DANTON_PARTICLE_NU_BAR_E,
        DANTON_PARTICLE_NU_E,
        DANTON_PARTICLE_NU_MU,
        DANTON_PARTICLE_NU_TAU,
        DANTON_PARTICLE_TAU_BAR,
        DANTON_PARTICLE_TAU,
        /** The total number of indices. */
        DANTON_PARTICLE_N
};

/* Handle for managing the events to sample. */
struct danton_sampler {
        double altitude[2];
        double cos_theta[2];
        double elevation[2];
        double energy[2];
        double weight[DANTON_PARTICLE_N];
};

/* Data structures relative to the sampled event. */
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

/* Handle for a simulation context. */
struct danton_context {
        int forward;
        int longitudinal;
        int decay;
        int grammage;
        struct danton_sampler * sampler;
        const char * output;
};

/* Callback for processing a sampled event. */
typedef void danton_event_cb(
    struct danton_context * context, struct danton_event * event);

/* Generic lock / unlock callback. */
typedef int danton_lock_cb(void);

int danton_initialise(
    const char * pdf, danton_lock_cb * lock, danton_lock_cb * unlock);
void danton_finalise(void);

/* Replace the sea layer of the PEM with Standard Rock. */
void danton_pem_dry(void);

/* Get the PDG particle number for a given particle index. */
int danton_particle_pdg(enum danton_particle particle);

/* Get the DANTON particle index for a given PDG number. */
enum danton_particle danton_particle_index(int pdg);

struct danton_sampler * danton_sampler_create(void);
void danton_sampler_destroy(struct danton_sampler ** sampler);
int danton_sampler_update(struct danton_sampler * sampler);

struct danton_context * danton_context_create(void);
void danton_context_destroy(struct danton_context ** context);

int danton_run(struct danton_context * context, long events);

const char * danton_error(struct danton_context * context);
