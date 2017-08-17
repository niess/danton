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

#ifndef danton_h
#define danton_h
#ifdef __cplusplus
extern "C" {
#endif

#ifndef DANTON_API
#define DANTON_API
#endif

/** Indices of DANTON particles. */
enum danton_particle {
        DANTON_PARTICLE_UNKNOWN = -1,
        DANTON_PARTICLE_NU_BAR_TAU = 0,
        DANTON_PARTICLE_NU_BAR_MU,
        DANTON_PARTICLE_NU_BAR_E,
        DANTON_PARTICLE_NU_E,
        DANTON_PARTICLE_NU_MU,
        DANTON_PARTICLE_NU_TAU,
        /** The total number of neutrino indices. */
        DANTON_PARTICLE_N_NU,
        DANTON_PARTICLE_TAU_BAR = DANTON_PARTICLE_N_NU,
        DANTON_PARTICLE_TAU,
        /** The total number of indices. */
        DANTON_PARTICLE_N
};

/* Data structures for the primary neutrino flux. */
struct danton_primary;

typedef double danton_primary_cb(
    struct danton_primary * primary, double energy);

struct danton_primary {
        danton_primary_cb * flux;
        double energy[2];
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

struct danton_event {
        long id;
        double weight;
        struct danton_state * primary;
        int generation;
        struct danton_state * vertex;
        struct danton_state * final;
        int n_products;
        struct danton_product * product;
};

/* Data structure for storing a grammage computation. */
struct danton_grammage {
        double elevation;
        double value;
};

/* Callback for recording a sampled event. */
struct danton_context;
struct danton_recorder;
typedef int danton_event_cb(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_event * event);

/* Callback for recording a grammage value. */
typedef int danton_grammage_cb(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_grammage * grammage);

/* Generic type for recording the sampled data. */
struct danton_recorder {
        danton_event_cb * record_event;
        danton_grammage_cb * record_grammage;
};

/* Available run modes. */
enum danton_mode {
        DANTON_MODE_BACKWARD = 0,
        DANTON_MODE_FORWARD,
        DANTON_MODE_GRAMMAGE,
        /** The total number of run modes. */
        DANTON_MODE_N
};

/* Handle for a simulation context. */
struct danton_context {
        enum danton_mode mode;
        int longitudinal;
        int decay;
        struct danton_primary * primary[DANTON_PARTICLE_N_NU];
        struct danton_sampler * sampler;
        struct danton_recorder * recorder;
};

/* Generic lock / unlock callback. */
typedef int danton_lock_cb(void);

DANTON_API int danton_initialise(
    const char * pdf, danton_lock_cb * lock, danton_lock_cb * unlock);
DANTON_API void danton_finalise(void);
DANTON_API void danton_destroy(void ** any);

/* Replace the sea layer of the PEM with Standard Rock. */
DANTON_API void danton_pem_dry(void);

/* Get the PDG particle number for a given particle index. */
DANTON_API int danton_particle_pdg(enum danton_particle particle);

/* Get the DANTON particle index for a given PDG number. */
DANTON_API enum danton_particle danton_particle_index(int pdg);

DANTON_API struct danton_sampler * danton_sampler_create(void);
DANTON_API int danton_sampler_update(struct danton_sampler * sampler);

DANTON_API struct danton_context * danton_context_create(void);
DANTON_API void danton_context_destroy(struct danton_context ** context);

DANTON_API int danton_run(struct danton_context * context, long events);

DANTON_API int danton_error_count(struct danton_context * context);
DANTON_API const char * danton_error_pop(struct danton_context * context);
DANTON_API int danton_error_push(
    struct danton_context * context, const char * format, ...);

#ifdef __cplusplus
}
#endif
#endif
