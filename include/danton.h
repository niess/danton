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
        /**  */
        DANTON_PARTICLE_UNKNOWN = -1,
        /** The anti tau neutrino. */
        DANTON_PARTICLE_NU_BAR_TAU = 0,
        /** The anti muon neutrino. */
        DANTON_PARTICLE_NU_BAR_MU,
        /** The anti electron neutrino. */
        DANTON_PARTICLE_NU_BAR_E,
        /** The electron neutrino. */
        DANTON_PARTICLE_NU_E,
        /** The muon neutrino. */
        DANTON_PARTICLE_NU_MU,
        /** The tau neutrino. */
        DANTON_PARTICLE_NU_TAU,
        /** The total number of neutrino indices. */
        DANTON_PARTICLE_N_NU,
        /** The anti tau. */
        DANTON_PARTICLE_TAU_BAR = DANTON_PARTICLE_N_NU,
        /** The tau. */
        DANTON_PARTICLE_TAU,
        /** The total number of indices. */
        DANTON_PARTICLE_N
};

struct danton_primary;
/**
 * Callback for a primary neutrino flux.
 *
 * @param  primary  Handle for the primary data.
 * @param  energy   The energy of the primary neutrino, in GeV.
 * @return          The corresponding flux.
 */
typedef double danton_primary_cb(
    struct danton_primary * primary, double energy);

/**
 * Base data for a primary flux model.
 *
 * This is the minimalistic representation of a primary neutrino flux. Any
 * model implementation should conform to this.
 */
struct danton_primary {
        /** The primary flux callback. */
        danton_primary_cb * flux;
        /** The energy range over which the model is defined, in GeV. */
        double energy[2];
};

/**
 * Data container for an event sampler.
 *
 * This structure specifies the final state of the events to sample. Once filled
 * one **must** call `danton_sampler_update` for the sampler to be ready to
 * use.
 */
struct danton_sampler {
        /** The altitude range over which events are sampled, in m. */
        double altitude[2];
        /** The elevation angle range of the final state, in deg. */
        double elevation[2];
        /** The energy range of the final state, in GeV. */
        double energy[2];
        /**
         * A weight vector specifying the sampling probabilities, in [0,1],
         * for all particles.
         */
        double weight[DANTON_PARTICLE_N];
};

/** Data container for a recorded particle state. */
struct danton_state {
        /** The PDG ID of the recorded particle.  */
        int pid;
        /** The total energy of the recorded particle, in GeV. */
        double energy;
        /** The ECEF position of the recorded particle, in m. */
        double position[3];
        /** The momentum's direction of the recorded particle, in ECEF frame. */
        double direction[3];
};

/** Data container for a tau decay product. */
struct danton_product {
        /** The PDG ID of the daughter. */
        int pid;
        /** The 3-momentum of the daughter, in GeV/c. */
        double momentum[3];
};

/** Data container for exposing a recorded event. */
struct danton_event {
        /** The Monte-Carlo index of the event. */
        long id;
        /** The Monte-Carlo weight for the event. */
        double weight;
        /** Pointer to the primary state. */
        struct danton_state * primary;
        /**
         * The generation index of the recorded tau. An index > 1 indicates
         * that the tau originates from a regenerated neutrino, not from a
         * primary one.
         */
        int generation;
        /** Pointer to the tau state at its creation vertex. */
        struct danton_state * vertex;
        /** Pointer to the final particle state, before decay for taus. */
        struct danton_state * final;
        /** The number of decay products. */
        int n_products;
        /** Array of decay products. */
        struct danton_product * product;
};

/** Data container for exposing a grammage computation. */
struct danton_grammage {
        /** The elevation angle of observation, in deg. */
        double elevation;
        /** The corresponding grammage value, in kg/m^3. */
        double value;
};

struct danton_context;
struct danton_recorder;
/** Callback for recording a sampled event.
 *
 * @param  context   Handle for the simulation context.
 * @param  recorder  Handle for the recorder data.
 * @param  event     Handle for the event data.
 * @return           `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
typedef int danton_event_cb(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_event * event);

/** Callback for recording a grammage value.
 *
 * @param  context   Handle for the simulation context.
 * @param  recorder  Handle for the recorder data.
 * @param  grammage  Handle for the grammage data.
 * @return          `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
typedef int danton_grammage_cb(struct danton_context * context,
    struct danton_recorder * recorder, const struct danton_grammage * grammage);

/**
 * Base data for recording the sampled events or grammage computations.
 *
 * This is the minimalistic representation of an event recorder. Any
 * concrete implementation should conform to this.
 */
struct danton_recorder {
        danton_event_cb * record_event;
        danton_grammage_cb * record_grammage;
};

/** The available run modes. */
enum danton_mode {
        /** Backward Monte-Carlo simulation. */
        DANTON_MODE_BACKWARD = 0,
        /** Classical forward Monte-Carlo. */
        DANTON_MODE_FORWARD,
        /** Grammage computation. */
        DANTON_MODE_GRAMMAGE,
        /** The total number of run modes. */
        DANTON_MODE_N
};

/** Handle for a simulation context.
 *
 * This structure is a proxy to thread specific simulation data. It exposes
 * some public data that the user may configure or alter directly.
 * However, it also encloses other opaque data. Therefore, it **must** be
 * initialised and released with the `danton_context` functions.
 */
struct danton_context {
        /** The run mode. */
        enum danton_mode mode;
        /**
         * Flag for controlling the tranverse transport.
         *
         * Set this flag to `1` in order to disable tranverse transport, i.e.
         * force straight trajectories. By default the transverse transport is
         * enabled.
         */
        int longitudinal;
        /**
         * Flag for controlling the decay of tau final states.
         *
         * Set this flag to `1` if the sampled tau final states must be decayed.
         * By default the taus are decayed.
         */
        int decay;
        /**
         * Array of pointers to the primary flux models for each neutrino
         * flavour.
         *
         * Set an entry to `ǸULL` in order to disable the corresponding
         * neutrino as a primary. **Note** that in Monte-Carlo modes (Forward,
         * or Backward) at least one entry **must** be non `ǸULL.
         */
        struct danton_primary * primary[DANTON_PARTICLE_N_NU];
        /**
         * Handle for the event sampler.
         *
         * Starts initialised to `ǸULL`. One must provide a sampler before
         * running.
         */
        struct danton_sampler * sampler;
        /**
         * Handle for the event recorder.
         *
         * Starts initialised to `ǸULL`. One must provide a recorder before
         * running.
         */
        struct danton_recorder * recorder;
};

/**
 * Generic lock or unlock callback.
 *
 * @return  `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * For multithreaded usage the user must supply a pair of lock and unlock
 * callbacks ensuring exclusive access to critical data by the simulation
 * contexts.
 */
typedef int danton_lock_cb(void);

/**
 * Initialise the DANTON library.
 *
 * @param  pdf     Path to the PDF tables for ENT, or `NULL`.
 * @param  lock    A locking callback, or `NULL`.
 * @param  unlock  An unlocking callback, or `NULL`.
 * @return         `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Initialise the DANTON library. If no *pdf* path is provided the CT14 NNLO
 * tables are used, shiped with ENT.
 *
 * __Warning__ : for multithreaded usage one **must** provide valid *lock* and
 * *unlock* callbacks.
 */
DANTON_API int danton_initialise(
    const char * pdf, danton_lock_cb * lock, danton_lock_cb * unlock);

/**
 * Finalise the danton library.
 *
 * Release the memory used by the Physics engines. Note that user created
 * objects, e.g. simulation *contexts*, are **not** dealocated.
 */
DANTON_API void danton_finalise(void);

/**
 * Properly dealocate a memory flat object.
 *
 * @param  any  The address of the object pointer.
 *
 * Note that at return the pointer is set to `NULL`.
 */
DANTON_API void danton_destroy(void ** any);

/**
 * Replace the sea layer of the PEM with Standard Rock.
 */
DANTON_API void danton_pem_dry(void);

/**
 * Get the PDG particle number for a given DANTON particle index.
 *
 * @param  particle  The DANTON particle index.
 * @return           The corresponding PDG number.
 */
DANTON_API int danton_particle_pdg(enum danton_particle particle);

/**
 * Get the DANTON particle index for a given PDG number.
 *
 * @param  pdg  The PDG number.
 * @return      The corresponding DANTON index.
 */
DANTON_API enum danton_particle danton_particle_index(int pdg);

/**
 * Create a particle sampler.
 * @return  A handle for the sampler or `NULL` on failure.
 */
DANTON_API struct danton_sampler * danton_sampler_create(void);

/**
 * Update a particle sampler according to its current settings.
 *
 * @param  sampler  The handler for the particle sampler.
 * @return          `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Note that a sampler must be updated **before** it can be used.
 */
DANTON_API int danton_sampler_update(struct danton_sampler * sampler);

/**
 * Create a simulation context.
 *
 * @return  A handle for the context or `NULL` on failure.
 *
 * The simulation context is initialised empty. It **must** be assigned a
 * *sampler* and a *recorder*. For Monte-Carlo runs a *primary* flux model is
 * also required.
 */
DANTON_API struct danton_context * danton_context_create(void);

/**
 * Destroy a simulation context.
 *
 * @param  context  A handle for the context.
 *
 * This function must be used instead of `danton_destroy` in order to properly
 * destroy a simulation context.
 */
DANTON_API void danton_context_destroy(struct danton_context ** context);

/**
 * Run a Monte-Carlo simulation or a grammage scan.
 *
 * @param  context      The simulation context to use.
 * @param  events       The maximum number of Monte-carlo events or scan points.
 * @param  requested    The number of requested events to log.
 * @return              `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Depending on the context *mode* a Monte-Carlo simulation or a grammage
 * scan is done. Note that setting *requested* to zero or less ignores this
 * option, resulting in all events to be processed.
 */
DANTON_API int danton_run(
    struct danton_context * context, long events, long requested);

/**
 * Get the current number of unprocessed errors.
 *
 * @param  context  The context of interest or `NULL`.
 * @return          The corresponding error count.
 *
 * If no *context* is provided the number of global errors is returned.
 */
DANTON_API int danton_error_count(struct danton_context * context);

/**
 * Retrieve the last error in the stack.
 *
 * @param  context  The context of interest or `NULL`.
 * @return          A string describing the last error.
 *
 * If no *context* is provided the error is taken from the global stack.
 */
DANTON_API const char * danton_error_pop(struct danton_context * context);

/**
 * Append an error message to the stack.
 *
 * @param  context  The context where the error occurred or `NULL`.
 * @param  format   A format string for the error message.
 * @param  VARARGS  Variable arguments matching the format.
 * @return          `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * The format string and variable arguments follow the printf syntax.
 */
DANTON_API int danton_error_push(
    struct danton_context * context, const char * format, ...);

#ifdef __cplusplus
}
#endif
#endif
