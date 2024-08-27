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

/** Indices of Danton particles. */
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
 * @return The corresponding flux.
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
        /** The latitude at which events are sampled, in deg. */
        double latitude;
        /** The longitude at which events are sampled, in deg. */
        double longitude;
        /** The altitude range over which events are sampled, in m. */
        double altitude[2];
        /** The azimuth angle range of the final state, in deg. */
        double azimuth[2];
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
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
typedef int danton_event_cb(
    struct danton_context * context,
    struct danton_recorder * recorder,
    const struct danton_event * event);

/** Callback for recording a grammage value.
 *
 * @param  context   Handle for the simulation context.
 * @param  recorder  Handle for the recorder data.
 * @param  grammage  Handle for the grammage data.
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
typedef int danton_grammage_cb(
    struct danton_context * context,
    struct danton_recorder * recorder,
    const struct danton_grammage * grammage);

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

/** Flags for events occuring during custom run action(s). */
enum danton_run_event {
        /** Flag for the start of a new event. */
        DANTON_RUN_EVENT_START = 0,
        /** Flag for the end of the current event. */
        DANTON_RUN_EVENT_STOP,
        /** Flag for a new step in the current event. */
        DANTON_RUN_EVENT_STEP,
        /** The total number of event flags. */
        DANTON_RUN_N
};

struct danton_run_action;
/**
 * Callback for custom run actions.
 *
 * @param  context Handle for the simulation context.
 * @param  event   Flag signing the occuring event.
 * @param  state   The current particle state.
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 */
typedef int danton_run_cb(
    struct danton_context * context,
    struct danton_run_action * run_action,
    enum danton_run_event event,
    int medium,
    struct danton_state * state);

/**
 * Base data for custom run actions.
 *
 * This is the minimalistic representation of a run action. Any
 * concrete implementation should conform to this.
 */
struct danton_run_action {
        danton_run_cb * call;
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
         * Set an entry to `NULL` in order to disable the corresponding
         * neutrino as a primary. **Note** that in Monte-Carlo modes (Forward,
         * or Backward) at least one entry **must** be non `NULL`.
         */
        struct danton_primary * primary[DANTON_PARTICLE_N_NU];
        /**
         * Handle for the event sampler.
         *
         * Starts initialised to `NULL`. One must provide a sampler before
         * running.
         */
        struct danton_sampler * sampler;
        /**
         * Handle for the event recorder.
         *
         * Starts initialised to `NULL`. One must provide a recorder before
         * running.
         */
        struct danton_recorder * recorder;
        /**
         * Handle for custom run action(s).
         *
         * Starts initialised to `NULL`, i.e. disabled. Providing additional
         * run action(s) is optionnal.
         */
        struct danton_run_action * run_action;
};

/**
 * Generic lock or unlock callback.
 *
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * For multithreaded usage the user must supply a pair of lock and unlock
 * callbacks ensuring exclusive access to critical data by the simulation
 * contexts.
 */
typedef int danton_lock_cb(void);

/**
 * Initialise the DANTON library.
 *
 * @param  prefix  Library installation prefix, or `NULL`.
 * @param  lock    A locking callback, or `NULL`.
 * @param  unlock  An unlocking callback, or `NULL`.
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Initialise the DANTON library. The installation *prefix* indicates the
 * location of physics data. If *prefix* is `NULL`, then the value provided at
 * compile time is used instead.
 *
 * See the `danton_finalise` callback for releasing the library temporary
 * memory.
 *
 * __Note__ : for multithreaded usage one **must** provide valid *lock* and
 * *unlock* callbacks. Otherwise those can be set to `NULL`.
 *
 * __Warning__ : this function must be called before using any other library
 * functions.
 */
DANTON_API int danton_initialise(
    const char * prefix,
    danton_lock_cb * lock,
    danton_lock_cb * unlock);

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
 * Set or update the global Earth model.
 *
 * @param geoid      The reference geoid for the sea level, or `NULL`.
 * @param topography Topography model, path to any topographic data, or `NULL`.
 * @param material   Material for the topography.
 * @param density    Density of the topography material in kg / m^(3).
 * @param sea        Pointer to a flag to enable or disable sea(s), or `NULL`.
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * __Warning__ : this function is **not** thread safe. It sets the Earth
 * model globally.
 *
 * The default Earth model is the Preliminary Reference Earth Model, i.e. a
 * spherical Earth fully covered with a 3km deep sea.
 *
 * The supported values for the *geoid* are:
 *
 *     - PREM  : spherical Earth
 *
 *     - WGS84 : GPS ellipsoid
 *
 *     - EGM96 : GPS ellipsoid + geoid undulations
 *
 * If *geoid* is left `NULL`, then the PREM model is assumed.
 *
 * The *topography* parameter can either specify a path to topographic data
 * (e.g. SRTMGL1 tiles) or encode a flat topography as `flat://${z}`, where
 * `${z}` is the constant altitude above sea level, in meters, e.g.
 * `flat://1000` for a 1km high cover.
 *
 * __Note__ : specifying a detailed topography requires the WGS84 or EGM96
 * reference system to be used.
 *
 * The topography *material* must match one of the materials defined in the xml
 * PUMAS Materials Definition File (MDF). If `NULL` is given the material is
 * left unchanged. It defaults to Rock.
 *
 * If a null or negative density is provided the material density is left
 * unchanged. It defaults to 2.650 kg / m^3.
 */
DANTON_API int danton_earth_model(
    const char * geoid,
    const char * topography,
    double density,
    int * sea);

/**
 * Get the current topography.
 * @return A `turtle_stack` pointer, or `NULL`.
 *
 * This routines provides access to the topography used by DANTON as
 * a`turtle_stack` object.  It can be used, e.g. for querying the topography
 * ground altitude with TURTLE.
 *
 * __Warnings__ : the stack should not be modified directly with TURTLE. Use
 * the DANTON API functions instead.
 */
void * danton_get_topography(void);

/**
 * Get the PDG particle number for a given DANTON particle index.
 *
 * @param  particle  The DANTON particle index.
 * @return The corresponding PDG number.
 */
DANTON_API int danton_particle_pdg(enum danton_particle particle);

/**
 * Get the DANTON particle index for a given PDG number.
 *
 * @param  pdg  The PDG number.
 * @return The corresponding DANTON index.
 */
DANTON_API enum danton_particle danton_particle_index(int pdg);

/**
 * Create a particle sampler.
 * @return A handle for the sampler or `NULL` on failure.
 */
DANTON_API struct danton_sampler * danton_sampler_create(void);

/**
 * Update a particle sampler according to its current settings.
 *
 * @param  sampler  The handler for the particle sampler.
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Note that a sampler must be updated **before** it can be used.
 */
DANTON_API int danton_sampler_update(struct danton_sampler * sampler);

/**
 * Create a simulation context.
 *
 * @return A handle for the context or `NULL` on failure.
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
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * Depending on the context *mode* a Monte-Carlo simulation or a grammage
 * scan is done. Note that setting *requested* to zero or less ignores this
 * option, resulting in all events to be processed.
 */
DANTON_API int danton_context_run(
    struct danton_context * context,
    long events,
    long requested);

/**
 * Get the current number of unprocessed errors.
 *
 * @param  context  The context of interest or `NULL`.
 * @return The corresponding error count.
 *
 * If no *context* is provided the number of global errors is returned.
 */
DANTON_API int danton_error_count(struct danton_context * context);

/**
 * Retrieve the last error in the stack.
 *
 * @param  context  The context of interest or `NULL`.
 * @return A string describing the last error.
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
 * @return `EXIT_SUCCESS` on success, `EXIT_FAILURE` otherwise.
 *
 * The format string and variable arguments follow the printf syntax.
 */
DANTON_API int danton_error_push(
    struct danton_context * context,
    const char * format,
    ...);


/* ============================================================================
 *
 * Rust interface.
 *
 * ============================================================================
 */

/* Reset the context on a physics change. */
DANTON_API void danton_context_reset(struct danton_context * context);

/* Set the temporary random context. */
struct danton_context_random;
DANTON_API void danton_context_random_set(
    struct danton_context * context, struct danton_context_random * random);

extern double danton_context_random_open01(
    struct danton_context_random * random);

/* Map the media materials. */
DANTON_API void danton_materials_set();

/* Expose the physics engines. */
struct danton_physics {
        void * ent;
        void * pumas;
};

extern struct danton_physics * danton_physics;

extern double * danton_tau_ctau0;
extern double * danton_tau_mass;

/* Expose the materials indices. */
struct danton_material_index {
        int rock;
        int water;
        int air;
        int topography;
};

extern struct danton_material_index * danton_material_index;

/* Tracer interface. */
struct danton_tracer;
DANTON_API struct danton_tracer * danton_tracer_create(void);
DANTON_API int danton_tracer_medium(
    struct danton_tracer * tracer, const double * position);
DANTON_API int danton_tracer_trace(
    struct danton_tracer * tracer,
    const double * position,
    const double * direction,
    double * distance,
    int * next_medium
);

#ifdef __cplusplus
}
#endif
#endif
