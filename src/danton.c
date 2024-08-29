/*
 * Copyright (C) 2017 Université Clermont Auvergne, CNRS/IN2P3, LPC
 * Author: Valentin NIESS (niess@in2p3.fr)
 *
 * This software is a C99 library for the simulation of the coupled transport
 * of ultra high energy taus and neutrinos through the Earth, by Monte-Carlo.
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

/* Standard library includes. */
#include <errno.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The various APIs. */
#include "alouette.h"
#include "danton.h"
#include "ent.h"
#include "pumas.h"
#include "turtle.h"

/* The spherical Earth radius, in m. */
#define PREM_EARTH_RADIUS 6371.E+03

/* The WGS 84 ellipsoid parameters, in m. */
#define WGS84_RADIUS_A 6378137.0
#define WGS84_RADIUS_B 6356752.314245

/* The radius of the geostationary orbit, in m. */
#define GEO_ORBIT 42164E+03

/* Avogadro's number. */
#define PHYS_NA 6.022E+23

/* Biasing factor for backward tau decays. */
#define DECAY_BIAS 6.

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Handles for Physic. */
static struct {
        struct ent_physics * ent;
        struct pumas_physics * pumas;
} physics = { NULL, NULL };

struct danton_physics * danton_physics = (void *)&physics;

/* Path for data sets. */
static char * geoid_path = NULL;

/* The tau lepton mass, in GeV / c^2. */
static double tau_mass = 0.0;
double * danton_tau_mass = &tau_mass;

/* The tau lifetime, in m / c. */
static double tau_ctau0 = 0.0;
double * danton_tau_ctau0 = &tau_ctau0;

/* Callbacks for locking and unlocking critical sections. */
static danton_lock_cb * lock = NULL;
static danton_lock_cb * unlock = NULL;

/* Callback for PUMAS or TURTLE errors */
#define ERROR_BUFFER_SIZE 2048
static char last_error[ERROR_BUFFER_SIZE] = { 0x0 };
static pumas_function_t * last_caller = NULL;

static void catch_pumas_error(enum pumas_return rc, pumas_function_t * caller,
    const char * message)
{
        last_caller = caller;
        strncpy(last_error, message, ERROR_BUFFER_SIZE - 1);
}

static void catch_turtle_error(enum turtle_return rc,
    turtle_function_t * caller, const char * message)
{
        last_caller = (pumas_function_t *)caller;
        strncpy(last_error, message, ERROR_BUFFER_SIZE - 1);
}

/* Shortcut for dumping a PUMAS error. */
#define ERROR_PUMAS(context)                                                   \
        danton_error_push(context, "%s (%d): error in %s, `%s`.", __FILE__,    \
            __LINE__, pumas_error_function(last_caller),                       \
            last_error)

/* Shortcut for dumping a TURTLE error. */
#define ERROR_TURTLE(context)                                                  \
        danton_error_push(context, "%s (%d): error in %s, `%s`.", __FILE__,    \
            __LINE__, turtle_error_function((turtle_function_t *)last_caller), \
            last_error)

/* Density according to the Preliminary Earth Model (PEM). */
static double pem_model0(double x, double * density)
{
        const double a2 = -8.8381E+03;
        *density = 13.0885E+03 + a2 * x * x;
        const double xg = (x <= 5E-02) ? 5E-02 : x;
        return PREM_EARTH_RADIUS / fabs(2. * a2 * xg);
}

static double pem_model1(double x, double * density)
{
        const double a = 1.2638E+03;
        *density = 12.58155E+03 + x * (-a + x * (-3.6426E+03 - x * 5.5281E+03));
        return PREM_EARTH_RADIUS / a;
}

static double pem_model2(double x, double * density)
{
        const double a = 6.4761E+03;
        *density = 7.9565E+03 + x * (-a + x * (5.5283E+03 - x * 3.0807E+03));
        return PREM_EARTH_RADIUS / a;
}

static double pem_model3(double x, double * density)
{
        const double a = 1.4836E+03;
        *density = 5.3197E+03 - a * x;
        return PREM_EARTH_RADIUS / a;
}

static double pem_model4(double x, double * density)
{
        const double a = 8.0298E+03;
        *density = 11.2494E+03 - a * x;
        return 0.01 * PREM_EARTH_RADIUS / a;
}

static double pem_model5(double x, double * density)
{
        const double a = 3.8045E+03;
        *density = 7.1089E+03 - a * x;
        return PREM_EARTH_RADIUS / a;
}

static double pem_model6(double x, double * density)
{
        const double a = 0.6924E+03;
        *density = 2.691E+03 + a * x;
        return PREM_EARTH_RADIUS / a;
}

static double pem_model7(double x, double * density)
{
        *density = 2.9E+03;
        return 0.;
}

static double pem_model8(double x, double * density)
{
        *density = 2.6E+03;
        return 0.;
}

static double pem_model9(double x, double * density)
{
        *density = 1.02E+03;
        return 0.;
}

/* The U.S. standard atmosphere model. */
#define USS_MODEL(INDEX, B, C)                                                 \
        static double uss_model##INDEX(double x, double * density)             \
        {                                                                      \
                *density = B / C * exp(-(x - 1.) * PREM_EARTH_RADIUS / C);     \
                return C;                                                      \
        }

USS_MODEL(0, 12226.562, 9941.8638)
USS_MODEL(1, 11449.069, 8781.5355)
USS_MODEL(2, 13055.948, 6361.4304)
USS_MODEL(3, 5401.778, 7721.7016)

/* Outer space density model. */
static double space_model0(double x, double * density)
{
        *density = 1.E-21; /* ~10^6 H per m^-3. */
        return 0.;
}

/* Structure for storing the data relative to the sampled events. */
struct event_record {
        struct danton_event api;
        struct danton_state primary;
        struct danton_state vertex;
        struct danton_state final;
        int buffer_size;
        struct danton_product product[];
};

/* Container for error message(s). */
struct error_stack {
        int count;
        int size;
#define ERROR_SIZE 512
        char data[ERROR_SIZE];
};

/* Container for contextual simulation data. */
struct simulation_context {
        /* Public API data. */
        struct danton_context api;

        /* Sub contexts for the transport engines. */
        struct pumas_context * pumas;
        struct ent_context ent;

        struct turtle_client * client;

        /* Recorded events. */
        struct event_record * record;

        /* The lower (upper) energy bound under (above) which all particles are
         * killed.
         */
        double energy_cut;

        /* Flag to check if a flux sampling is requested. */
        int sample_flux;

        /* Flag to check if the event was published. */
        int published;

        /* Temporary random context. */
        struct danton_context_random * random;

        struct error_stack error;
};

DANTON_API void danton_context_random_set(
    struct danton_context * context, struct danton_context_random * random)
{
        struct simulation_context * context_ = (void *)context;
        context_->random = random;
}

static double random_uniform01(struct simulation_context * context)
{
        return danton_context_random_open01(context->random);
}

/* Get the full simulation context from the PUMAS' one. */
static struct simulation_context * pumas2context(struct pumas_context * context)
{
        return (struct simulation_context *)context->user_data;
}

/* Get the full simulation context from the ENT's one. */
static struct simulation_context * ent2context(struct ent_context * context)
{
        char * p = (char *)context;
        p -= offsetof(struct simulation_context, ent);
        return (struct simulation_context *)p;
}

/* Generic Monte-Carlo state with stepping data. */
struct generic_state {
        union {
                struct ent_state ent;
                struct pumas_state pumas;
        } base;
        struct simulation_context * context;
        int medium;
        double density;
        double x;
        int is_tau;
        int is_inside;
        int has_crossed;
        int flux_direction;
};

/* Supported reference systems for the Earth model */
enum earth_reference { EARTH_PREM = 0, EARTH_WGS84, EARTH_EGM96 };

/* Global settings for the Earth model. */
static struct {
        enum earth_reference reference;
        struct turtle_stack * stack;
        struct turtle_map * undulations;
        double z0;
        int flat_topography;
        double density;
        int sea;
} earth = { EARTH_PREM, NULL, NULL, 0., 0, 2.65E+03, 1 };

/* API function for accessing the turtle topography stack. */
void * danton_get_topography(void) { return earth.stack; }

/* ENT density callback for the topography. */
static double density_topography(
    struct ent_medium * medium, struct ent_state * state, double * density)
{
        struct generic_state * s = (struct generic_state *)state;
        *density = s->density = earth.density;
        return 0.;
}

/* PUMAS locals callbacks for the topography. */
static double locals_topography(struct pumas_medium * medium,
    struct pumas_state * state, struct pumas_locals * locals)
{
        memset(locals->magnet, 0x0, sizeof(locals->magnet));
        struct generic_state * s = (struct generic_state *)state;
        s->density = locals->density = earth.density;
        return 0.;
}

/* Density callbacks for ENT. */
#define DENSITY(MODEL, INDEX)                                                  \
        static double density_##MODEL##INDEX(struct ent_medium * medium,       \
            struct ent_state * state, double * density)                        \
        {                                                                      \
                struct generic_state * s = (struct generic_state *)state;      \
                const double step = MODEL##_model##INDEX(s->x, density);       \
                s->density = *density;                                         \
                return 0.01 * step;                                            \
        }

DENSITY(pem, 0)
DENSITY(pem, 1)
DENSITY(pem, 2)
DENSITY(pem, 3)
DENSITY(pem, 4)
DENSITY(pem, 5)
DENSITY(pem, 6)
DENSITY(pem, 7)
DENSITY(pem, 8)
DENSITY(pem, 9)
DENSITY(uss, 0)
DENSITY(uss, 1)
DENSITY(uss, 2)
DENSITY(uss, 3)
DENSITY(space, 0)

/* Local callbacks for PUMAS. */
#define LOCALS(MODEL, INDEX)                                                   \
        static double locals_##MODEL##INDEX(struct pumas_medium * medium,      \
            struct pumas_state * state, struct pumas_locals * locals)          \
        {                                                                      \
                struct generic_state * s = (struct generic_state *)state;      \
                const double step =                                            \
                    MODEL##_model##INDEX(s->x, &locals->density);              \
                memset(locals->magnet, 0x0, sizeof(locals->magnet));           \
                s->density = locals->density;                                  \
                return step;                                                   \
        }

LOCALS(pem, 0)
LOCALS(pem, 1)
LOCALS(pem, 2)
LOCALS(pem, 3)
LOCALS(pem, 4)
LOCALS(pem, 5)
LOCALS(pem, 6)
LOCALS(pem, 7)
LOCALS(pem, 8)
LOCALS(pem, 9)
LOCALS(uss, 0)
LOCALS(uss, 1)
LOCALS(uss, 2)
LOCALS(uss, 3)
LOCALS(space, 0)

/* Helper function for computing geodetic coordinates from ECEF ones. */
static double compute_geodetic(
    double x, const double * position, double * latitude, double * longitude)
{
        if (earth.reference == EARTH_PREM)
                return (x - 1.) * PREM_EARTH_RADIUS;

        double altitude, trash;
        if (latitude == NULL) latitude = &trash;
        if (longitude == NULL) longitude = &trash;
        turtle_ecef_to_geodetic(
            (double *)position, latitude, longitude, &altitude);
        return altitude;
}

/* Helper function for computing ECEF coordinates. */
static void compute_ecef_position(
    double latitude, double longitude, double altitude, double * ecef)
{
        if (earth.reference == EARTH_PREM) {
                const double deg = M_PI / 180.;
                const double theta = (90. - latitude) * deg;
                const double phi = longitude * deg;
                const double s = sin(theta);
                const double r = PREM_EARTH_RADIUS + altitude;
                ecef[0] = r * s * cos(phi);
                ecef[1] = r * s * sin(phi);
                ecef[2] = r * cos(theta);
        } else {
                turtle_ecef_from_geodetic(
                    latitude, longitude, altitude, ecef);
        }
}

/* Helper function for computing ECEF direction from horizontal angular
 * coordinates.
 */
static void compute_ecef_direction(
    double latitude, double longitude, double azimuth, double c, double * ecef)
{
        if (earth.reference == EARTH_PREM) {
                /* Compute the rotation matrix from local to ECEF. */
                const double deg = M_PI / 180.;
                const double theta = (90. - latitude) * deg;
                const double phi = longitude * deg;
                const double st = sin(theta);
                const double ct = cos(theta);
                const double sp = sin(phi);
                const double cp = cos(phi);

                const double R[3][3] = { { -sp, cp, 0. },
                        { -ct * cp, -ct * sp, st }, { st * cp, st * sp, ct } };

                /* Compute the local direction. */
                const double p = (90. - azimuth) * deg;
                const double s = sqrt(1. - c * c);
                const double u[3] = { s * cos(p), s * sin(p), c };

                /* Apply the rotation. */
                ecef[0] = R[0][0] * u[0] + R[1][0] * u[1] + R[2][0] * u[2];
                ecef[1] = R[0][1] * u[0] + R[1][1] * u[1] + R[2][1] * u[2];
                ecef[2] = R[0][2] * u[0] + R[1][2] * u[1] + R[2][2] * u[2];
        } else {
                const double elevation = 90. - acos(c) / M_PI * 180.;
                turtle_ecef_from_horizontal(
                    latitude, longitude, azimuth, elevation, ecef);
        }
}

/* Get the parameters for computing the intersection with the ellipsoid. */
static void ellipsoid_parameters_intersection(const double * position,
    const double * direction, double * a, double * b, double * r2)
{
        double ai, bi;
        if (earth.reference == EARTH_PREM) {
                ai = bi = 1. / PREM_EARTH_RADIUS;
        } else {
                ai = 1. / WGS84_RADIUS_A;
                bi = 1. / WGS84_RADIUS_B;
        }

        const double rx = ai * position[0];
        const double ry = ai * position[1];
        const double rz = bi * position[2];
        *r2 = rx * rx + ry * ry + rz * rz;

        const double ux = ai * direction[0];
        const double uy = ai * direction[1];
        const double uz = bi * direction[2];
        *a = ux * ux + uy * uy + uz * uz;
        *b = rx * ux + ry * uy + rz * uz;
}

/* Generic medium callback. */
static double medium(const double * position, const double * direction,
    struct generic_state * state)
{
        state->medium = -1;
        double step = 0.;

        double a, b, r2;
        ellipsoid_parameters_intersection(position, direction, &a, &b, &r2);
        const double rmax = GEO_ORBIT / PREM_EARTH_RADIUS;
        if (r2 > rmax * rmax) return step;
        const double r = sqrt(r2);
        state->x = r;

        double latitude, longitude, altitude = -DBL_MAX;
        if (!state->context->api.decay && (state->has_crossed >= 0)) {
                /* Check the flux boundary in forward MC. */
                altitude = compute_geodetic(r, position, &latitude, &longitude);
                const double zi = state->context->api.sampler->altitude[0];
                if (state->is_inside < 0)
                        state->is_inside = (altitude < zi) ? 1 : 0;
                else if (((state->flux_direction >= 0) && state->is_inside &&
                    (altitude >= zi)) || ((state->flux_direction <= 0) &&
                    (!state->is_inside) && (altitude <= zi))) {
                        state->has_crossed = 1;
                        return step;
                }
        };

        const double ri[] = { 1221.5E+03 / PREM_EARTH_RADIUS,
                3480.E+03 / PREM_EARTH_RADIUS, 5701.E+03 / PREM_EARTH_RADIUS,
                5771.E+03 / PREM_EARTH_RADIUS, 5971.E+03 / PREM_EARTH_RADIUS,
                6151.E+03 / PREM_EARTH_RADIUS, 6346.6E+03 / PREM_EARTH_RADIUS,
                6356.E+03 / PREM_EARTH_RADIUS, 6368.E+03 / PREM_EARTH_RADIUS,
                1., 1. + 4.E+03 / PREM_EARTH_RADIUS,
                1. + 1.E+04 / PREM_EARTH_RADIUS,
                1. + 4.E+04 / PREM_EARTH_RADIUS,
                1. + 1.E+05 / PREM_EARTH_RADIUS,
                GEO_ORBIT / PREM_EARTH_RADIUS };

        /* Kill neutrinos that exit the atmosphere.  */
        if ((!state->is_tau) && (state->x > ri[13])) return step;

        int i;
        for (i = 0; i < sizeof(ri) / sizeof(*ri) - 1; i++) {
                if (r <= ri[i]) {
                        state->medium = i;

                        /* Outgoing intersection. */
                        const double d2 = b * b + a * (ri[i] * ri[i] - r2);
                        const double d = (d2 <= 0.) ? 0. : sqrt(d2);
                        step = (d - b) / a;

                        if ((i > 0) && (b < 0.)) {
                                /* This is a downgoing trajectory. Let
                                 * us compute the intersection with the lower
                                 * radius.
                                 */
                                const double r1 = ri[i - 1];
                                const double d2 = b * b + a * (r1 * r1 - r * r);
                                if (d2 > 0.) {
                                        const double d = sqrt(d2);
                                        double s = -(b + d) / a;
                                        if ((s > 0.) && (s < step)) step = s;
                                }
                        }
#define STEP_MIN 1E-03
                        if (step < STEP_MIN) step = STEP_MIN;
                        break;
                }
        }

        /* Check for any topography. */
#define MEDIUM_TOPOGRAPHY 100
        if ((i < 8) || (i > 11) ||
            ((earth.undulations == NULL) && !earth.flat_topography))
                return step;

        if (altitude == -DBL_MAX)
                altitude = compute_geodetic(r, position, &latitude, &longitude);

        double z_sea = danton_geoid_undulation(latitude, longitude);

        /* Let us compute the ground altitude. */
        double zg = -DBL_MAX;
        if (earth.flat_topography)
                zg = earth.z0;
        else if (earth.stack != NULL) {
                int inside;
                if (lock != NULL) {
                        turtle_client_elevation(state->context->client,
                            latitude, longitude, &zg, &inside);
                } else {
                        turtle_stack_elevation(
                            earth.stack, latitude, longitude, &zg, &inside);
                }
                if (!inside) zg = 0.;
        }
        zg += z_sea;

        /* Let us update the medium index accordingly. */
        double s = altitude - zg;
        if (s <= 0.) {
                if (i >= 9) state->medium = MEDIUM_TOPOGRAPHY;
        } else if (earth.sea && (altitude < z_sea)) {
                state->medium = 9;
        } else {
                if (i < 10) state->medium = 10;
        }

        /* Finally let us update the step length. */
        s = 0.4 * fabs(s);
        if (s < step) step = s;
        if (step < STEP_MIN) step = STEP_MIN;
        return step;

#undef STEP_MIN
}

/* Media table for ENT. */

static struct ent_medium media_ent[] = { { 0., 0., &density_pem0 },
        { 0., 0., &density_pem1 }, { 0., 0., &density_pem2 },
        { 0., 0., &density_pem3 }, { 0., 0., &density_pem4 },
        { 0., 0., &density_pem5 }, { 0., 0., &density_pem6 },
        { 0., 0., &density_pem7 }, { 0., 0., &density_pem8 },
        { 0., 0., &density_pem9 }, { 0., 0., &density_uss0 },
        { 0., 0., &density_uss1 }, { 0., 0., &density_uss2 },
        { 0., 0., &density_uss3 }, { 0., 0., &density_space0 } };

static struct ent_medium topography_ent = { 0., 0., &density_topography };

/* Medium callback encapsulation for ENT. */
static double medium_ent(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium_ptr)
{
        double direction[3];
        if (context->ancestor) {
                direction[0] = -state->direction[0];
                direction[1] = -state->direction[1];
                direction[2] = -state->direction[2];
        } else {
                direction[0] = state->direction[0];
                direction[1] = state->direction[1];
                direction[2] = state->direction[2];
        }
        struct generic_state * g = (struct generic_state *)state;
        const double step = medium(state->position, direction, g);
        if (g->medium == MEDIUM_TOPOGRAPHY)
                *medium_ptr = &topography_ent;
        else if (g->medium >= 0)
                *medium_ptr = media_ent + g->medium;
        else
                *medium_ptr = NULL;
        return step;
}

/* Media indices for PUMAS (mapped from Rust). */
static struct danton_material_index material_index = { 0 };
struct danton_material_index * danton_material_index = &material_index;

/* Media table for PUMAS. */
static struct pumas_medium media_pumas[] = { { 0, &locals_pem0 },
        { 0, &locals_pem1 }, { 0, &locals_pem2 }, { 0, &locals_pem3 },
        { 0, &locals_pem4 }, { 0, &locals_pem5 }, { 0, &locals_pem6 },
        { 0, &locals_pem7 }, { 0, &locals_pem8 }, { 0, &locals_pem9 },
        { 0, &locals_uss0 }, { 0, &locals_uss1 }, { 0, &locals_uss2 },
        { 0, &locals_uss3 }, { 0, &locals_space0 } };

static struct pumas_medium topography_pumas = { 0, &locals_topography };

/* Convert material data from PUMAS to ENT */
static void convert_material(int material, double * Z_ptr, double * A_ptr)
{
        /* Get atomic elements from PUMAS */
        int n;
        pumas_physics_material_properties(physics.pumas, material, &n,
            NULL, NULL, NULL, NULL);
        int components[n];
        double fractions[n];
        pumas_physics_material_properties(physics.pumas, material, NULL,
            NULL, NULL, components, fractions);

        double Z = 0., A = 0., w = 0.;
        int i;
        for (i = 0; i < n; i++) {
                double Zi, Ai;
                pumas_physics_element_properties(
                    physics.pumas, components[i], &Zi, &Ai, NULL);
                const double wi = fractions[i] / Ai;
                w += wi;
                Z += Zi * wi;
                A += Ai * wi;
        }

        *Z_ptr = Z / w;
        *A_ptr = A / w;
}

/* Initialise the material data */
DANTON_API void danton_materials_set()
{
        /* Set PUMAS material indices */
#define N_LAYERS 15
        media_pumas[ 0].material = material_index.rock;
        media_pumas[ 1].material = material_index.rock;
        media_pumas[ 2].material = material_index.rock;
        media_pumas[ 3].material = material_index.rock;
        media_pumas[ 4].material = material_index.rock;
        media_pumas[ 5].material = material_index.rock;
        media_pumas[ 6].material = material_index.rock;
        media_pumas[ 8].material = material_index.rock;
        media_pumas[ 9].material = material_index.water;
        media_pumas[10].material = material_index.air;
        media_pumas[11].material = material_index.air;
        media_pumas[12].material = material_index.air;
        media_pumas[13].material = material_index.air;
        media_pumas[14].material = material_index.air;

        topography_pumas.material = material_index.topography;

        /* Set the outer PREM layer to rock or water */
        if ((earth.stack != NULL) || earth.flat_topography ||
            (earth.sea == 0)) {
                media_pumas[9].material = material_index.rock;
                media_pumas[9].locals = &locals_pem8;
        } else {
                media_pumas[9].material = material_index.water;
                media_pumas[9].locals = &locals_pem9;
        }

        /* Set material properties for ENT */
        int i;
        struct ent_medium * m;
        for (i = 0, m = media_ent; i < N_LAYERS; i++, m++) {
                /* Get atomic elements from PUMAS */
                int material = media_pumas[i].material;
                convert_material(material, &m->Z, &m->A);
        }
        convert_material(
            material_index.topography, &topography_ent.Z, &topography_ent.A);
#undef N_LAYERS
}

/* Medium callback encapsulation for PUMAS. */
static enum pumas_step medium_pumas(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr,
    double * step_ptr)
{
        double direction[3];
        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                direction[0] = state->direction[0];
                direction[1] = state->direction[1];
                direction[2] = state->direction[2];
        } else {
                direction[0] = -state->direction[0];
                direction[1] = -state->direction[1];
                direction[2] = -state->direction[2];
        }
        struct generic_state * g = (struct generic_state *)state;
        const double step = medium(state->position, direction, g);

        if (medium_ptr != NULL) {
                if (g->medium == MEDIUM_TOPOGRAPHY)
                        *medium_ptr = &topography_pumas;
                else if (g->medium >= 0)
                        *medium_ptr = media_pumas + g->medium;
                else
                        *medium_ptr = NULL;
        }

        if (step_ptr != NULL) {
                *step_ptr = step;
        }

        return PUMAS_STEP_CHECK;
}

/* Encapsulation of the random engine for PUMAS. */
static double random_pumas(struct pumas_context * context)
{
        struct simulation_context * c = pumas2context(context);
        return random_uniform01(c);
}

/* Encapsulation of the random engine for ENT. */
static double random_ent(struct ent_context * context)
{
        struct simulation_context * c = ent2context(context);
        return random_uniform01(c);
}

/* Encapsulation of the random engine for ALOUETTE. */
static struct simulation_context * alouette_current_context = NULL;

static float random_alouette(void)
{
        return random_uniform01(alouette_current_context);
}

/* Set a neutrino state from a tau decay product. */
static void copy_neutrino(struct simulation_context * context,
    const struct pumas_state * tau, int pid, const double * momentum,
    struct ent_state * neutrino, double * direction)
{
        neutrino->pid = pid;
        neutrino->energy = sqrt(momentum[0] * momentum[0] +
            momentum[1] * momentum[1] + momentum[2] * momentum[2]);
        memcpy(neutrino->position, tau->position, sizeof(neutrino->position));
        if (context->api.longitudinal) {
                neutrino->direction[0] = direction[0];
                neutrino->direction[1] = direction[1];
                neutrino->direction[2] = direction[2];
        } else {
                neutrino->direction[0] = momentum[0] / neutrino->energy;
                neutrino->direction[1] = momentum[1] / neutrino->energy;
                neutrino->direction[2] = momentum[2] / neutrino->energy;
        }
        neutrino->distance = tau->distance;
        neutrino->grammage = tau->grammage;
        neutrino->weight = tau->weight;
}

/* Copy an ENT state to the event record. */
static void record_copy_ent(struct danton_state * dst, struct ent_state * src)
{
        if (dst == NULL) return;
        dst->pid = src->pid;
        dst->energy = src->energy;
        memcpy(dst->position, src->position, sizeof(dst->position));
        memcpy(dst->direction, src->direction, sizeof(dst->direction));
}

/* Copy a PUMAS state to the event record. */
static void record_copy_pumas(
    struct danton_state * dst, struct pumas_state * src)
{
        if (dst == NULL) return;
        dst->pid = (src->charge > 0.) ? ENT_PID_TAU_BAR : ENT_PID_TAU;
        dst->energy = src->energy + tau_mass;
        memcpy(dst->position, src->position, sizeof(dst->position));
        memcpy(dst->direction, src->direction, sizeof(dst->direction));
}

/* Copy an ALOUETTE decay product to the event record. */
static void record_copy_product(
    struct simulation_context * context, int pid, double * momentum)
{
        /* Manage the memory. */
        struct event_record * record = context->record;
        if (record->api.n_products == record->buffer_size) {
                const int n = 2 * record->buffer_size;
                record = realloc(
                    record, sizeof(*record) + n * sizeof(*record->product));
                if (record == NULL) return;
                context->record = record;
                record->buffer_size = n;
        }

        /* Append the decay product to the stack. */
        struct danton_product * product =
            record->product + record->api.n_products;
        product->pid = pid;
        memcpy(product->momentum, momentum, sizeof(product->momentum));
        record->api.n_products++;
}

/* Publish the event record to the recorder. */
static int record_publish(struct simulation_context * context)
{
        /* Check and prune the products. */
        struct event_record * record = context->record;
        if (context->api.decay && (record->api.n_products == 0))
                return EXIT_SUCCESS;
        if (record->api.n_products > 0)
                record->api.product = record->product;
        else
                record->api.product = NULL;

        /* Check the sampler in forward mode */
        if (context->api.mode == DANTON_MODE_FORWARD) {
                struct danton_sampler * sampler = context->api.sampler;
                const int pid = record->final.pid;
                int selected = 0;
                int j;
                for (j = 0; j < DANTON_PARTICLE_N; j++) {
                        if (sampler->weight[j] > 0.) {
                                if (pid == danton_particle_pdg(j)) {
                                        const double energy =
                                            record->final.energy;
                                        if ((energy >= sampler->energy[0]) &&
                                            (energy <= sampler->energy[1])) {
                                                selected = 1;
                                                break;
                                        }
                                }
                        }
                }
                if (!selected) {
                        return EXIT_SUCCESS;
                }
        }

        /* Check and prune the vertex */
        struct danton_state * vertex = record->api.vertex;
        if (record->vertex.pid != record->final.pid) {
                record->api.vertex = NULL;
        }

        /* Call the event processor. */
        int rc = context->api.recorder->record_event(
            &context->api, context->api.recorder, &record->api);

        /* Restore the vertex */
        record->api.vertex = vertex;

        /* Reset the record for new data. */
        record->api.n_products = 0;

        /* Flag the event as published. */
        context->published = 1;

        return rc;
}

static void stepping_pumas_flux(struct pumas_context * context,
    struct pumas_state * state_, struct pumas_medium * medium,
    enum pumas_event event)
{
        struct generic_state * state = (void *)state_;
        if (state->context->api.decay || (state->has_crossed < 0)) {
                return;
        }

        /* Update the last location in forward MC, for the flux crossing. */
        const double * const position = state_->position;
        const double * const direction = state_->direction;

        double a, b, r2;
        ellipsoid_parameters_intersection(position, direction, &a, &b, &r2);
        const double r = sqrt(r2);

        double altitude = -DBL_MAX;
        altitude = compute_geodetic(r, position, NULL, NULL);
        const double zi = state->context->api.sampler->altitude[0];
        state->is_inside = (altitude < zi) ? 1 : 0;
}

static void stepping_pumas(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium * medium_,
    enum pumas_event event)
{
        if (context->recorder->user_data != EXIT_SUCCESS)
                return;

        if (context->mode.direction == PUMAS_MODE_FORWARD) {
                /* Forward MC case */
                stepping_pumas_flux(context, state, medium_, event);
        }

        struct danton_state s;
        record_copy_pumas(&s, state);

        int medium;
        if (medium_ == NULL)
                medium = -1;
        else if (medium_ == &topography_pumas)
                medium = MEDIUM_TOPOGRAPHY;
        else
                medium = (int)(medium_ - media_pumas);

        struct generic_state * g = (void *)state;
        struct danton_context * c = (void *)g->context;
        struct danton_run_action * a = c->run_action;
        const int rc = a->call(c, a, DANTON_RUN_EVENT_STEP, medium, &s);

        if (rc != EXIT_SUCCESS) {
                context->recorder->user_data = (void *)EXIT_FAILURE;
        }
}

/* Encapsulation of pumas' calls with run action(s). */
static int call_pumas(
    struct simulation_context * context_, struct generic_state * state)
{
        struct danton_context * context = (struct danton_context *)context_;

        if (context->run_action != NULL) {
                context_->pumas->recorder->user_data = (void *)EXIT_SUCCESS;
        }

        /* Call PUMAS. */
        if (pumas_context_transport(
            context_->pumas, &state->base.pumas, NULL, NULL) !=
            PUMAS_RETURN_SUCCESS) {
                ERROR_PUMAS(context);
                return EXIT_FAILURE;
        }

        if (context->run_action != NULL) {
                if (context_->pumas->recorder->user_data != EXIT_SUCCESS)
                        return EXIT_FAILURE;
        }

        return EXIT_SUCCESS;
}

/* Custom stepping action for ENT. */
static enum ent_return stepping_ent_flux(struct ent_context * context,
    struct ent_medium * medium, struct ent_state * state_)
{
        struct generic_state * state = (void *)state_;
        if (state->context->api.decay || (state->has_crossed < 0)) {
                return ENT_RETURN_SUCCESS;
        }

        /* Update the last location in forward MC, for the flux crossing. */
        const double * const position = state_->position;
        const double * const direction = state_->direction;

        double a, b, r2;
        ellipsoid_parameters_intersection(position, direction, &a, &b, &r2);
        const double r = sqrt(r2);

        double altitude = -DBL_MAX;
        altitude = compute_geodetic(r, position, NULL, NULL);
        const double zi = state->context->api.sampler->altitude[0];
        state->is_inside = (altitude < zi) ? 1 : 0;

        return ENT_RETURN_SUCCESS;
}

static enum ent_return stepping_ent(struct ent_context * context,
    struct ent_medium * medium, struct ent_state * state)
{
        if (context->ancestor == NULL) {
                /* Forward MC case */
                stepping_ent_flux(context, medium, state);
        }

        struct danton_state s;
        record_copy_ent(&s, state);
        int m;
        if (medium == NULL)
                m = -1;
        else if (medium == &topography_ent)
                m = MEDIUM_TOPOGRAPHY;
        else
                m = (int)(medium - media_ent);

        struct generic_state * g = (void *)state;
        struct danton_context * c = (void *)g->context;
        struct danton_run_action * a = c->run_action;
        int rc = a->call(c, a, DANTON_RUN_EVENT_STEP, m, &s);

        if (rc == EXIT_SUCCESS)
                return ENT_RETURN_SUCCESS;
        else
                return ENT_RETURN_DOMAIN_ERROR;
}

/* Shortcut for dumping an ENT error. */
#define ERROR_ENT(context, rc, function)                                       \
        danton_error_push(context, "%s (%d): error in %s, `%s`.", __FILE__,    \
            __LINE__, ent_error_function((ent_function_t *)function),          \
            ent_error_string(rc))

/* Forward transport routine, recursive. */
static int transport_forward(struct simulation_context * context,
    struct ent_state * neutrino, int generation)
{
        if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
            (abs(neutrino->pid) != ENT_PID_NU_TAU)) {
                danton_error_push(&context->api, "%s (%s): invalid pid (%d).",
                    __FILE__, __LINE__, neutrino->pid);
                return EXIT_FAILURE;
        }

        /* Backup the initial direction if the transverse transport is
         * disabled.
         */
        double direction[3];
        if (context->api.longitudinal)
                memcpy(direction, neutrino->direction, sizeof(direction));

        int flux_direction;
        {
                struct generic_state * g_state =
                    (struct generic_state *)neutrino;
                flux_direction = g_state->flux_direction;
        }

        struct ent_state product;
        enum ent_event event;
        for (;;) {
                /* Neutrino transport with ENT. */
                enum ent_return rc;
                if ((rc = ent_transport(physics.ent, &context->ent, neutrino,
                         &product, &event)) != ENT_RETURN_SUCCESS) {
                        ERROR_ENT(&context->api, rc, ent_transport);
                        return EXIT_FAILURE;
                }
                if (neutrino->energy <= context->energy_cut + FLT_EPSILON)
                        break;
                if (context->sample_flux && (event == ENT_EVENT_EXIT)) {
                        /* Check for a flux crossing condition. */
                        struct generic_state * g_state =
                            (struct generic_state *)neutrino;
                        if (g_state->has_crossed) {
                                context->record->api.generation = generation;
                                record_copy_ent(
                                    context->record->api.final, neutrino);
                                if (record_publish(context) != EXIT_SUCCESS)
                                        return EXIT_FAILURE;
                                break;
                        }
                }
                if (event == ENT_EVENT_EXIT) break;
                if (context->api.longitudinal) {
                        memcpy(neutrino->direction, direction,
                            sizeof(neutrino->direction));
                        memcpy(product.direction, direction,
                            sizeof(product.direction));
                }
                if (abs(neutrino->pid) == ENT_PID_TAU) {
                        /* Exchange the lepton and the product state. */
                        struct ent_state tmp;
                        memcpy(&tmp, neutrino, sizeof(tmp));
                        memcpy(neutrino, &product, sizeof(*neutrino));
                        memcpy(&product, &tmp, sizeof(product));
                }

                if (abs(product.pid) == ENT_PID_TAU) {
                        /* Tau transport with PUMAS. */
                        const double charge = (product.pid > 0) ? -1. : 1.;
                        const double kinetic = product.energy - tau_mass;
                        struct generic_state tau_data = {
                                .base.pumas = { charge, kinetic,
                                    product.distance, product.grammage, 0.,
                                    product.weight, .decayed = 0 },
                                .context = context,
                                .medium = -1,
                                .density = 0.,
                                .x = 0.,
                                .is_tau = 1,
                                .is_inside = -1,
                                .has_crossed =
                                    (context->sample_flux) ? 0 : -1,
                                .flux_direction = flux_direction
                        };
                        struct pumas_state * tau = &tau_data.base.pumas;
                        memcpy(&tau->position, &product.position,
                            sizeof(tau->position));
                        memcpy(&tau->direction, &product.direction,
                            sizeof(tau->direction));
                        context->record->api.vertex = &context->record->vertex;
                        record_copy_pumas(context->record->api.vertex, tau);
                        if (call_pumas(context, &tau_data) != EXIT_SUCCESS)
                                return EXIT_FAILURE;
                        if (tau->decayed) {
                                /* Tau decay with ALOUETTE/TAUOLA. */
                                const double p = sqrt(tau->energy *
                                    (tau->energy + 2. * tau_mass));
                                double momentum[3] = { p * tau->direction[0],
                                        p * tau->direction[1],
                                        p * tau->direction[2] };
                                const double sgn = (product.pid > 0) ? -1. : 1.;
                                double polarisation[3] = {
                                    sgn * tau->direction[0],
                                    sgn * tau->direction[1],
                                    sgn * tau->direction[2]};
                                int trials;
                                if (lock != NULL) lock();
                                alouette_current_context = context;
                                struct alouette_products products;
                                for (trials = 0; trials < 20; trials++) {
                                        if (alouette_decay(0, product.pid,
                                                momentum, polarisation,
                                                &products) ==
                                            ALOUETTE_RETURN_SUCCESS)
                                                break;
                                }
                                struct generic_state nu_e_data, nu_t_data;
                                memset(&nu_e_data, 0x0, sizeof(nu_e_data));
                                memset(&nu_t_data, 0x0, sizeof(nu_t_data));
                                struct ent_state *nu_e = NULL, *nu_t = NULL;
                                int j;
                                for (j = 0; j < products.size; j++) {
                                        const int pid = products.pid[j];
                                        momentum[0] = products.P[j][0];
                                        momentum[1] = products.P[j][1];
                                        momentum[2] = products.P[j][2];

                                        if (abs(pid) == 16) {
                                                /* Update the neutrino state
                                                 * with the nu_tau daughter.
                                                 */
                                                nu_t = &nu_t_data.base.ent;
                                                copy_neutrino(context, tau, pid,
                                                    momentum, nu_t, direction);
                                                continue;
                                        } else if (pid == -12) {
                                                nu_e = &nu_e_data.base.ent;
                                                copy_neutrino(context, tau, pid,
                                                    momentum, nu_e, direction);
                                                continue;
                                        } else if (!context->api.decay ||
                                            (pid == 12) || (abs(pid) == 13) ||
                                            (abs(pid) == 14))
                                                continue;

                                        /* Log the decay if in air. */
                                        if ((tau_data.medium < 10) ||
                                            (tau_data.medium ==
                                                MEDIUM_TOPOGRAPHY))
                                                continue;
                                        if (context->record->api.n_products ==
                                            0)
                                                record_copy_pumas(
                                                    context->record->api.final,
                                                    tau);
                                        record_copy_product(
                                            context, pid, momentum);
                                }
                                if (unlock != NULL) unlock();
                                if (context->record->api.n_products > 0) {
                                        context->record->api.generation =
                                            generation;
                                        if (record_publish(context) !=
                                            EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }
                                generation++;

                                /* Process any additional nu_e~ or nu_tau. */
                                if (nu_e != NULL) {
                                        nu_e_data.context = context;
                                        if (context->sample_flux) {
                                                nu_e_data.is_inside = -1;
                                                nu_e_data.has_crossed =
                                                    (context->sample_flux) ?
                                                    0 : -1,
                                                nu_e_data.flux_direction =
                                                    flux_direction;
                                        } else {
                                                nu_e_data.has_crossed = -1;
                                                nu_e_data.flux_direction = 0;
                                        }
                                        if (transport_forward(context, nu_e,
                                                generation) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }
                                if (nu_t != NULL) {
                                        nu_t_data.context = context;
                                        if (context->sample_flux) {
                                                nu_t_data.is_inside = -1;
                                                nu_t_data.has_crossed =
                                                    (context->sample_flux) ?
                                                    0 : -1,
                                                nu_t_data.flux_direction =
                                                    flux_direction;
                                        } else {
                                                nu_t_data.has_crossed = -1;
                                                nu_t_data.flux_direction = 0;
                                        }
                                        if (transport_forward(context, nu_t,
                                                generation) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }
                        } else if (tau_data.has_crossed == 1) {
                                record_copy_pumas(
                                    context->record->api.final, tau);
                                context->record->api.generation = generation;
                                if (record_publish(context) != EXIT_SUCCESS)
                                        return EXIT_FAILURE;
                        }
                }
                if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
                    (abs(neutrino->pid) != ENT_PID_NU_TAU))
                        break;
        }

        return EXIT_SUCCESS;
}

/* Ancestor callback for taus. */
static double ancestor_tau(
    struct ent_context * context, struct ent_state * state)
{
        struct generic_state * g = (struct generic_state *)state;
        return 1.63E-17 * pow(state->energy, 1.363) * g->density;
}

/* Main ancestor callback, for ENT backward sampling of interaction
 * processes.
 */
static double ancestor_cb(struct ent_context * context, enum ent_pid ancestor,
    struct ent_state * daughter)
{
        if (daughter->pid == ENT_PID_NU_BAR_E) {
                if (ancestor == ENT_PID_NU_BAR_E) return 1.;
                return 0.;
        } else if (daughter->pid == ENT_PID_NU_TAU) {
                if (ancestor == ENT_PID_NU_TAU)
                        return 1.;
                else if (ancestor == ENT_PID_TAU)
                        return ancestor_tau(context, daughter);
                return 0.;
        } else if (daughter->pid == ENT_PID_NU_BAR_TAU) {
                if (ancestor == ENT_PID_NU_BAR_TAU)
                        return 1.;
                else if (ancestor == ENT_PID_TAU_BAR)
                        return ancestor_tau(context, daughter);
                return 0.;
        } else if (daughter->pid == ENT_PID_TAU) {
                if (ancestor == ENT_PID_NU_TAU) return 1.;
                return 0.;
        } else if (daughter->pid == ENT_PID_TAU_BAR) {
                if (ancestor == ENT_PID_NU_BAR_TAU) return 1.;
                return 0.;
        } else
                return 0.;
}

/* Polarisation callback for ALOUETTE in backward mode. A 100% longitudinal
 * polarisation is assumed.
 */
static void polarisation_cb(
    int pid, const double momentum[3], double * polarisation)
{
        double nrm = momentum[0] * momentum[0] + momentum[1] * momentum[1] +
            momentum[2] * momentum[2];
        if (nrm <= 0.) {
                polarisation[0] = 0.;
                polarisation[1] = 0.;
                polarisation[2] = 0.;
                return;
        }
        nrm = 1. / sqrt(nrm);
        if (pid > 0) {
                nrm = -nrm;
        }
        polarisation[0] = momentum[0] * nrm;
        polarisation[1] = momentum[1] * nrm;
        polarisation[2] = momentum[2] * nrm;
}

/* Backward transport routine. */
static int transport_backward(
    struct simulation_context * context, struct generic_state * current)
{
        /* Backup the final state, e.g. the tau at decay. */
        if (context->record->api.generation == 1) {
                if (current->is_tau)
                        record_copy_pumas(
                            context->record->api.final, &current->base.pumas);
                else
                        record_copy_ent(
                            context->record->api.final, &current->base.ent);
        }

        struct generic_state g_state = {.context = current->context};
        struct ent_state * state = NULL;
        struct pumas_state * tau = NULL;
        double direction[3];
        if (current->is_tau) {
                /* Apply the BMC weight for the tau decay. */
                tau = &current->base.pumas;
                if (context->api.decay ||
                    (context->record->api.generation > 1)) {
                        const double Pf =
                            sqrt(tau->energy * (tau->energy + 2. * tau_mass));
                        tau->weight *= tau_mass / (tau_ctau0 * Pf);
                }

                /* Backward propagate the tau state. */
                memcpy(direction, tau->direction, sizeof(direction));
                const double lambda0 = 3E+07;
                const double p1 = 0.1;
                double x0;
                for (;;) {
                        x0 = tau->grammage;
                        context->pumas->limit.grammage =
                            x0 - lambda0 * log(random_uniform01(context));
                        tau->decayed = 0;
                        if (call_pumas(context, (struct generic_state *)tau) !=
                            EXIT_SUCCESS)
                                return EXIT_FAILURE;
                        if ((!tau->decayed &&
                                (tau->grammage <
                                    context->pumas->limit.grammage -
                                    FLT_EPSILON)) ||
                            (tau->energy + tau_mass >=
                                context->energy_cut - FLT_EPSILON) ||
                            (tau->weight <= 0.))
                                return EXIT_SUCCESS;
                        if (context->record->api.generation > 1) break;

                        /* Check that the tau is **not** emerging from the
                         * Earth.
                         */
                        double a, b, r2;
                        ellipsoid_parameters_intersection(
                            tau->position, tau->direction, &a, &b, &r2);
                        b = -b;
                        const double d2 = b * b + a * (1. - r2);
                        if ((d2 <= 0.) || (sqrt(d2) > -b)) break;

                        /* Check that the proposed vertex is **not** in air. */
                        struct generic_state * g = (struct generic_state *)tau;
                        if ((g->medium < 10) ||
                            (g->medium == MEDIUM_TOPOGRAPHY) ||
                            (g->density <= 0.))
                                break;

                        /* If upgoing and in air, randomly recycle the event
                         * by biasing the decay probability at the vertex.
                         * First let us compute the true decay probability.
                         */
                        const double Pf =
                            sqrt(tau->energy * (tau->energy + 2. * tau_mass));
                        const double lD = tau_ctau0 * Pf / tau_mass;
                        const double lB = lambda0 / g->density;
                        const double pD = lB / (lB + lD);
                        const double pB = lD / (lB + lD); /* For underflow. */
                        if ((pD <= 0.) || (pB <= 0.)) break;

                        /* Draw over the biased probability and re-weight
                         * accordingly.
                         */
                        if (random_uniform01(context) < p1) {
                                tau->weight *= pD / p1;
                                break;
                        } else
                                tau->weight *= pB / (1. - p1);
                }

                /* Let us record the tau state at production. */
                if (context->record->api.generation == 1) {
                        context->record->api.vertex = &context->record->vertex;
                        record_copy_pumas(
                            context->record->api.vertex, &current->base.pumas);
                }

                /* Backward generate the production vertex. */
                enum ent_pid pid =
                    (tau->charge < 0.) ? ENT_PID_TAU : ENT_PID_TAU_BAR;
                state = &g_state.base.ent;
                state->pid = pid;
                state->energy = tau->energy + tau_mass;
                state->distance = tau->distance;
                state->grammage = tau->grammage;
                state->weight = tau->weight;
                memcpy(state->position, tau->position, sizeof(state->position));
                memcpy(
                    state->direction, tau->direction, sizeof(state->direction));
                g_state.x = 0.;
                g_state.is_tau = 0;
                g_state.is_inside = -1;
                g_state.has_crossed = -1;
                g_state.flux_direction = 0;

                struct ent_medium * medium;
                medium_ent(&context->ent, state, &medium);
                if (medium == NULL) return EXIT_SUCCESS;
                enum ent_return re;
                if ((re = ent_vertex(physics.ent, &context->ent, state, medium,
                         ENT_PROCESS_NONE, NULL)) != ENT_RETURN_SUCCESS) {
                        ERROR_ENT(&context->api, re, ent_vertex);
                        return EXIT_FAILURE;
                }

                /* Append the effective BMC weight for the transport, in order
                 * to recover a flux convention.
                 */
                double cs;
                ent_physics_cross_section(physics.ent, pid, state->energy,
                    medium->Z, medium->A, ENT_PROCESS_NONE, &cs);
                double density;
                medium->density(medium, state, &density);
                const double lP = 1E-03 * medium->A / (cs * PHYS_NA * density);
                const double Pi =
                    sqrt(tau->energy * (tau->energy + 2. * tau_mass));
                const double lD = tau_ctau0 * Pi / tau_mass;
                const double lB = lambda0 / density;
                const double p0 = exp(-(tau->grammage - x0) / lambda0);
                state->weight *= lB * lD / ((lB + lD) * lP * p0);

                /* Reset the initial direction if the transverse transport is
                 * disabled. */
                if (context->api.longitudinal)
                        memcpy(state->direction, direction,
                            sizeof(state->direction));
        } else {
                state = &current->base.ent;
                tau = &g_state.base.pumas;
                memcpy(direction, state->direction, sizeof(direction));
        }

        /* Backward propagate the neutrino. */
        enum ent_event event = ENT_EVENT_NONE;
        while ((event != ENT_EVENT_EXIT) &&
            (state->energy < context->energy_cut - FLT_EPSILON)) {
                enum ent_return re;
                if ((re = ent_transport(physics.ent, &context->ent, state, NULL,
                         &event)) != ENT_RETURN_SUCCESS) {
                        ERROR_ENT(&context->api, re, ent_transport);
                        return EXIT_FAILURE;
                }
                if (state->weight <= 0.) return EXIT_SUCCESS;
                if (context->api.longitudinal)
                        memcpy(state->direction, direction,
                            sizeof(state->direction));

                if (event == ENT_EVENT_DECAY_TAU) {
                        /* Backward randomise the tau decay. */
                        double momentum[3] = { state->energy *
                                    state->direction[0],
                                state->energy * state->direction[1],
                                state->energy * state->direction[2] };
                        double weight;
                        int trials;
                        if (lock != NULL) lock();
                        alouette_current_context = context;
                        struct alouette_products products;
                        int mother;
                        double bias;
                        if (state->pid > 0) {
                                mother = ENT_PID_TAU;
                                bias = -1;
                        } else {
                                mother = ENT_PID_TAU_BAR;
                                bias = 1.;
                        }
                        for (trials = 0; trials < 20; trials++) {
                                if (alouette_undecay(0, state->pid, mother,
                                    momentum, &polarisation_cb, bias,
                                    &products)
                                    == ALOUETTE_RETURN_SUCCESS)
                                        break;
                        }
                        weight = products.weight;

                        const int pid1 = products.pid[0];
                        if ((products.size <= 0) ||
                            (abs(pid1) != ENT_PID_TAU)) {
                                if (unlock != NULL) unlock();
                                return EXIT_SUCCESS;
                        }
                        momentum[0] = products.P[0][0];
                        momentum[1] = products.P[0][1];
                        momentum[2] = products.P[0][2];
                        if (unlock != NULL) unlock();
                        const double p12 = momentum[0] * momentum[0] +
                            momentum[1] * momentum[1] +
                            momentum[2] * momentum[2];
                        const double E1 = sqrt(p12 + tau_mass * tau_mass);
                        if (E1 >= context->energy_cut - FLT_EPSILON)
                                return EXIT_SUCCESS;

                        /* Update the tau state and start again the BMC
                         * transport.
                         */
                        tau->charge = (pid1 > 0) ? -1. : 1.;
                        tau->energy = E1 - tau_mass;
                        tau->distance = state->distance;
                        tau->grammage = state->grammage;
                        tau->time = 0.;
                        tau->weight = state->weight * weight * state->energy *
                            state->energy / p12;
                        memcpy(tau->position, state->position,
                            sizeof(tau->position));
                        struct generic_state * g = (struct generic_state *)tau;
                        if (context->api.longitudinal) {
                                if (g != current)
                                        memcpy(tau->direction, direction,
                                            sizeof(tau->direction));
                        } else {
                                const double d = 1. / sqrt(p12);
                                tau->direction[0] = momentum[0] * d;
                                tau->direction[1] = momentum[1] * d;
                                tau->direction[2] = momentum[2] * d;
                        }
                        tau->decayed = 0;
                        g->x = 0.;
                        g->is_tau = 1;
                        g->is_inside = -1;
                        g->has_crossed = -1;
                        g->flux_direction = 0;
                        context->record->api.generation++;
                        return transport_backward(context, g);
                }
        }
        if (event != ENT_EVENT_EXIT) return EXIT_SUCCESS;

        /* Let us sample the primary flux. */
        enum danton_particle index = danton_particle_index(state->pid);
        struct danton_primary * primary = context->api.primary[index];
        if (primary == NULL) return EXIT_SUCCESS;
        const double flux = primary->flux(primary, state->energy);
        if (flux <= 0.) return EXIT_SUCCESS;
        state->weight *= flux;

        /* This is a valid event. Let us record the primary and set the
         * weight.
         */
        context->record->api.weight = state->weight;
        record_copy_ent(context->record->api.primary, state);

        /* In flux mode let us publish the record and then return. */
        if (!context->api.decay) {
                if (record_publish(context) != EXIT_SUCCESS)
                        return EXIT_FAILURE;
                return EXIT_SUCCESS;
        }

        /* In full mode let us perform the tau decay with ALOUETTE/TAUOLA. */
        enum ent_pid pid = context->record->final.pid;
        const double p = sqrt((context->record->final.energy + tau_mass) *
            (context->record->final.energy - tau_mass));
        double momentum[3] = { p * context->record->final.direction[0],
                p * context->record->final.direction[1],
                p * context->record->final.direction[2] };
        int trials;
        if (lock != NULL) lock();
        alouette_current_context = context;
        struct alouette_products products;
        const double sgn = (pid > 0) ? -1. : 1.;
        const double * const dir = context->record->final.direction;
        double polarisation[3] = {sgn * dir[0], sgn * dir[1], sgn * dir[2]};
        for (trials = 0; trials < 20; trials++) {
                if (alouette_decay(
                        0, pid, momentum, polarisation, &products) ==
                    ALOUETTE_RETURN_SUCCESS)
                        break;
        }

        int j;
        for (j = 0; j < products.size; j++) {
                const int pid1 = products.pid[j];
                if (abs(pid1 == 12) || (abs(pid1) == 13) || (abs(pid1) == 14) ||
                    (abs(pid1) == 16))
                        continue;
                momentum[0] = products.P[j][0];
                momentum[1] = products.P[j][1];
                momentum[2] = products.P[j][2];
                record_copy_product(context, pid1, momentum);
        }
        if (unlock != NULL) unlock();
        if (record_publish(context) != EXIT_SUCCESS) return EXIT_FAILURE;
        return EXIT_SUCCESS;
}

/* Concatenate a dirname and a filename. */
char * path_concat(const char * dirname, const char * filename)
{
        const int n0 = strlen(dirname);
        const int n1 = strlen(filename);

        char * s = malloc(n0 + n1 + 2);
        if (s == NULL) return NULL;
        memcpy(s, dirname, n0);
        memcpy(s + n0 + 1, filename, n1 + 1);
        s[n0] = '/';

        return s;
}

static int copy_config_string(
    const char * prefix, const char * suffix, char ** config)
{
        if (*config != NULL) free(*config);
        const int m = strlen(prefix);
        const int n = m + strlen(suffix) + 2;
        *config = malloc(n);
        if (*config == NULL) {
                danton_error_push(NULL,
                    "%s (%d): could not allocate memory.", __FILE__, __LINE__);
                return EXIT_FAILURE;
        }
        char * s = *config;
        memcpy(s, prefix, m);
        s[m] = '/';
        memcpy(s + m + 1, suffix, n - m - 1);

        return EXIT_SUCCESS;
}

/* Initialise the DANTON library. */
int danton_initialise(
    const char * prefix, danton_lock_cb * lock_, danton_lock_cb * unlock_)
{
        /* Clear the error status. */
        errno = 0;

        /* Set the data paths for later use. */
        if (prefix == NULL)
                prefix = DANTON_PREFIX;

        if (copy_config_string(prefix, DANTON_GEOID, &geoid_path) !=
            EXIT_SUCCESS)
                return EXIT_FAILURE;

        /* Check and update the lock callbacks. */
        if (((lock_ == NULL) && (unlock_ != NULL)) ||
            ((lock_ != NULL) && (unlock_ == NULL))) {
                danton_error_push(NULL,
                    "%s (%d): incomplete (un)lock function(s).", __FILE__,
                    __LINE__);
                return EXIT_FAILURE;
        }
        lock = lock_;
        unlock = unlock_;

        /* Set the error handlers for PUMAS and TURTLE */
        pumas_error_handler_set(&catch_pumas_error);
        turtle_error_handler_set(&catch_turtle_error);

        return EXIT_SUCCESS;
}

/* Finalise the DANTON library. */
void danton_finalise(void)
{
        free(geoid_path);
        geoid_path = NULL;
        turtle_stack_destroy(&earth.stack);
        turtle_map_destroy(&earth.undulations);
}

/* Destroy any memory flat object. */
void danton_destroy(void ** any)
{
        free(*any);
        *any = NULL;
}

/* Set the global Earth model. */
int danton_earth_model(const char * reference, const char * topography,
    double density, int * sea)
{
        /* Parse the reference. */
        if (reference != NULL) {
                if (strcmp(reference, "PREM") == 0) {
                        earth.reference = EARTH_PREM;
                } else if (strcmp(reference, "WGS84") == 0) {
                        earth.reference = EARTH_WGS84;
                } else if (strcmp(reference, "EGM96") == 0) {
                        earth.reference = EARTH_EGM96;
                        if (earth.undulations == NULL) {
                                if (turtle_map_load(
                                    &earth.undulations, geoid_path) !=
                                    TURTLE_RETURN_SUCCESS) {
                                        ERROR_TURTLE(NULL);
                                        return EXIT_FAILURE;
                                }
                        }
                } else {
                        danton_error_push(NULL,
                            "%s (%d): unknown reference `%s`", __FILE__,
                            __LINE__, reference);
                        return EXIT_FAILURE;
                }
        }

        /* Destroy any previous topography data. */
        if (earth.stack != NULL) {
                turtle_stack_destroy(&earth.stack);
        }

        if ((earth.reference != EARTH_EGM96) && (earth.undulations != NULL)) {
                turtle_map_destroy(&earth.undulations);
        }

        /* Parse the topography. */
        if (topography != NULL) {
                if (strncmp(topography, "flat://", 7) == 0) {
                        const char * nptr = topography + 7;
                        char * endptr;
                        errno = 0;
                        earth.z0 = strtod(nptr, &endptr);
                        earth.flat_topography = 1;
                        if ((errno != 0) || (endptr == nptr)) {
                                danton_error_push(NULL,
                                    "%s (%d): invalid topography `%s`",
                                    __FILE__, __LINE__, topography);
                                return EXIT_FAILURE;
                        }
                } else {
                        if ((reference != NULL) &&
                            (earth.reference == EARTH_PREM)) {
                                danton_error_push(NULL,
                                    "%s (%d): reference must be `WGS84` "
                                    "or `EGM96` when specifying a detailed "
                                    "topography",
                                    __FILE__, __LINE__);
                                return EXIT_FAILURE;
                        }
                        if (turtle_stack_create(&earth.stack, topography,
                             -1, lock, unlock) !=
                             TURTLE_RETURN_SUCCESS) {
                                ERROR_TURTLE(NULL);
                                return EXIT_FAILURE;
                        }
                        earth.flat_topography = 0;
                }
        } else {
                earth.flat_topography = 0;
        }

        /* Parse the topography density. */
        if (density > 0.) earth.density = density;

        /* Set the sea flag. */
        if (sea != NULL) earth.sea = *sea;

        return EXIT_SUCCESS;
}

/* Low level data structure for an event sampler. */
struct event_sampler {
        struct danton_sampler api;
        char endstr;
        double neutrino_weight;
        double total_weight;
        unsigned long hash;
};

/* Create a new event sampler. */
struct danton_sampler * danton_sampler_create(void)
{
        struct event_sampler * sampler;
        sampler = malloc(sizeof(*sampler));
        if (sampler == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory.",
                    __FILE__, __LINE__);
                return NULL;
        }
        memset(sampler, 0x0, sizeof(*sampler));
        sampler->hash--;
        return &sampler->api;
}

/* Bernstein's djb2 hash function, from
 * http://www.cse.yorku.ca/~oz/hash.html.
 */
static unsigned long hash(unsigned char * str)
{
        unsigned long hash = 5381;
        int c;

        while ((c = *str++)) {
                hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
        }
        return hash;
}

/* Update an event sampler. */
int danton_sampler_update(struct danton_sampler * sampler)
{
        struct event_sampler * sampler_ = (struct event_sampler *)sampler;

        /* Check the latitude. */
        if ((sampler->latitude < -90.) || (sampler->latitude > 90.)) {
                danton_error_push(NULL,
                    "%s (%d): invalid latitude value `%.5g`.", __FILE__,
                    __LINE__, sampler->latitude);
                return EXIT_FAILURE;
        }

        /* Check the longitude. */
        if ((sampler->longitude < -180.) || (sampler->longitude > 1800.)) {
                danton_error_push(NULL,
                    "%s (%d): invalid longitude value `%.5g`.", __FILE__,
                    __LINE__, sampler->longitude);
                return EXIT_FAILURE;
        }

        /* Check the altitude. */
        if (sampler->altitude[0] > sampler->altitude[1]) {
                danton_error_push(NULL, "%s (%d): invalid altitude value(s).",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        /* Check the azimuth angle. */
        if ((sampler->azimuth[0] > sampler->azimuth[1]) ||
            (sampler->azimuth[1] - sampler->azimuth[0] > 360.)) {
                danton_error_push(NULL, "%s (%d): invalid azimuth value(s).",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        /* Check the elevation angle. */
        if ((sampler->elevation[0] < -90.) ||
            (sampler->elevation[0] > sampler->elevation[1]) ||
            (sampler->elevation[1] > 90.)) {
                danton_error_push(NULL, "%s (%d): invalid elevation value(s).",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        /* Check the energy. */
        if ((sampler->energy[0] < 1E+02) ||
            (sampler->energy[0] > sampler->energy[1]) ||
            (sampler->energy[1] > 1E+12)) {
                danton_error_push(NULL, "%s (%d): invalid energy values.",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        /* Compute the particle weights. */
        sampler_->neutrino_weight = 0.;
        int i;
        for (i = 0; i < DANTON_PARTICLE_N_NU; i++)
                sampler_->neutrino_weight +=
                    (sampler->weight[i] <= 0.) ? 0. : sampler->weight[i];
        sampler_->total_weight = sampler_->neutrino_weight;
        if (sampler->weight[DANTON_PARTICLE_TAU_BAR] > 0.)
                sampler_->total_weight +=
                    sampler->weight[DANTON_PARTICLE_TAU_BAR];
        if (sampler->weight[DANTON_PARTICLE_TAU] > 0.)
                sampler_->total_weight += sampler->weight[DANTON_PARTICLE_TAU];

        /* Update the hash and return. */
        sampler_->hash = hash((unsigned char *)sampler);
        return EXIT_SUCCESS;
}

/* Create a new simulation context for DANTON. */
struct danton_context * danton_context_create(void)
{
        struct simulation_context * context;
        context = malloc(sizeof(*context));
        if (context == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory.",
                    __FILE__, __LINE__);
                return NULL;
        }

        /* Initialise the random context. */
        context->random = NULL;

        /* Initialise the Monte-Carlo contexts and the recorder.
         */
        context->ent.medium = &medium_ent;
        context->ent.random = &random_ent;
        context->ent.ancestor = NULL;
        context->ent.distance_max = 0.;
        context->ent.grammage_max = 0.;

        context->pumas = NULL;
        context->client = NULL;
        context->record = NULL;
        context->error.count = 0;
        context->error.size = 0;
        context->error.data[0] = 0;

        /* Initialise the public API data. */
        context->api.run_action = NULL;
        context->api.mode = DANTON_MODE_BACKWARD;
        context->api.longitudinal = 0;
        context->api.decay = 1;
        int i;
        for (i = 0; i < DANTON_PARTICLE_N_NU; i++)
                context->api.primary[i] = NULL;
        context->api.sampler = NULL;
        context->api.recorder = NULL;

        /* The lower (upper) energy bound under (above) which
         * all
         * particles are
         * killed.
         */
        context->energy_cut = -1;

        /* Flag to check if the neutrino flux is requested. */
        context->sample_flux = 0;

        return &context->api;
}

/* Destroy a DANTON simulation context. */
void danton_context_destroy(struct danton_context ** context)
{
        if (*context == NULL) return;
        struct simulation_context * context_ =
            (struct simulation_context *)(*context);
        if (context_->pumas != NULL) {
                pumas_recorder_destroy(&context_->pumas->recorder);
                pumas_context_destroy(&context_->pumas);
        }
        turtle_client_destroy(&context_->client);
        free(context_->record);
        free(context_);
        *context = NULL;
}

/* Get the PDG number corresponding to the given table index. */
int danton_particle_pdg(enum danton_particle index)
{
        int pdg[DANTON_PARTICLE_N] = { ENT_PID_NU_BAR_TAU, ENT_PID_NU_BAR_MU,
                ENT_PID_NU_BAR_E, ENT_PID_NU_E, ENT_PID_NU_MU, ENT_PID_NU_TAU,
                ENT_PID_TAU_BAR, ENT_PID_TAU };
        if ((index < 0) || (index >= DANTON_PARTICLE_N)) return 0;
        return pdg[index];
}

/* Get the table index corresponding to the given PDG number. */
enum danton_particle danton_particle_index(int pdg)
{
        if (pdg > 0) {
                if (pdg == ENT_PID_NU_E)
                        return DANTON_PARTICLE_NU_E;
                else if (pdg == ENT_PID_NU_MU)
                        return DANTON_PARTICLE_NU_MU;
                else if (pdg == ENT_PID_NU_TAU)
                        return DANTON_PARTICLE_NU_TAU;
                else if (pdg == ENT_PID_TAU)
                        return DANTON_PARTICLE_TAU;
        } else {
                if (pdg == ENT_PID_NU_BAR_TAU)
                        return DANTON_PARTICLE_NU_BAR_TAU;
                else if (pdg == ENT_PID_NU_BAR_MU)
                        return DANTON_PARTICLE_NU_BAR_MU;
                else if (pdg == ENT_PID_NU_BAR_E)
                        return DANTON_PARTICLE_NU_BAR_E;
                else if (pdg == ENT_PID_TAU_BAR)
                        return DANTON_PARTICLE_TAU_BAR;
        }
        return DANTON_PARTICLE_UNKNOWN;
}

/* Sample a parameter uniformly over a range. */
static double sample_linear(struct simulation_context * context,
    const double x[2], long i, long n, double * weight)
{
        double xi;
        if (x[0] < x[1]) {
                const double dx = x[1] - x[0];
                double u;
                if ((context->api.mode == DANTON_MODE_GRAMMAGE) && (n > 0))
                        u = (n > 1) ? i / (n - 1.) : 0.;
                else {
                        u = random_uniform01(context);
                        if (weight != NULL) *weight *= dx;
                }
                xi = dx * u + x[0];
        } else
                xi = x[0];
        return xi;
}

/* Sample a parameter log-uniformly over a range. */
static double sample_log_or_linear(
    struct simulation_context * context, double x[2], double * weight)
{
        double xi;
        if (x[0] < x[1]) {
                if ((x[0] > 0.) || (x[1] < 0.)) {
                        const double r = log(x[1] / x[0]);
                        xi = x[0] * exp(r * random_uniform01(context));
                        if (weight != NULL) *weight *= fabs(r) * xi;
                } else {
                        const double dx = x[1] - x[0];
                        xi = x[0] + dx * random_uniform01(context);
                        if (weight != NULL) *weight *= dx;
                }
        } else
                xi = x[0];
        return xi;
}

/* Helper function for creating a PUMAS recorder */
static int create_pumas_recorder(
    struct simulation_context * context, pumas_recorder_cb * record)
{
        struct pumas_recorder * recorder;
        if (pumas_recorder_create(&recorder, 0)
            != PUMAS_RETURN_SUCCESS) {
                ERROR_PUMAS(&context->api);
                return EXIT_FAILURE;
        }
        recorder->period = 0;
        recorder->record = record;
        context->pumas->recorder = recorder;

        return EXIT_SUCCESS;
}

/* Reset the context on a physics change. */
DANTON_API void danton_context_reset(struct danton_context * context)
{
        struct simulation_context * context_ =
            (struct simulation_context *)context;
        pumas_context_destroy(&context_->pumas);
}

/* Run a DANTON simulation. */
int danton_context_run(
    struct danton_context * context, long events, long requested)
{
        /* Unpack the various objects. */
        struct simulation_context * context_ =
            (struct simulation_context *)context;
        struct danton_sampler * sampler = context->sampler;
        struct event_sampler * sampler_ = (struct event_sampler *)sampler;

        /* Check and configure the context according to the API.
         */
        if (sampler == NULL) {
                danton_error_push(context, "%s (%d): no sampler was provided.",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }

        if (sampler_->hash != hash((unsigned char *)sampler)) {
                danton_error_push(context,
                    "%s (%d): sampler has not been updated.", __FILE__,
                    __LINE__);
                return EXIT_FAILURE;
        }

        if (context->recorder == NULL) {
                danton_error_push(context, "%s (%d): no recorder was provided.",
                    __FILE__, __LINE__);
                return EXIT_FAILURE;
        }
        if ((context->mode < 0) || (context->mode >= DANTON_MODE_N)) {
                danton_error_push(context, "%s (%d): invalid run mode (%d).",
                    __FILE__, __LINE__, context->mode);
                return EXIT_FAILURE;
        }

        if (context->mode == DANTON_MODE_GRAMMAGE) {
                if (sampler->elevation[0] == sampler->elevation[1])
                        events = 1;
                else if (events < 2) {
                        danton_error_push(context,
                            "%s (%d): numbers of bins must be 2 or more.",
                            __FILE__, __LINE__);
                        return EXIT_FAILURE;
                }
        } else {
                if (sampler_->total_weight <= 0.) {
                        danton_error_push(context,
                            "%s (%d): no particle to sample.", __FILE__,
                            __LINE__);
                        return EXIT_FAILURE;
                }

                if (context->mode == DANTON_MODE_FORWARD) {
                        context_->sample_flux =
                            sampler->altitude[1] == sampler->altitude[0];
                } else {
                        context_->sample_flux =
                            sampler_->neutrino_weight > 0.;
                }

                if (context->decay) {
                        if (sampler_->neutrino_weight ==
                            sampler_->total_weight) {
                                danton_error_push(context,
                                    "%s (%d): no tau(s) target to decay.",
                                    __FILE__, __LINE__);
                                return EXIT_FAILURE;
                        }
                        if (context->mode == DANTON_MODE_FORWARD) {
                                if (context_->sample_flux) {
                                        danton_error_push(context,
                                            "%s (%d): combining flux "
                                            "and tau decays sampling is not "
                                            "supported.",
                                            __FILE__, __LINE__);
                                        return EXIT_FAILURE;
                                }
                        }
                }
        }

        /* Set the initial projectile in backward mode. */
        int projectile = ENT_PID_NU_TAU;
        if ((context->mode == DANTON_MODE_BACKWARD) &&
            (context->mode != DANTON_MODE_GRAMMAGE)) {
                int j;
                double w = 0.;
                for (j = 0; j < DANTON_PARTICLE_N; j++) {
                        if (sampler->weight[j] > 0.) {
                                w = sampler->weight[j];
                                projectile = danton_particle_pdg(j);
                                break;
                        }
                }

                if (w < sampler_->total_weight) {
                        danton_error_push(context,
                            "%s (%d): sampling multiple particles is not "
                            "allowed in backward mode.", __FILE__, __LINE__);
                        return EXIT_FAILURE;
                }
        }

        if ((!earth.flat_topography) && (lock != NULL)) {
                turtle_client_destroy(&context_->client);
                if (turtle_client_create(&context_->client,
                         earth.stack) != TURTLE_RETURN_SUCCESS) {
                        ERROR_TURTLE(context);
                        return EXIT_FAILURE;
                }
        }

        /* Check for any custom run action and configure accordingly. */
        if (context->run_action != NULL) {
                context_->ent.stepping_action = &stepping_ent;
                if (context_->pumas != NULL) {
                        if (create_pumas_recorder(context_, &stepping_pumas) !=
                            EXIT_SUCCESS)
                                return EXIT_FAILURE;
                }
        } else if (context->mode == DANTON_MODE_FORWARD) {
                context_->ent.stepping_action = &stepping_ent_flux;
                if (context_->pumas != NULL) {
                        if (create_pumas_recorder(
                            context_, &stepping_pumas_flux) != EXIT_SUCCESS)
                                return EXIT_FAILURE;
                }
        } else {
                context_->ent.stepping_action = NULL;
                if (context_->pumas != NULL) {
                        pumas_recorder_destroy(&context_->pumas->recorder);
                }
        }

        if (context->mode != DANTON_MODE_GRAMMAGE) {
                /* Check that there is a valid primary flux. */
                int is_primary = 0;
                int j;
                struct danton_primary ** p;
                for (j = 0, p = context->primary; j < DANTON_PARTICLE_N_NU;
                     j++, p++) {
                        if (*p != NULL) {
                                is_primary = 1;
                                if ((*p)->energy[0] > (*p)->energy[1]) {
                                        danton_error_push(context,
                                            "%s (%d): invalid "
                                            "energy "
                                            "range for "
                                            "primary flux (pid "
                                            "= %d).",
                                            __FILE__, __LINE__,
                                            danton_particle_pdg(j));
                                        return EXIT_FAILURE;
                                }
                        }
                }
                if (!is_primary) {
                        danton_error_push(context, "%s (%d): no primary flux.",
                            __FILE__, __LINE__);
                        return EXIT_FAILURE;
                }

                if (context_->pumas == NULL) {
                        /* Create PUMAS context, if not already done. */
                        if (pumas_context_create(
                            &context_->pumas, physics.pumas, 0) !=
                            PUMAS_RETURN_SUCCESS) {
                                ERROR_PUMAS(context);
                                return EXIT_FAILURE;
                        }
                        context_->pumas->medium = &medium_pumas;
                        context_->pumas->random = &random_pumas;
                        context_->pumas->user_data = context_;
                }
                context_->pumas->mode.energy_loss = PUMAS_MODE_STRAGGLED;
                context_->pumas->mode.decay = PUMAS_MODE_RANDOMISED;
                context_->pumas->mode.scattering = context->longitudinal ?
                        PUMAS_MODE_DISABLED : PUMAS_MODE_MIXED;

                /* Configure any custom run action. */
                if (context->run_action != NULL) {
                        if (create_pumas_recorder(context_, &stepping_pumas) !=
                            EXIT_SUCCESS)
                                return EXIT_FAILURE;
                } else if (context->mode == DANTON_MODE_FORWARD) {
                        if (create_pumas_recorder(
                            context_, &stepping_pumas_flux) != EXIT_SUCCESS)
                                return EXIT_FAILURE;
                }

                /* Create the event record, if not already done.
                 */
                if (context_->record == NULL) {
                        const int buffer_size = 10;
                        context_->record = malloc(sizeof(*context_->record) +
                            buffer_size * sizeof(*context_->record->product));
                        if (context_->record == NULL) {
                                danton_error_push(context,
                                    "%s (%d): could not allocate memory.",
                                    __FILE__, __LINE__);
                                exit(EXIT_FAILURE);
                        }
                        context_->record->buffer_size = buffer_size;
                        context_->record->api.vertex = NULL;
                }

                /* Configure the event record. */
                context_->record->api.primary = &context_->record->primary;
                context_->record->api.final = &context_->record->final;
        }

        /* Configure the event count. */
        if ((context->mode == DANTON_MODE_GRAMMAGE) || (requested <= 0))
                requested = events;
        long n_published = 0;

        /* Compute the generation cosine. */
        double cos_theta[2];
        int l;
        for (l = 0; l < 2; l++)
                cos_theta[l] = cos((90. - sampler->elevation[l]) * M_PI / 180.);

        /* Run the simulation. */
        if (context->mode == DANTON_MODE_FORWARD) {
                /* Check the primary flux and pre-compute some sampling
                 * parameters.
                 */
                double primary_p[DANTON_PARTICLE_N_NU];
                int j;
                struct danton_primary ** p;
                for (j = 0, p = context->primary; j < DANTON_PARTICLE_N_NU;
                     j++, p++) {
                        double d;
                        if (*p == NULL)
                                d = 0.;
                        else {
                                d = (*p)->flux(*p, 0.);
                                if (d < 0.) d = 0.;
                        }
                        if (j == 0)
                                primary_p[0] = d;
                        else
                                primary_p[j] = primary_p[j - 1] + d;
                }
                if (primary_p[DANTON_PARTICLE_N_NU - 1] <= 0.) {
                        danton_error_push(context,
                            "%s (%d): null primary flux.", __FILE__, __LINE__);
                        return EXIT_FAILURE;
                }

                /* Run a bunch of forward Monte-Carlo events. */
                context_->ent.ancestor = NULL;
                context_->energy_cut = context->sampler->energy[0];
                if (context_->pumas != NULL) {
                        context_->pumas->mode.direction = PUMAS_MODE_FORWARD;
                        context_->pumas->limit.energy =
                            context_->energy_cut - tau_mass;
                        context_->pumas->event = PUMAS_EVENT_LIMIT_ENERGY;
                }

                long i;
                for (i = 0; (i < events) && (n_published < requested); i++) {
                        /* Sample the projection of the primary state
                         * uniformly.
                         */
                        double weight = 1.;
                        const double ct =
                            sample_linear(context_, cos_theta, i, 0, &weight);
                        const double azimuth = sample_linear(
                            context_, sampler->azimuth, i, 0, &weight);
                        if (sampler->azimuth[1] > sampler->azimuth[0])
                                weight *= M_PI / 180.;
                        const double z0 = sampler->altitude[0];
                        double ecef0[3], u0[3];
                        compute_ecef_position(
                            sampler->latitude, sampler->longitude, z0, ecef0);
                        compute_ecef_direction(sampler->latitude,
                            sampler->longitude, azimuth, ct, u0);

                        int flux_direction;
                        if (fabs(ct) <= FLT_EPSILON) {
                                flux_direction = 0;
                        } else {
                                flux_direction = (ct > 0.) ? 1 : -1;
                        }

                        /* Backward translate the primary state. */
                        double a, b, r2;
                        ellipsoid_parameters_intersection(
                            ecef0, u0, &a, &b, &r2);
                        b = -b;
                        const double ri = 1. + 1.E+05 / PREM_EARTH_RADIUS;
                        const double d2 = b * b + a * (ri * ri - r2);
                        const double d = (d2 <= 0.) ? 0. : sqrt(d2);
                        const double ds = (d - b) / a;
                        ecef0[0] -= ds * u0[0];
                        ecef0[1] -= ds * u0[1];
                        ecef0[2] -= ds * u0[2];

                        /* Sample the primary flavour and its
                         * energy. */
                        int j;
                        double we = 0., energy;
                        while (we == 0.) {
                                const double u = random_uniform01(context_) *
                                    primary_p[DANTON_PARTICLE_N_NU - 1];
                                struct danton_primary ** p;
                                for (j = 0, p = context->primary;
                                     j < DANTON_PARTICLE_N_NU - 1; j++, p++)
                                        if (u <= primary_p[j]) break;
                                we = 1.;
                                energy = sample_log_or_linear(
                                    context_, (*p)->energy, &we);
                                if ((we > 0.) &&
                                    ((*p)->energy[0] < (*p)->energy[1]))
                                        we *= (*p)->flux(*p, energy);
                        }
                        weight *= we;
                        const int pid = danton_particle_pdg(j);

                        /* Configure the primary state. */
                        const int crossed = context_->sample_flux ? 0 : -1;
                        struct generic_state state = {
                                .base.ent = { pid, energy, 0., 0., weight,
                                    { ecef0[0], ecef0[1], ecef0[2] },
                                    { u0[0], u0[1], u0[2] } },
                                .context = context_,
                                .medium = -1,
                                .density = 0.,
                                .x = 0.,
                                .is_tau = 0,
                                .is_inside = -1,
                                .has_crossed = crossed,
                                .flux_direction = flux_direction
                        };

                        /* Initialise the event record. */
                        context_->published = 0;
                        context_->record->api.id = i;
                        context_->record->api.weight = weight;
                        context_->record->api.vertex = NULL;
                        context_->record->api.n_products = 0;
                        record_copy_ent(
                            context_->record->api.primary, &state.base.ent);

                        /* Call any custom initial run action. */
                        if (context->run_action != NULL) {
                                medium(state.base.ent.position,
                                    state.base.ent.direction, &state);
                                struct danton_run_action * action =
                                    context->run_action;
                                if (action->call(context, action,
                                        DANTON_RUN_EVENT_START, state.medium,
                                        context_->record->api.primary) !=
                                    EXIT_SUCCESS)
                                        return EXIT_FAILURE;
                        }

                        /* Do the Monte-Carlo simulation. */
                        if (transport_forward(context_,
                                (struct ent_state *)&state, 1) != EXIT_SUCCESS)
                                return EXIT_FAILURE;

                        /* Update the published count. */
                        if (context_->published) n_published++;

                        /* Call any custom final run action. */
                        if (context->run_action != NULL) {
                                struct danton_run_action * action =
                                    context->run_action;
                                if (action->call(context, action,
                                        DANTON_RUN_EVENT_STOP, -1,
                                        context_->record->api.final) !=
                                    EXIT_SUCCESS)
                                        return EXIT_FAILURE;
                        }
                }
        } else {
                /* Run a bunch of backward Monte-Carlo events.
                 */
                context_->ent.ancestor = &ancestor_cb;
                context_->energy_cut = context->sampler->energy[1];
                int j;
                struct danton_primary ** p;
                for (j = 0, p = context->primary; j < DANTON_PARTICLE_N_NU;
                     j++, p++) {
                        if ((*p != NULL) &&
                            ((*p)->energy[1] > context_->energy_cut))
                                context_->energy_cut = (*p)->energy[1];
                }
                if (context_->pumas != NULL) {
                        context_->pumas->mode.direction = PUMAS_MODE_BACKWARD;
                        context_->pumas->limit.energy =
                            context_->energy_cut - tau_mass;
                        context_->pumas->event = PUMAS_EVENT_LIMIT_ENERGY |
                            PUMAS_EVENT_LIMIT_GRAMMAGE;
                }

                long i;
                for (i = 0; (i < events) && (n_published < requested); i++) {
                        double weight = 1.;
                        const double ct = sample_linear(
                            context_, cos_theta, i, events, &weight);
                        const double azimuth = sample_linear(
                            context_, sampler->azimuth, i, 0, &weight);
                        if (sampler->azimuth[1] > sampler->azimuth[0])
                                weight *= M_PI / 180.;
                        const double energy = sample_log_or_linear(
                            context_, sampler->energy, &weight);
                        const double z0 = sample_log_or_linear(
                            context_, sampler->altitude, &weight);

                        context_->published = 0;

                        if (context->mode != DANTON_MODE_GRAMMAGE) {
                                context_->record->api.id = i;
                                context_->record->api.generation = 1;
                                context_->record->api.vertex = NULL;
                                context_->record->api.n_products = 0;
                        }
                        if ((context->mode != DANTON_MODE_GRAMMAGE) &&
                            !context_->sample_flux) {
                                /* This is a particle Monte-Carlo. */
                                const double charge =
                                    (projectile > 0) ? -1. : 1.;
                                double ecef0[3], u0[3];
                                compute_ecef_position(sampler->latitude,
                                    sampler->longitude, z0, ecef0);
                                compute_ecef_direction(sampler->latitude,
                                    sampler->longitude, azimuth, ct, u0);
                                struct generic_state state = {
                                        .base.pumas = { charge,
                                            energy - tau_mass, 0., 0., 0.,
                                            weight,
                                            { ecef0[0], ecef0[1], ecef0[2] },
                                            { u0[0], u0[1], u0[2] }, 0 },
                                        .context = context_,
                                        .medium = -1,
                                        .density = 0.,
                                        .x = 0.,
                                        .is_tau = 1,
                                        .is_inside = -1,
                                        .has_crossed = -1,
                                        .flux_direction = 0
                                };

                                /* Call any custom initial run action. */
                                if (context->run_action != NULL) {
                                        medium(state.base.pumas.position,
                                            state.base.pumas.direction, &state);
                                        struct danton_state s;
                                        record_copy_pumas(
                                            &s, &state.base.pumas);
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_START,
                                                state.medium,
                                                &s) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                                if (transport_backward(context_, &state) !=
                                    EXIT_SUCCESS)
                                        return EXIT_FAILURE;

                                /* Update the published count. */
                                if (context_->published) n_published++;

                                /* Call any custom final run action. */
                                if (context->run_action != NULL) {
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_STOP, -1,
                                                context_->record->api
                                                    .primary) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                        } else if ((context->mode != DANTON_MODE_GRAMMAGE) &&
                            context_->sample_flux) {
                                double ecef0[3], u0[3];
                                compute_ecef_position(sampler->latitude,
                                    sampler->longitude, z0, ecef0);
                                compute_ecef_direction(sampler->latitude,
                                    sampler->longitude, azimuth, ct, u0);
                                struct generic_state state = {
                                        .base.ent = { projectile, energy, 0.,
                                            0., weight,
                                            { ecef0[0], ecef0[1], ecef0[2] },
                                            { u0[0], u0[1], u0[2] } },
                                        .context = context_,
                                        .medium = -1,
                                        .density = 0.,
                                        .x = 0.,
                                        .is_tau = 0,
                                        .is_inside = -1,
                                        .has_crossed = -1,
                                        .flux_direction = 0
                                };

                                /* Call any custom initial run action. */
                                if (context->run_action != NULL) {
                                        medium(state.base.ent.position,
                                            state.base.ent.direction, &state);
                                        struct danton_state s;
                                        record_copy_ent(&s, &state.base.ent);
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_START,
                                                state.medium,
                                                &s) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                                if (transport_backward(context_, &state) !=
                                    EXIT_SUCCESS)
                                        return EXIT_FAILURE;

                                /* Update the published count. */
                                if (context_->published) n_published++;

                                /* Call any custom final run action. */
                                if (context->run_action != NULL) {
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_STOP, -1,
                                                context_->record->api
                                                    .primary) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                        } else {
                                /* This is a grammage scan using
                                 * a non-interacting neutrino. First
                                 * let us initialise the neutrino state.
                                 */
                                double ecef0[3], u0[3];
                                compute_ecef_position(sampler->latitude,
                                    sampler->longitude, z0, ecef0);
                                compute_ecef_direction(sampler->latitude,
                                    sampler->longitude, azimuth, ct, u0);
                                struct generic_state g_state = {
                                        .base.ent = { projectile, energy, 0.,
                                            0., weight,
                                            { ecef0[0], ecef0[1], ecef0[2] },
                                            { u0[0], u0[1], u0[2] } },
                                        .context = context_,
                                        .medium = -1,
                                        .density = 0.,
                                        .x = 0.,
                                        .is_tau = 0,
                                        .is_inside = -1,
                                        .has_crossed = -1,
                                        .flux_direction = 0
                                };

                                /* Call any custom initial run action. */
                                struct ent_state * state = &g_state.base.ent;
                                if (context->run_action != NULL) {
                                        medium(g_state.base.ent.position,
                                            g_state.base.ent.direction,
                                            &g_state);
                                        struct danton_state s;
                                        record_copy_ent(&s, state);
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_START,
                                                g_state.medium,
                                                &s) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                                /* Then let us do the transport
                                 * with ENT.
                                 */
                                enum ent_event event = ENT_EVENT_NONE;
                                while (event != ENT_EVENT_EXIT) {
                                        ent_transport(NULL, &context_->ent,
                                            state, NULL, &event);
                                }

                                /* Call any custom final run action. */
                                if (context->run_action != NULL) {
                                        struct danton_state s;
                                        record_copy_ent(&s, state);
                                        struct danton_run_action * action =
                                            context->run_action;
                                        if (action->call(context, action,
                                                DANTON_RUN_EVENT_STOP, -1,
                                                &s) != EXIT_SUCCESS)
                                                return EXIT_FAILURE;
                                }

                                /* Finally, let us publish the
                                 * result.
                                 */
                                struct danton_grammage g = { 90. -
                                            acos(ct) / M_PI * 180.,
                                        state->grammage };
                                if (context->recorder->record_grammage(context,
                                        context->recorder, &g) != EXIT_SUCCESS)
                                        return EXIT_FAILURE;
                        }
                }
        }

        return EXIT_SUCCESS;
}

/* Global error buffer. */
static struct error_stack g_error = { 0, 0, "" };

/* Helper function for getting the relevant error stack. */
static struct error_stack * error_stack_get(struct danton_context * context)
{
        if (context == NULL)
                return &g_error;
        else {
                struct simulation_context * c =
                    (struct simulation_context *)context;
                return &c->error;
        }
}

/* API function for getting the number of error(s) in the stack.
 */
int danton_error_count(struct danton_context * context)
{
        struct error_stack * error = error_stack_get(context);
        return error->count;
}

/* API function for pushing an error message to the stack. */
int danton_error_push(struct danton_context * context, const char * format, ...)
{
        struct error_stack * error = error_stack_get(context);
        const int left = ERROR_SIZE - error->size;

        va_list args;
        va_start(args, format);
        const int n = vsnprintf(error->data + error->size, left, format, args);
        va_end(args);
        if (n >= left) return EXIT_FAILURE;
        error->size += n + 1;
        error->count++;
        return EXIT_SUCCESS;
}

/* API function for getting the last error message from the
 * stack. */
const char * danton_error_pop(struct danton_context * context)
{
        struct error_stack * error = error_stack_get(context);
        if (error->size == 0) return NULL;

        char * s;
        int n;
        for (s = error->data + error->size - 2, n = error->size - 2;
             (n > 0) && (*s != 0); s--, n--)
                ;
        if (n > 0) s++, n++;
        error->size = n;
        error->count--;
        return s;
}

/* Topography interface. */

DANTON_API double danton_geoid_undulation(
    double latitude, double longitude)
{
        if (earth.undulations == NULL) {
                return 0.0;
        } else {
                double undulation = 0.0;
                int inside = 0;
                while (longitude < 0.0) {
                        longitude += 360.0;
                }
                while (longitude > 360.0) {
                        longitude -= 360.0;
                }
                turtle_map_elevation(
                    earth.undulations,
                    longitude,
                    latitude,
                    &undulation,
                    &inside
                );
                return inside ? undulation : 0.0;
        }
}

DANTON_API double danton_topography_elevation(
    double latitude, double longitude)
{
        if (earth.flat_topography) {
                return earth.z0;
        } else if (earth.stack == NULL) {
                return 0.0;
        } else {
                double elevation = 0.0;
                int inside = 0;
                turtle_stack_elevation(
                    earth.stack,
                    longitude,
                    latitude,
                    &elevation,
                    &inside
                );
                return inside ? elevation : 0.0;
        }
}

/* Tracer interface. */

struct danton_tracer {
        struct generic_state state;
        struct danton_context context;
};

DANTON_API struct danton_tracer * danton_tracer_create(void) {
        struct danton_tracer * tracer;
        tracer = malloc(sizeof(*tracer));
        if (tracer == NULL) {
                danton_error_push(NULL, "%s (%d): could not allocate memory.",
                    __FILE__, __LINE__);
                return NULL;
        }
        tracer->state.is_tau = 0;
        tracer->state.context = (void *)&tracer->context;
        tracer->context.decay = 1;
        return tracer;
}

DANTON_API int danton_tracer_medium(
    struct danton_tracer * tracer, const double * position)
{
        double direction[3] = { 0.0, 0.0, 1.0 };
        medium(position, direction, &tracer->state);
        return tracer->state.medium;
}

static int tracer_medium(
    struct danton_tracer * tracer,
    const double * r,
    const double * u,
    double * step)
{
        double s = medium(r, u, &tracer->state);
        if (step != NULL) *step = s;
        return tracer->state.medium;
}

static int tracer_step(
    struct danton_tracer * tracer,
    const double * ri,
    const double * u,
    double * step,
    int * medium)
{
        double si;
        int mi = tracer_medium(tracer, ri, u, &si);
        if (mi < 0) {
                *step = 0.0;
                *medium = mi;
                return mi;
        }

        double rf[3] = {
            ri[0] + si * u[0],
            ri[1] + si * u[1],
            ri[2] + si * u[2]
        };
        double sf;
        int mf = tracer_medium(tracer, rf, u, &sf);
        if (mf == mi) {
                *step = si;
                *medium = mf;
                return mi;
        }

#define STEP_MIN 1E-03
        const double step_min = (si < STEP_MIN) ? si : STEP_MIN;
        double s1 = 0., s2 = -step_min;
        double r[3];
        r[0] = rf[0] + s2 * u[0];
        r[1] = rf[1] + s2 * u[1];
        r[2] = rf[2] + s2 * u[2];
        int m = tracer_medium(tracer, r, u, NULL);
        if (m != mi) {
                /* Locate the medium change by dichotomy. */
                if (m != mf) mf = m;
                s1 = s2;
                s2 = -si;
                while (fabs(s1 - s2) > STEP_MIN) {
                        double s3 = 0.5 * (s1 + s2);
                        r[0] = rf[0] + s3 * u[0];
                        r[1] = rf[1] + s3 * u[1];
                        r[2] = rf[2] + s3 * u[2];
                        m = tracer_medium(tracer, r, u, NULL);
                        if (m == mi) {
                                s2 = s3;
                        } else {
                                s1 = s3;
                                /* Update the end medium if required. */
                                if (m != mf) mf = m;
                        }
                }
                si += s1;
        }
        *step = si;
        *medium = mf;
        return mi;
#undef STEP_MIN
}

DANTON_API int danton_tracer_trace(
    struct danton_tracer * tracer,
    const double * position,
    const double * direction,
    double * distance,
    int * next_medium)
{
        double r[3];
        memcpy(r, position, sizeof r);
        double d = 0.0;
        int mi;
        for (;;) {
                int mf;
                double s;
                mi = tracer_step(tracer, r, direction, &s, &mf);
                if ((mi < 0) || (s <= 0.0)) break;
                d += s;
                if (mf != mi) {
                        *next_medium = mf;
                        break;
                }
                r[0] += s * direction[0];
                r[1] += s * direction[1];
                r[2] += s * direction[2];
        }
        *distance = d;
        return mi;
}
