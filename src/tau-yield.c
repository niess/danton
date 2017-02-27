/* Standard library includes. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* The Physics APIs. */
#include "ent.h"
#include "pumas.h"
#include "tauola-c.h"

/* The spherical Earth radius, in m. */
#define EARTH_RADIUS 6371.E+03

/* The radius of the geostationary orbit, in m. */
#define GEO_ORBIT 42164E+03

/* The lower energy bound. */
#define ENERGY_MIN 1E+03

/* Handles for the transport engines. */
static struct ent_physics * physics = NULL;
static struct pumas_context * ctx_pumas = NULL;

/* The tau lepton mass, in GeV / c^2. */
static double tau_mass;

/* Finalise and exit to the OS. */
static void gracefully_exit(int rc)
{
        /* Finalise and exit to the OS. */
        ent_physics_destroy(&physics);
        pumas_context_destroy(&ctx_pumas);
        pumas_finalise();
        tauola_finalise();
        exit(rc);
}

/* Error handler for ENT. */
static void handle_ent(enum ent_return rc, ent_function_t * caller)
{
        /* Dump an error message. */
        ent_error_print(stderr, rc, caller, "\t", "\n");
        fprintf(stderr, "\n");

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_FAILURE);
}

/* Error handler for PUMAS. */
static void handle_pumas(enum pumas_return rc, pumas_function_t * caller,
	struct pumas_error * error)
{
	/* Dump the error summary. */
	fprintf(stderr, "error : ");
	pumas_error_print(stderr, rc, caller, error);
	fprintf(stderr, "\n");

	 /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_FAILURE);
}

/* Density according to the Preliminary Earth Model (PEM). */
static double pem_model0(double r, double * density)
{
	*density = 1.02E+03;
	return 0.;
}

static double pem_model1(double r, double * density)
{
	*density = 2.6E+03;
	return 0.;
}

static double pem_model2(double r, double * density)
{
	*density = 2.9E+03;
	return 0.;
}

static double pem_model3(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 0.6924E+03;
	*density = 2.691E+03 + a * x;
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model4(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 3.8045E+03;
	*density = 7.1089E+03 - a * x;
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model5(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 8.0298E+03;
	*density = 11.2494E+03 - a * x;
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model6(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 1.4836E+03;
	*density = 5.3197E+03 - a * x;
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model7(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 6.4761E+03;
	*density = 7.9565E+03 + x * (-a + x * (2.5283E+03 - x * 3.0807E+03));
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model8(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a = 1.2638E+03;
	*density = 12.58155E+03 + x * (-a + x * (-3.6426E+03 - x * 5.5281E+03));
	return 0.01 * EARTH_RADIUS / a;
}

static double pem_model9(double r, double * density)
{
	const double x = r / EARTH_RADIUS;
	const double a2 = -8.8381E+03;
	*density = 13.0885E+03 + a2 * x * x;
	const double xg = (x <= 5E-02) ? 5E-02 : x;
	return 0.01 * EARTH_RADIUS / fabs(2. * a2 * xg);
}

/* The U.S. standard atmosphere model. */
#define USS_MODEL(INDEX, B, C)\
static double uss_model ## INDEX (double r, double * density)\
{\
	*density = B / C * exp(-r / C);\
	return 0.01 * C;\
}

USS_MODEL(0, 12226.562, 9941.8638)
USS_MODEL(1, 11449.069, 8781.5355)
USS_MODEL(2, 13055.948, 6361.4304)
USS_MODEL(3, 5401.778, 7721.7016)

/* Outer space density model. */
static double space_model0(double r, double * density)
{
	*density = 1.E-21; /* ~10^6 H per m^-3. */
	return 0.;
}

/* Generic Monte-Carlo state with stepping data. */
struct generic_state {
	union {
		struct ent_state ent;
		struct pumas_state pumas;
	} base;
	double r;
};

/* Density callbacks for ENT. */
#define DENSITY(MODEL, INDEX)\
static double density_ ## MODEL ## INDEX (\
    struct ent_medium * medium, struct ent_state * state, double * density)\
{\
        struct generic_state * s = (struct generic_state *)state;\
        return MODEL ## _model ## INDEX (s->r, density);\
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
#define LOCALS(MODEL, INDEX)\
static double locals_ ## MODEL ## INDEX (const struct pumas_state * state,\
    struct pumas_locals * locals)\
{\
        struct generic_state * s = (struct generic_state *)state;\
        const double step = MODEL ## _model ## INDEX (s->r, &locals->density);\
        memset(locals->magnet, 0x0, sizeof(locals->magnet));\
        return step;\
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

/* Generic medium callback. */
static double medium(const double * position, const double * direction,
    int * index, struct generic_state * state)
{
        *index = -1;
        double step = 0.;

        const double r2 = position[0] * position[0] + position[1] * position[1]
            + position[2] * position[2];
        if (r2 > GEO_ORBIT * GEO_ORBIT) return step;
	const double r = sqrt(r2);
	state->r = r;

	const double ri[] = {1221.5E+03, 3480.E+03, 5701.E+03,
	    5771.E+03, 5971.E+03, 6151.E+03, 6346.6E+03, 6356.E+03, 6368.E+03,
	    EARTH_RADIUS, EARTH_RADIUS + 4.E+03, EARTH_RADIUS + 1.E+04,
	    EARTH_RADIUS + 4.E+04, EARTH_RADIUS + 1.E+05, GEO_ORBIT,
	    2*GEO_ORBIT};
	int i;
	for (i = 0; i < sizeof(ri) / sizeof(*ri) - 1; i++) {
		if (r <= ri[i]) {
			*index = i;

			/* Outgoing intersection. */
			const double b = position[0] * direction[0] +
			    position[1] * direction[1] +
			    position[2] * direction[2];
			const double r1 = ri[i + 1];
			const double d2 = b * b + r1 * r1 - r * r;
			const double d = (d2 <= 0.) ? 0. : sqrt(d2);
			step = d - b - 1.;

			if ((i > 0) && (b < 0.)) {
				/* This is a downgoing trajectory. First, let
				 * us compute the intersection with the lower
				 * radius.
				 */
				const double r1 = ri[i - 1];
				const double d2 = b * b + r1 * r1 - r * r;
				if (d2 > 0.) {
					const double d = sqrt(d2);
					const double s = d - b - 1.;
					if (s < step) step = s;
				}

				if (i > 1) {
					/* Let us check for an intersection with
					 * the below lower radius.
					 */
					const double r1 = ri[i - 2];
					const double d2 =
					    b * b + r1 * r1 - r * r;
					if (d2 > 0.) {
						const double d = sqrt(d2);
						const double s = d - b - 1.;
						if (s < step) step = s;
					}
				}
			}
			if (step < 1.) step = 1.;
			break;
		}
	}
	return step;
}

/* Medium callback encapsulation for ENT. */
static double medium_ent(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium_ptr)
{
#define ZR 13.
#define AR 26.
#define ZA 7.26199
#define AA 14.5477

        static struct ent_medium media[] = {
		{ ZR, AR, &density_pem0 },
		{ ZR, AR, &density_pem1 },
		{ ZR, AR, &density_pem2 },
		{ ZR, AR, &density_pem3 },
		{ ZR, AR, &density_pem4 },
		{ ZR, AR, &density_pem5 },
		{ ZR, AR, &density_pem6 },
		{ ZR, AR, &density_pem7 },
		{ ZR, AR, &density_pem8 },
		{ ZR, AR, &density_pem9 },
		{ ZA, AA, &density_uss0 },
		{ ZA, AA, &density_uss1 },
		{ ZA, AA, &density_uss2 },
		{ ZA, AA, &density_uss3 },
		{ ZA, AA, &density_space0 }};
        int index;
        const double step = medium(state->position, state->direction, &index,
            (struct generic_state *)state);
        if (index >= 0) *medium_ptr = media + index;
        else *medium_ptr = NULL;
        return step;

#undef  ZR
#undef  AR
#undef  ZA
#undef  AA
}

/* Medium callback encapsulation for PUMAS. */
double medium_pumas(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr)
{
        static struct pumas_medium media[] = {
		{ 0, &locals_pem0 },
		{ 0, &locals_pem1 },
		{ 0, &locals_pem2 },
		{ 0, &locals_pem3 },
		{ 0, &locals_pem4 },
		{ 0, &locals_pem5 },
		{ 0, &locals_pem6 },
		{ 0, &locals_pem7 },
		{ 0, &locals_pem8 },
		{ 0, &locals_pem9 },
		{ 1, &locals_uss0 },
		{ 1, &locals_uss1 },
		{ 1, &locals_uss2 },
		{ 1, &locals_uss3 },
		{ 1, &locals_space0 }};
        int index;
        const double step = medium(state->position, state->direction, &index,
            (struct generic_state *)state);
        if (index >= 0) *medium_ptr = media + index;
        else *medium_ptr = NULL;
        return step;
}

/* Uniform distribution over [0,1]. */
static double random(void * context)
{
        return rand() / (double)RAND_MAX;
}

/* Loader for PUMAS. */
void load_pumas()
{
	const enum pumas_particle particle = PUMAS_PARTICLE_TAU;
	const char * dump = "materials.b";

	/* First, attempt to load any binary dump. */
	FILE * stream = fopen(dump, "rb");
        if (stream != NULL) {
                pumas_error_catch(1); /* catch any error. */
                pumas_load(stream);
                fclose(stream);
                pumas_error_raise();
		pumas_particle(NULL, NULL, &tau_mass);
                return;
        }

        /* If no binary dump, initialise from the MDF and dump. */
        pumas_initialise(particle, NULL, NULL, NULL);
        pumas_particle(NULL, NULL, &tau_mass);

        /* Dump the library configuration. */
        stream = fopen(dump, "wb+");
        if (stream == NULL) handle_pumas(PUMAS_RETURN_IO_ERROR, NULL, NULL);
        pumas_dump(stream);
}

/* Set a neutrino state from a tau decay product. */
static void copy_neutrino(const struct pumas_state * tau, int pid,
    const double * momentum, struct ent_state * neutrino)
{
	neutrino->pid = pid;
	neutrino->energy = sqrt(momentum[0] * momentum[0] +
	    momentum[1] * momentum[1] + momentum[2] * momentum[2]);
	memcpy(neutrino->position, tau->position, sizeof(neutrino->position));
	neutrino->direction[0] = momentum[0] / neutrino->energy;
	neutrino->direction[1] = momentum[1] / neutrino->energy;
	neutrino->direction[2] = momentum[2] / neutrino->energy;
	neutrino->distance = tau->distance;
	neutrino->grammage = tau->grammage;
	neutrino->weight = tau->weight;
}

/* Transport routine, recursive. */
static void transport(struct ent_context * ctx_ent, struct ent_state * neutrino,
    int eventid, int generation, FILE * stream)
{
	if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
	    (abs(neutrino->pid) != ENT_PID_NU_TAU)) return;

	struct ent_state product;
	enum ent_event event;
	for (;;) {
		/* Neutrino transport with ENT. */
		ent_transport(physics, ctx_ent, neutrino, &product, &event);
		if ((event == ENT_EVENT_EXIT) ||
		    (neutrino->energy <= ENERGY_MIN)) break;

		if (abs(product.pid) == ENT_PID_TAU) {
			/* Tau transport with PUMAS. */
			const double charge = (product.pid > 0) ? -1. : 1.;
			const double kinetic = product.energy - tau_mass;
			struct generic_state tau_data = {.base.pumas = {charge,
			    kinetic, product.distance, product.grammage, 0.,
			    product.weight}, .r = 0.};
			struct pumas_state * tau = &tau_data.base.pumas;
			memcpy(&tau->position, &product.position,
			    sizeof(tau->position));
			memcpy(&tau->direction, &product.direction,
			    sizeof(tau->direction));
			pumas_transport(ctx_pumas, tau);
			if (tau->decayed) {
				/* Tau decay with TAUOLA. */
				const double p = sqrt(tau->kinetic * (
				    tau->kinetic + 2. * tau_mass));
				double momentum[3] = {p * tau->direction[0],
				    p * tau->direction[1],
				    p * tau->direction[2]};
				 tauola_decay(product.pid, momentum,
				    tau->direction);
				int pid, nprod = 0;
				struct generic_state nu_e_data, nu_t_data;
				struct ent_state * nu_e = NULL, * nu_t = NULL;
				while (tauola_product(&pid, momentum)) {
				    if (abs(pid) == 16) {
					    /* Update the neutrino state with
					     * the nu_tau daughter.
					     */
					    if (neutrino->pid == ENT_PID_HADRON)
					            copy_neutrino(tau, pid,
					                momentum, neutrino);
					    else {
						    nu_t = &nu_t_data.base.ent;
					            copy_neutrino(tau, pid,
					                momentum, nu_t);
				            }
					    continue;
				    }
				    else if (pid == -12) {
					    nu_e = &nu_e_data.base.ent;
					    copy_neutrino(tau, pid, momentum,
					        nu_e);
					    continue;
				    }
				    else if ((pid == 12) || (abs(pid) == 13) ||
				        (abs(pid) == 14)) continue;
				    if (nprod == 0) fprintf(stream,
				        "%5d %2d %4d %12.5lE %12.3lf %12.3lf "
				        "%12.3lf\n", eventid, generation,
				        product.pid, tau->kinetic,
				        tau->position[0], tau->position[1],
				        tau->position[2]);
				    fprintf(stream, "    %4d %12.5lE %12.5lE "
				    "%12.5lE\n", pid, momentum[0], momentum[1],
				    momentum[2]);
				    nprod++;
				}
				generation++;

				/* Process any additional nu_e~ or nu_tau. */
				if (nu_e != NULL) transport(ctx_ent, nu_e,
				    eventid, generation, stream);
				if (nu_t != NULL) transport(ctx_ent, nu_t,
				    eventid, generation, stream);
			}
		}
		if ((neutrino->pid != ENT_PID_NU_BAR_E) &&
		    (abs(neutrino->pid) != ENT_PID_NU_TAU)) break;
	}
}

int main(int nargc, char * argv[])
{
        /* Set the input arguments. */
        double energy = 1E+08;
        enum ent_pid projectile = ENT_PID_NU_TAU;
        int events = 100;

        /* Register the error handlers. */
        ent_error_handler_set(&handle_ent);
        pumas_error_handler_set(&handle_pumas);

        /* Create a new neutrino Physics environment. */
        ent_physics_create(&physics, "ent/data/pdf/CT14nnlo_0000.dat");

        /* Initialise the PUMAS transport engine. */
        load_pumas();

        /* Initialise TAUOLA. */
        tauola_initialise(1, NULL);

        /* Initialise the Monte-Carlo contexts. */
        struct ent_context ctx_ent = { &medium_ent,
	    (ent_random_cb *)&random, 1};
        pumas_context_create(0, &ctx_pumas);
        ctx_pumas->medium = &medium_pumas;
        ctx_pumas->random = (pumas_random_cb *)&random;
        ctx_pumas->kinetic_limit = ENERGY_MIN - tau_mass;

        /* Run a batch of Monte-Carlo events. */
        FILE * stream = fopen("transport.dat", "w+");
        int i;
        for (i = 0; i < events; i++) {
                struct generic_state state = {.base.ent = { projectile,
		    energy, 0., 0., 1., { 0., 0., -EARTH_RADIUS - 1E+05 },
		    { 0., 0., 1. }}, .r = 0.};
                transport(&ctx_ent, (struct ent_state *)&state, i, 1, stream);

        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_SUCCESS);
}
