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
#define EARTH_RADIUS 6.4E+06

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

/* Density callback for ENT, with a uniform density. */
static double density(
    struct ent_medium * medium, struct ent_state * state, double * density)
{
        *density = 2.65E+03;
        return 0.;
}

/* Locals callback for PUMAs. */
static double locals(const struct pumas_state * state,
    struct pumas_locals * locals)
{
        static struct pumas_locals locals_ = {2.65E+03, {0., 0., 0.}};
        memcpy(locals, &locals_, sizeof(locals));
        return 0.;
}

/* Generic medium callback, with a single medium. */
static double medium(const double * r, int * index)
{
        const double r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        if (r2 <= EARTH_RADIUS * EARTH_RADIUS) *index = 0;
	else *index = -1;
        return EARTH_RADIUS;
}

/* Medium callback encapsulation for ENT. */
static double medium_ent(struct ent_context * context, struct ent_state * state,
    struct ent_medium ** medium_ptr)
{
        static struct ent_medium medium_ = { 13., 26., &density };
        int index;
        const double step = medium(state->position, &index);
        if (index >= 0) *medium_ptr = &medium_;
        else *medium_ptr = NULL;
        return step;
}

/* Medium callback encapsulation for PUMAS. */
double medium_pumas(struct pumas_context * context,
    struct pumas_state * state, struct pumas_medium ** medium_ptr)
{
        static struct pumas_medium medium_ = { 0, &locals };
        int index;
        const double step = medium(state->position, &index);
        if (index >= 0) *medium_ptr = &medium_;
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
			struct pumas_state tau = {charge, kinetic,
			    product.distance, product.grammage, 0.,
			        product.weight};
			memcpy(&tau.position, &product.position,
			    sizeof(tau.position));
			memcpy(&tau.direction, &product.direction,
			    sizeof(tau.direction)); 
			pumas_transport(ctx_pumas, &tau);
			if (tau.decayed) {
				/* Tau decay with TAUOLA. */
				const double p = sqrt(tau.kinetic * (
				    tau.kinetic + 2. * tau_mass));
				double momentum[3] = {p * tau.direction[0],
				    p * tau.direction[1], p * tau.direction[2]};
				 tauola_decay(product.pid, momentum,
				    tau.direction);
				int pid, is_nue = 0, is_nut = 0, nprod = 0;
				struct ent_state nu_e, nu_t;
				while (tauola_product(&pid, momentum)) {
				    if (abs(pid) == 16) {
					    /* Update the neutrino state with
					     * the nu_tau daughter.
					     */
					    if (neutrino->pid == ENT_PID_HADRON)
					            copy_neutrino(&tau, pid,
					                momentum, neutrino);
					    else {
					            copy_neutrino(&tau, pid,
					                momentum, &nu_t);
					            is_nut = 1;
				            }
					    continue;
				    }
				    else if (pid == -12) {
					    is_nue = 1;
					    copy_neutrino(&tau, pid, momentum,
					        &nu_e);
					    continue;
				    }
				    else if ((pid == 12) || (abs(pid) == 13) ||
				        (abs(pid) == 14)) continue;
				    if (nprod == 0) fprintf(stream,
				        "%5d %2d %4d %12.5lE %12.3lf %12.3lf "
				        "%12.3lf\n", eventid, generation,
				        product.pid, tau.kinetic,
				        tau.position[0], tau.position[1],
				        tau.position[2]);
				    fprintf(stream, "    %4d %12.5lE %12.5lE "
				    "%12.5lE\n", pid, momentum[0], momentum[1],
				    momentum[2]);
				    nprod++;
				}
				generation++;
				
				/* Process any additional nu_e~ or nu_tau. */
				if (is_nue) transport(ctx_ent, &nu_e,
				    eventid, generation, stream);
				if (is_nut) transport(ctx_ent, &nu_t,
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
                struct ent_state neutrino = { projectile, energy, 0., 0., 1.,
                        { 0., 0., -EARTH_RADIUS }, { 0., 0., 1. } };
                transport(&ctx_ent, &neutrino, i, 1, stream);
               
        }
        fclose(stream);

        /* Finalise and exit to the OS. */
        gracefully_exit(EXIT_SUCCESS);
}
