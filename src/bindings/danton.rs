#![allow(unused)]

use crate::bindings::{ent, pumas};
use crate::simulation::random::RandomContext;
use ::std::ffi::{c_char, c_int, c_long, c_uint, c_ulong, c_void};

pub const SUCCESS: c_int = 0;
pub const FAILURE: c_int = 1;

pub const BACKWARD: c_uint = 0;
pub const FORWARD: c_uint = 1;
pub const GRAMMAGE: c_uint = 2;

pub const UNKNOWN: c_int = -1;
pub const NU_BAR_TAU: c_int = 0;
pub const NU_BAR_MU: c_int = 1;
pub const NU_BAR_E: c_int = 2;
pub const NU_E: c_int = 3;
pub const NU_MU: c_int = 4;
pub const NU_TAU: c_int = 5;
pub const TAU_BAR: c_int = 6;
pub const TAU: c_int = 7;

pub const START: c_uint = 0;
pub const STOP: c_uint = 1;
pub const STEP: c_uint = 2;

pub const PREM_EARTH_RADIUS: f64 = 6371.0E+03;

#[repr(C)]
#[derive(Debug, Copy, Clone)]
pub struct Primary {
    pub flux: PrimaryFlux,
    pub energy: [f64; 2],
}

pub type PrimaryFlux = Option<
    unsafe extern "C" fn(primary: *mut Primary, energy: f64) -> f64,
>;

#[repr(C)]
pub struct Sampler {
    pub latitude: f64,
    pub longitude: f64,
    pub altitude: [f64; 2],
    pub azimuth: [f64; 2],
    pub elevation: [f64; 2],
    pub energy: [f64; 2],
    pub weight: [f64; 8],
}

#[repr(C)]
pub struct State {
    pub pid: c_int,
    pub energy: f64,
    pub position: [f64; 3],
    pub direction: [f64; 3],
}

#[repr(C)]
pub struct Product {
    pub pid: c_int,
    pub momentum: [f64; 3],
}

#[repr(C)]
pub struct Event {
    pub id: c_long,
    pub weight: f64,
    pub primary: *mut State,
    pub generation: c_int,
    pub vertex: *mut State,
    pub secondary: *mut State,
    pub n_products: c_int,
    pub product: *mut Product,
}

#[derive(Clone, Copy, Default)]
#[repr(C)]
pub struct MaterialIndex {
    pub rock: c_int,
    pub water: c_int,
    pub air: c_int,
    pub topography: c_int,
}

#[derive(Clone, Copy)]
#[repr(C)]
pub struct Physics {
    pub ent: *mut ent::Physics,
    pub pumas: *mut pumas::Physics,
}

#[repr(C)]
pub struct Grammage {
    pub elevation: f64,
    pub value: f64,
}

pub type RecordEvent = Option<
    unsafe extern "C" fn(
        context: *mut Context,
        recorder: *mut Recorder,
        event: *const Event,
    ) -> c_int,
>;

pub type RecordGrammage = Option<
    unsafe extern "C" fn(
        context: *mut Context,
        recorder: *mut Recorder,
        grammage: *const Grammage,
    ) -> c_int,
>;

#[repr(C)]
pub struct Recorder {
    pub record_event: RecordEvent,
    pub record_grammage: RecordGrammage,
}

pub type Call = Option<
    unsafe extern "C" fn(
        context: *mut Context,
        run_action: *mut RunAction,
        event: c_uint,
        medium: c_int,
        state: *mut State,
    ) -> c_int,
>;

#[repr(C)]
pub struct RunAction {
    pub call: Call,
}

#[repr(C)]
pub struct Context {
    pub mode: c_uint,
    pub longitudinal: c_int,
    pub decay: c_int,
    pub primary: [*mut Primary; 6],
    pub sampler: *mut Sampler,
    pub recorder: *mut Recorder,
    pub run_action: *mut RunAction,
}

pub type Lock = Option<
    unsafe extern "C" fn() -> c_int,
>;

#[repr(C)]
pub struct Tracer {
    _unused: [u8; 0],
}

pub const MERGE: c_uint = 0;
pub const RESOLVE: c_uint = 1;

#[link(name = "danton-c")]
extern "C" {
    #[link_name="danton_initialise"]
    pub fn initialise(prefix: *const c_char, lock: Lock, unlock: Lock) -> c_int;

    #[link_name="danton_finalise"]
    pub fn finalise();

    #[link_name="danton_destroy"]
    pub fn destroy(any: *mut *mut c_void);

    #[link_name="danton_earth_model"]
    pub fn earth_model(
        geoid: *const c_char,
        topography: *const c_char,
        density: f64,
        ocean: *mut c_int,
    ) -> c_int;

    #[link_name="danton_get_topography"]
    pub fn get_topography() -> *mut c_void;

    #[link_name="danton_particle_pdg"]
    pub fn particle_pdg(particle: c_int) -> c_int;

    #[link_name="danton_particle_index"]
    pub fn particle_index(pdg: c_int) -> c_int;

    #[link_name="danton_sampler_create"]
    pub fn sampler_create() -> *mut Sampler;

    #[link_name="danton_sampler_update"]
    pub fn sampler_update(sampler: *mut Sampler) -> c_int;

    #[link_name="danton_context_create"]
    pub fn context_create() -> *mut Context;

    #[link_name="danton_context_destroy"]
    pub fn context_destroy(context: *mut *mut Context);

    #[link_name="danton_context_run"]
    pub fn context_run(
        context: *mut Context,
        events: c_long,
        requested: c_long,
    ) -> c_int;

    #[link_name="danton_error_count"]
    pub fn error_count(context: *mut Context) -> c_int;

    #[link_name="danton_error_pop"]
    pub fn error_pop(context: *mut Context) -> *const c_char;

    #[link_name="danton_error_push"]
    pub fn error_push(
        context: *mut Context,
        format: *const c_char,
        ...
    ) -> c_int;

    // Rust interface.

    #[allow(improper_ctypes)]
    #[link_name="danton_context_random_set"]
    pub fn context_random_set(context: *mut Context, random: *mut RandomContext);

    #[link_name="danton_context_reset"]
    pub fn context_reset(context: *mut Context);

    #[link_name="danton_materials_set"]
    pub fn materials_set();

    #[link_name="danton_material_index"]
    pub static mut material_index: *mut MaterialIndex;

    #[link_name="danton_physics"]
    pub static mut physics: *mut Physics;

    #[link_name="danton_tau_ctau0"]
    pub static mut tau_ctau0: *mut f64;

    #[link_name="danton_tau_mass"]
    pub static mut tau_mass: *mut f64;

    #[link_name="danton_geoid_undulation"]
    pub fn geoid_undulation(latitude: f64, longitude: f64) -> f64;

    #[link_name="danton_topography_elevation"]
    pub fn topography_elevation(latitude: f64, longitude: f64) -> f64;

    #[link_name="danton_tracer_create"]
    pub fn tracer_create(mode: c_uint) -> *mut Tracer;

    #[link_name="danton_tracer_medium"]
    pub fn tracer_medium(tracer: *mut Tracer, position: *const f64) -> c_int;

    #[link_name="danton_tracer_trace"]
    pub fn tracer_trace(
        tracer: *mut Tracer,
        position: *const f64,
        direction: *const f64,
        distance_max: f64,
        distance: *mut f64,
        next_medium: *mut c_int,
    ) -> c_int;
}
