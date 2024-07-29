#![allow(unused)]

use ::std::ffi::{c_char, c_int, c_long, c_uint};

pub const MUON: c_uint = 0;
pub const TAU: c_uint = 1;

pub const SUCCESS: c_uint = 0;
pub const ACCURACY_ERROR: c_uint = 1;
pub const END_OF_FILE: c_uint = 2;
pub const DECAY_ERROR: c_uint = 3;
pub const DENSITY_ERROR: c_uint = 4;
pub const DIRECTION_ERROR: c_uint = 5;
pub const INCOMPLETE_FILE: c_uint = 6;
pub const INDEX_ERROR: c_uint = 7;
pub const PHYSICS_ERROR: c_uint = 8;
pub const INTERNAL_ERROR: c_uint = 9;
pub const IO_ERROR: c_uint = 10;
pub const FORMAT_ERROR: c_uint = 11;
pub const MEDIUM_ERROR: c_uint = 12;
pub const MEMORY_ERROR: c_uint = 13;
pub const MODEL_ERROR: c_uint = 14;
pub const MISSING_LIMIT: c_uint = 15;
pub const MISSING_RANDOM: c_uint = 16;
pub const PATH_ERROR: c_uint = 17;
pub const RAISE_ERROR: c_uint = 18;
pub const TOO_LONG: c_uint = 19;
pub const UNDEFINED_MDF: c_uint = 20;
pub const UNKNOWN_ELEMENT: c_uint = 21;
pub const UNKNOWN_MATERIAL: c_uint = 22;
pub const UNKNOWN_PARTICLE: c_uint = 23;
pub const VALUE_ERROR: c_uint = 24;

pub const BREMSSTRAHLUNG: c_uint = 0;
pub const PAIR_PRODUCTION: c_uint = 1;
pub const PHOTONUCLEAR: c_uint = 2;

#[repr(C)]
pub struct Physics {
    _unused: [u8; 0],
}

#[repr(C)]
pub struct PhysicsSettings {
    pub cutoff: f64,
    pub elastic_ratio: f64,
    pub bremsstrahlung: *const c_char,
    pub pair_production: *const c_char,
    pub photonuclear: *const c_char,
    pub n_energies: c_int,
    pub energy: *mut f64,
    pub update: c_int,
    pub dry: c_int,
}

pub type Dcs = Option<
    unsafe extern "C" fn(Z: f64, A: f64, m: f64, K: f64, q: f64) -> f64,
>;

#[link(name = "danton-c")]
extern "C" {
    #[link_name="pumas_physics_create"]
    pub fn physics_create(
        physics: *mut *mut Physics,
        particle: c_uint,
        mdf_path: *const c_char,
        dedx_path: *const c_char,
        settings: *const PhysicsSettings,
    ) -> c_uint;

    #[link_name="pumas_physics_dcs"]
    pub fn physics_dcs(
        physics: *const Physics,
        process: c_uint,
        model: *mut *const c_char,
        dcs: *mut Dcs,
    ) -> c_uint;

    #[link_name="pumas_physics_destroy"]
    pub fn physics_destroy(physics: *mut *mut Physics);

    #[link_name="pumas_physics_dump"]
    pub fn physics_dump(physics: *const Physics, stream: *mut libc::FILE) -> c_uint;

    #[link_name="pumas_physics_load"]
    pub fn physics_load(physics: *mut *mut Physics, stream: *mut libc::FILE) -> c_uint;

    #[link_name="pumas_physics_material_index"]
    pub fn physics_material_index(
        physics: *const Physics,
        name: *const c_char,
        index: *mut c_int
    ) -> c_uint;

    #[link_name="pumas_physics_particle"]
    pub fn physics_particle(
        physics: *const Physics,
        particle: *mut c_uint,
        ctau: *mut f64,
        mass: *mut f64
    ) -> c_uint;
}
