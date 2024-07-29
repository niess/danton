#![allow(unused)]

use ::std::ffi::{c_char, c_int, c_long, c_uint};

pub const SUCCESS: c_uint = 0;
pub const BAD_ADDRESS: c_uint = 1;
pub const DOMAIN_ERROR: c_uint = 2;
pub const FORMAT_ERROR: c_uint = 3;
pub const IO_ERROR: c_uint = 4;
pub const MEMORY_ERROR: c_uint = 5;
pub const PATH_ERROR: c_uint = 6;

#[repr(C)]
pub struct Physics {
    _unused: [u8; 0],
}

#[link(name = "danton-c")]
extern "C" {
    #[link_name="ent_physics_create"]
    pub fn physics_create(
        physics: *mut *mut Physics,
        pdf: *const c_char,
        cs: *const c_char,
    ) -> c_uint;

    #[link_name="ent_physics_destroy"]
    pub fn physics_destroy(physics: *mut *mut Physics);
}
