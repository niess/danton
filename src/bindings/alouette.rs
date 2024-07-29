#![allow(unused)]

use ::std::ffi::c_uint;


pub const SUCCESS: c_uint = 0;
pub const VALUE_ERROR: c_uint = 1;
pub const FLOATING_ERROR: c_uint = 2;
pub const TAUOLA_ERROR: c_uint = 3;

#[link(name = "danton-c")]
extern "C" {
    #[link_name="alouette_initialise"]
    pub fn initialise(xk0dec: *mut f64) -> c_uint;
}
