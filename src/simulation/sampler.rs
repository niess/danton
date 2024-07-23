use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::NotImplementedError;
use pyo3::prelude::*;
use ::std::ffi::{c_int, c_void};


impl danton::Sampler {
    pub fn destroy(sampler: &mut *mut danton::Sampler) {
        let sampler = sampler as *mut *mut Self;
        unsafe {
            danton::destroy(sampler as *mut *mut c_void);
        }
    }

    pub fn set(&mut self, particle: &Particle) -> PyResult<()> {
        self.latitude = particle.latitude;
        self.longitude = particle.longitude;
        self.altitude = [particle.altitude; 2];
        self.azimuth = [particle.azimuth; 2];
        self.elevation = [particle.elevation; 2];
        self.energy = [particle.energy; 2];
        let particle_index = unsafe { danton::particle_index(particle.pid) };
        if particle_index == danton::UNKNOWN {
            let why = format!("{}", particle.pid);
            let err = Error::new(NotImplementedError)
                .what("pid")
                .why(&why);
            return Err(err.to_err());
        }
        for i in 0..self.weight.len() {
            let p = if (i as c_int) == particle_index { 1.0 } else { 0.0 };
            self.weight[i] = p;
        }
        to_result(
            unsafe { danton::sampler_update(self as *mut Self) },
            None,
        )
    }
}
