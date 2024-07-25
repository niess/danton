use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::convert::Mode;
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::NotImplementedError;
use pyo3::prelude::*;
use ::std::ffi::{c_int, c_void};


impl danton::Sampler {
    const ALTITUDE_MAX: f64 = 1E+05;

    pub fn destroy(sampler: &mut *mut danton::Sampler) {
        let sampler = sampler as *mut *mut Self;
        unsafe {
            danton::destroy(sampler as *mut *mut c_void);
        }
    }

    pub fn set(&mut self, mode: Mode, decay: bool, particle: &Particle) -> PyResult<()> {
        match mode {
            Mode::Backward => self.set_backward(particle),
            Mode::Forward => self.set_forward(decay, particle),
            Mode::Grammage => Err(Error::new(NotImplementedError).to_err()),
        }
    }

    fn set_backward(&mut self, particle: &Particle) -> PyResult<()> {
        self.latitude = particle.latitude;
        self.longitude = particle.longitude;
        self.altitude = [particle.altitude; 2];
        self.azimuth = [particle.azimuth; 2];
        self.elevation = [particle.elevation; 2];
        self.energy = [particle.energy; 2];
        let particle_index = particle.index()?;
        for i in 0..self.weight.len() {
            let p = if (i as c_int) == particle_index { 1.0 } else { 0.0 };
            self.weight[i] = p;
        }
        to_result(
            unsafe { danton::sampler_update(self as *mut Self) },
            None,
        )
    }

    fn set_forward(&mut self, decay: bool, particle: &Particle) -> PyResult<()> {
        self.latitude = particle.latitude;
        self.longitude = particle.longitude;
        self.altitude = if decay {
            [particle.altitude, Self::ALTITUDE_MAX]
        } else {
            [particle.altitude, particle.altitude]
        };
        self.azimuth = [particle.azimuth; 2];
        self.elevation = [particle.elevation; 2];
        self.energy = [danton::Primary::ENERGY_MIN, danton::Primary::ENERGY_MAX];
        for i in 0..self.weight.len() {
            self.weight[i] = 1.0;
        }
        to_result(
            unsafe { danton::sampler_update(self as *mut Self) },
            None,
        )
    }
}
