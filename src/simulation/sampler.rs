use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::convert::Mode;
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::ValueError;
use pyo3::prelude::*;
use ::std::ffi::{c_int, c_void};


impl danton::Sampler {
    const ALTITUDE_MIN: f64 = -1E+05;
    const ALTITUDE_MAX: f64 = 1E+05;

    pub fn destroy(sampler: &mut *mut danton::Sampler) {
        let sampler = sampler as *mut *mut Self;
        unsafe {
            danton::destroy(sampler as *mut *mut c_void);
        }
    }

    pub fn set(&mut self, mode: Mode, decay: bool, particle: &Particle) -> PyResult<()> {
        if (particle.altitude <= Self::ALTITUDE_MIN) || (particle.altitude >= Self::ALTITUDE_MAX) {
            let why = format!(
                "expected a value in ({}, {}), found {}",
                Self::ALTITUDE_MIN,
                Self::ALTITUDE_MAX,
                particle.altitude,
            );
            let err = Error::new(ValueError)
                .what("altitude")
                .why(&why);
            return Err(err.to_err());
        }
        match mode {
            Mode::Backward => self.set_backward(particle),
            Mode::Forward => self.set_forward(decay, particle),
            Mode::Grammage => self.set_grammage(particle),
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

    fn set_grammage(&mut self, particle: &Particle) -> PyResult<()> {
        self.latitude = particle.latitude;
        self.longitude = particle.longitude;
        self.altitude = [particle.altitude; 2];
        self.azimuth = [particle.azimuth + 180.0; 2];
        self.elevation = [-particle.elevation; 2];
        self.energy = [danton::Primary::ENERGY_MIN, danton::Primary::ENERGY_MAX];
        for i in 0..self.weight.len() {
            self.weight[i] = 0.0;
        }
        to_result(
            unsafe { danton::sampler_update(self as *mut Self) },
            None,
        )
    }
}
