use crate::bindings::danton;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::NotImplementedError;
use crate::utils::numpy::{PyArray, ShapeArg};
use pyo3::prelude::*;
use ::std::ffi::c_int;


#[repr(C)]
#[derive(Clone, Copy)]
pub struct Particle {
    pub pid: c_int,
    pub energy: f64,
    pub latitude: f64,
    pub longitude: f64,
    pub altitude: f64,
    pub azimuth: f64,
    pub elevation: f64,
}

impl Particle {
    pub fn index(&self) -> PyResult<c_int> {
        let particle_index = unsafe { danton::particle_index(self.pid) };
        if particle_index == danton::UNKNOWN {
            let why = format!("{}", self.pid);
            let err = Error::new(NotImplementedError)
                .what("pid")
                .why(&why);
            Err(err.to_err())
        } else {
            Ok(particle_index)
        }
    }
}

/// Create an array of Monte Carlo particles.
#[pyfunction]
pub fn particles(
    py: Python,
    shape: ShapeArg,
) -> PyResult<PyObject> {
    let shape: Vec<usize> = shape.into();
    let array: &PyAny = PyArray::<Particle>::zeros(py, &shape)?;
    Ok(array.into())
}
