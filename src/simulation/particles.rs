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
