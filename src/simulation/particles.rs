use crate::bindings::danton;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{KeyError, NotImplementedError, ValueError};
use crate::utils::numpy::{Dtype, PyArray, ShapeArg};
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


// ===============================================================================================
//
// Iterators over particles like arrays.
//
// ===============================================================================================

pub struct ParticlesIterator<'a> {
    pid: Option<&'a PyArray<i32>>,
    energy: Option<&'a PyArray<f64>>,
    latitude: &'a PyArray<f64>,
    longitude: &'a PyArray<f64>,
    altitude: &'a PyArray<f64>,
    azimuth: &'a PyArray<f64>,
    elevation: &'a PyArray<f64>,
    size: usize,
    index: usize,
}

impl<'a> ParticlesIterator<'a> {
    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        let pid = extract(elements, "pid")?;
        let energy = extract(elements, "energy")?;
        let latitude = extract(elements, "latitude")?;
        let longitude = extract(elements, "longitude")?;
        let altitude = extract(elements, "altitude")?;
        let azimuth = extract(elements, "azimuth")?;
        let elevation = extract(elements, "elevation")?;
        let size = pid.size();
        let others = [
            energy.size(), latitude.size(), longitude.size(), altitude.size(), azimuth.size(),
            elevation.size()
        ];
        if others.iter().any(|x| *x != size) {
            let err = Error::new(ValueError)
                .what("particles")
                .why("differing arrays sizes")
                .to_err();
            return Err(err);
        }
        let index = 0;
        let pid = Some(pid);
        let energy = Some(energy);
        let iter = Self {
            pid, energy, latitude, longitude, altitude, azimuth, elevation, size, index
        };
        Ok(iter)
    }

    pub fn coordinates<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        let latitude = extract(elements, "latitude")?;
        let longitude = extract(elements, "longitude")?;
        let altitude = extract(elements, "altitude")?;
        let azimuth = extract(elements, "azimuth")?;
        let elevation = extract(elements, "elevation")?;
        let size = latitude.size();
        let others = [
            longitude.size(), altitude.size(), azimuth.size(), elevation.size()
        ];
        if others.iter().any(|x| *x != size) {
            let err = Error::new(ValueError)
                .what("coordinates")
                .why("differing arrays sizes")
                .to_err();
            return Err(err);
        }
        let index = 0;
        let pid = None;
        let energy = None;
        let iter = Self {
            pid, energy, latitude, longitude, altitude, azimuth, elevation, size, index
        };
        Ok(iter)
    }

    fn get(&self, index: usize) -> PyResult<Particle> {
        let (pid, energy) = match self.pid {
            None => (0, 0.0),
            Some(pid) => (pid.get(index)?, self.energy.unwrap().get(index)?),
        };
        let particle = Particle {
            pid,
            energy,
            latitude: self.latitude.get(index)?,
            longitude: self.longitude.get(index)?,
            altitude: self.altitude.get(index)?,
            azimuth: self.azimuth.get(index)?,
            elevation: self.elevation.get(index)?,
        };
        Ok(particle)
    }
}

fn extract<'a, 'py, T>(elements: &'a Bound<'py, PyAny>, key: &str) -> PyResult<&'a PyArray<T>>
where
    'py: 'a,
    T: Dtype,
{
    let py = elements.py();
    elements.get_item(key)
        .map_err(|err| {
            Error::new(KeyError)
                .what("particles")
                .why(&err.value_bound(py).to_string()).to_err()
        })?
        .extract()
}

impl<'a> Iterator for ParticlesIterator<'a> {
    type Item = PyResult<Particle>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.index < self.size {
            let item = Some(self.get(self.index));
            self.index += 1;
            item
        } else {
            None
        }
    }
}
