use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{TypeError, ValueError};
use crate::utils::numpy::PyArray;
use pyo3::prelude::*;


// ===============================================================================================
//
// Generic geodetic direction.
//
// ===============================================================================================

pub struct Direction<'a> {
    pub azimuth: &'a PyArray<f64>,
    pub elevation: &'a PyArray<f64>,
}

impl<'a> Direction<'a> {
    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
        let py = elements.py();
        let azimuth = extract(py, elements, "azimuth")?;
        let elevation = extract(py, elements, "elevation")?;

        let all = azimuth.is_some() && elevation.is_some();
        let any = azimuth.is_some() || elevation.is_some();
        if all {
            let azimuth = azimuth.unwrap();
            let elevation = elevation.unwrap();
            if azimuth.size() != elevation.size() {
                let err = Error::new(ValueError)
                    .what("direction")
                    .why("differing arrays sizes")
                    .to_err();
                return Err(err);
            }
            let direction = Direction { azimuth, elevation };
            Ok(Some(direction))
        } else if any {
            let why = if azimuth.is_some() {
                "missing 'elevation'"
            } else {
                "missing 'azimuth'"
            };
            let err = Error::new(ValueError).what("direction").why(why);
            Err(err.to_err())
        } else {
            Ok(None)
        }
    }

    pub fn shape(&self) -> Vec<usize> {
        self.azimuth.shape()
    }

    pub fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    pub fn size(&self) -> usize {
        self.azimuth.size()
    }
}


// ===============================================================================================
//
// Generic geodetic position.
//
// ===============================================================================================

pub struct Position<'a> {
    pub latitude: &'a PyArray<f64>,
    pub longitude: &'a PyArray<f64>,
    pub altitude: &'a PyArray<f64>,
}

impl<'a> Position<'a> {
    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
        let py = elements.py();
        let latitude = extract(py, elements, "latitude")?;
        let longitude = extract(py, elements, "longitude")?;
        let altitude = extract(py, elements, "altitude")?;

        let all = latitude.is_some() && longitude.is_some() && altitude.is_some();
        let any = latitude.is_some() || longitude.is_some() || altitude.is_some();
        if all {
            let latitude = latitude.unwrap();
            let longitude = longitude.unwrap();
            let altitude = altitude.unwrap();
            let size = latitude.size();
            if (longitude.size() != size) || (altitude.size() != size) {
                let err = Error::new(ValueError)
                    .what("position")
                    .why("differing arrays sizes")
                    .to_err();
                return Err(err);
            }
            let position = Self { latitude, longitude, altitude };
            Ok(Some(position))
        } else if any {
            let mut missing: Vec<&str> = Vec::new();
            if latitude.is_none() { missing.push("latitude") };
            if longitude.is_none() { missing.push("longitude") };
            if altitude.is_none() { missing.push("altitude") };
            let why = if missing.len() == 2 {
                format!("missing '{}' and '{}'", missing[0], missing[1])
            } else {
                format!("missing '{}'", missing[0])
            };
            let err = Error::new(ValueError).what("position").why(&why);
            Err(err.to_err())
        } else {
            Ok(None)
        }
    }

    pub fn shape(&self) -> Vec<usize> {
        self.latitude.shape()
    }

    pub fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    pub fn size(&self) -> usize {
        self.latitude.size()
    }
}


// ===============================================================================================
//
// Generic extraction.
//
// ===============================================================================================

fn extract<'py>(
    py: Python<'py>,
    elements: &Bound<'py, PyAny>,
    key: &str
) -> PyResult<Option<&'py PyArray<f64>>> {
    let value: Option<&PyArray<f64>> = elements
        .get_item(key).ok()
        .and_then(|a| Some(a.extract())).transpose()
        .map_err(|err| {
            Error::new(TypeError)
                .what(key)
                .why(&err.value_bound(py).to_string()).to_err()
        })?;
    Ok(value)
}
