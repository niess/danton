use crate::utils::convert::Array;
use crate::utils::coordinates::{GeodeticCoordinates, HorizontalCoordinates};
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{TypeError, ValueError};
use crate::utils::numpy::{PyArray, PyUntypedArray};
use pyo3::prelude::*;


// ===============================================================================================
//
// Generic geodetic direction.
//
// ===============================================================================================

pub struct Direction<'a> {
    pub azimuth: &'a PyArray<f64>,
    pub elevation: &'a PyArray<f64>,
    size: Size,
}

impl<'a> Direction<'a> {
    pub fn maybe_new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
        let py = elements.py();
        let azimuth = extract(py, elements, "azimuth")?;
        let elevation = extract(py, elements, "elevation")?;

        let all = azimuth.is_some() && elevation.is_some();
        let any = azimuth.is_some() || elevation.is_some();
        if all {
            let azimuth = azimuth.unwrap();
            let elevation = elevation.unwrap();
            let azimuth_size = Size::new(azimuth);
            let elevation_size = Size::new(elevation);
            let size = azimuth_size.common(&elevation_size)
                .ok_or_else(|| Error::new(ValueError)
                    .what("azimuth and elevation")
                    .why("inconsistent arrays sizes")
                    .to_err()
                )?
                .clone();
            let direction = Direction { azimuth, elevation, size };
            Ok(Some(direction))
        } else if any {
            let why = if azimuth.is_some() {
                "missing 'elevation'"
            } else {
                "missing 'azimuth'"
            };
            let err = Error::new(TypeError).what("direction").why(why);
            Err(err.to_err())
        } else {
            Ok(None)
        }
    }

    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        Self::maybe_new(elements)?
            .ok_or_else(|| {
                Error::new(ValueError)
                    .what("direction")
                    .why("missing azimuth and elevation")
                    .to_err()
            })
    }

    pub fn get(&self, i: usize) -> PyResult<HorizontalCoordinates> {
        let horizontal = HorizontalCoordinates {
            azimuth: self.azimuth.get(i)?,
            elevation: self.elevation.get(i)?,
        };
        Ok(horizontal)
    }

    pub fn shape(&self) -> Vec<usize> {
        match &self.size {
            Size::Scalar => Vec::new(),
            Size::Array { shape, .. } => shape.clone(),
        }
    }

    pub fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    pub fn size(&self) -> usize {
        match &self.size {
            Size::Scalar => 1,
            Size::Array { size, .. } => *size,
        }
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
    size: Size,
}

impl<'a> Position<'a> {
    pub fn maybe_new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
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
            let latitude_size = Size::new(latitude);
            let longitude_size = Size::new(longitude);
            let altitude_size = Size::new(altitude);
            let size = latitude_size
                .common(&longitude_size)
                .and_then(|size| size.common(&altitude_size))
                .ok_or_else(|| Error::new(ValueError)
                    .what("latitude, longitude and altitude")
                    .why("inconsistent arrays sizes")
                    .to_err()
                )?
                .clone();
            let position = Self { latitude, longitude, altitude, size };
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
            let err = Error::new(TypeError).what("position").why(&why);
            Err(err.to_err())
        } else {
            Ok(None)
        }
    }

    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        Self::maybe_new(elements)?
            .ok_or_else(|| {
                Error::new(ValueError)
                    .what("direction")
                    .why("missing latitude, longitude and altitude")
                    .to_err()
            })
    }

    pub fn common(&self, direction: &Direction) -> PyResult<(usize, Vec<usize>)> {
        let size = self.size.common(&direction.size)
            .ok_or_else(|| Error::new(ValueError)
                .what("position and direction")
                .why("inconsistent arrays size")
            )?;
        let result = match size {
            Size::Scalar => (1, Vec::new()),
            Size::Array { size, shape } => (*size, shape.clone()),
        };
        Ok(result)
    }

    pub fn get(&self, i: usize) -> PyResult<GeodeticCoordinates> {
        let geodetic = GeodeticCoordinates {
            latitude: self.latitude.get(i)?,
            longitude: self.longitude.get(i)?,
            altitude: self.altitude.get(i)?,
        };
        Ok(geodetic)
    }

    pub fn shape(&self) -> Vec<usize> {
        match &self.size {
            Size::Scalar => Vec::new(),
            Size::Array { shape, .. } => shape.clone(),
        }
    }

    pub fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    pub fn size(&self) -> usize {
        match &self.size {
            Size::Scalar => 1,
            Size::Array { size, .. } => *size,
        }
    }
}


// ===============================================================================================
//
// Generic geodetic projection.
//
// ===============================================================================================

pub struct Projection<'a> {
    pub latitude: &'a PyArray<f64>,
    pub longitude: &'a PyArray<f64>,
    size: Size,
}

impl<'a> Projection<'a> {
    pub fn maybe_new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
        let py = elements.py();
        let latitude = extract(py, elements, "latitude")?;
        let longitude = extract(py, elements, "longitude")?;

        let all = latitude.is_some() && longitude.is_some();
        let any = latitude.is_some() || longitude.is_some();
        if all {
            let latitude = latitude.unwrap();
            let longitude = longitude.unwrap();
            let latitude_size = Size::new(latitude);
            let longitude_size = Size::new(longitude);
            let size = latitude_size.common(&longitude_size)
                .ok_or_else(|| Error::new(ValueError)
                    .what("latitude and longitude")
                    .why("inconsistent arrays sizes")
                    .to_err()
                )?
                .clone();
            let projection = Projection { latitude, longitude, size };
            Ok(Some(projection))
        } else if any {
            let why = if latitude.is_some() {
                "missing 'longitude'"
            } else {
                "missing 'latitude'"
            };
            let err = Error::new(TypeError).what("projection").why(why);
            Err(err.to_err())
        } else {
            Ok(None)
        }
    }

    pub fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Self> {
        Self::maybe_new(elements)?
            .ok_or_else(|| {
                Error::new(ValueError)
                    .what("projection")
                    .why("missing latitude and longitude")
                    .to_err()
            })
    }

    pub fn shape(&self) -> Vec<usize> {
        match &self.size {
            Size::Scalar => Vec::new(),
            Size::Array { shape, .. } => shape.clone(),
        }
    }

    pub fn size(&self) -> usize {
        match &self.size {
            Size::Scalar => 1,
            Size::Array { size, .. } => *size,
        }
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
        })?
        .map(|array: Array<f64>| array.resolve());
    Ok(value)
}


// ===============================================================================================
//
// Managed array size.
//
// ===============================================================================================

#[derive(Clone)]
enum Size {
    Scalar,
    Array { size: usize, shape: Vec<usize> },
}

impl Size {
    fn new(array: &PyUntypedArray) -> Self {
        if array.ndim() == 0 {
            Self::Scalar
        } else {
            Self::Array { size: array.size(), shape: array.shape() }
        }
    }

    fn common<'a>(&'a self, other: &'a Self) -> Option<&'a Self> {
        match self {
            Self::Scalar => Some(other),
            Self::Array { size, .. } => match other {
                Self::Scalar => Some(self),
                Self::Array { size: other_size, .. } => if size == other_size {
                    Some(self)
                } else {
                    None
                }
            }
        }
    }
}
