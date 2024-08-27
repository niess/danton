use crate::bindings::danton;
use crate::utils::convert::{Ellipsoid, Geoid, Medium};
use crate::utils::coordinates::{Coordinates, CoordinatesExport, GeodeticCoordinates,
    GeodeticsExport, HorizontalCoordinates};
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::{NotImplementedError, TypeError, ValueError};
use crate::utils::export::Export;
use crate::utils::extract::{Direction, Position};
use crate::utils::numpy::PyArray;
use crate::utils::tuple::NamedTuple;
use pyo3::prelude::*;
use pyo3::types::PyDict;
use ::std::ffi::{c_int, CString, c_void};
use ::std::ptr::null;
use ::std::sync::atomic::{AtomicUsize, Ordering};


// XXX Forward geoid undulations.

static INSTANCES: AtomicUsize = AtomicUsize::new(0);
static CURRENT: AtomicUsize = AtomicUsize::new(0);

#[pyclass(module="danton")]
pub struct Geometry {
    #[pyo3(get)]
    /// Reference geoid for the sea level.
    pub geoid: Geoid,
    #[pyo3(get)]
    /// Topography elevation data.
    topography: Option<Topography>,
    #[pyo3(get)]
    /// The topography composition.
    pub material: String,
    #[pyo3(get)]
    /// The topography density.
    density: f64,
    #[pyo3(get)]
    /// Flag to enable/disable the ocean.
    pub ocean: bool,

    instance: usize,
    modified: bool,
}

#[derive(Clone, FromPyObject, PartialEq)]
enum Topography {
    #[pyo3(transparent, annotation = "str")]
    Dem(String),
    #[pyo3(transparent, annotation = "float")]
    Flat(f64),
}

impl IntoPy<PyObject> for Topography {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::Dem(path) => path.into_py(py),
            Self::Flat(z) => z.into_py(py),
        }
    }
}

#[pymethods]
impl Geometry {
    #[new]
    pub fn new() -> Self {
        Self {
            geoid: Geoid::default(),
            topography: None,
            material: "Rock".to_string(),
            density: 2.65E+03,
            ocean: true,
            instance: INSTANCES.fetch_add(1, Ordering::SeqCst),
            modified: true,
        }
    }

    #[setter]
    fn set_density(&mut self, value: f64) {
        if value != self.density {
            self.modified = true;
            self.density = value;
        }
    }

    /// Reference ellipsoid for geodetic coordinates.
    #[getter]
    fn get_ellipsoid(&self) -> Ellipsoid {
        self.geoid.into()
    }

    #[setter]
    fn set_geoid(&mut self, value: Geoid) {
        if value != self.geoid {
            self.modified = true;
            self.geoid = value;
        }
    }

    #[setter]
    fn set_material(&mut self, value: String) {
        if value != self.material {
            self.modified = true;
            self.material = value;
        }
    }

    #[setter]
    fn set_ocean(&mut self, value: bool) {
        if value != self.ocean {
            self.modified = true;
            self.ocean = value;
        }
    }

    #[setter]
    fn set_topography(&mut self, value: Option<Topography>) {
        if value != self.topography {
            self.modified = true;
            self.topography = value;
        }
    }

    /// Convert ECEF coordinates to geodetic ones.
    fn from_ecef(
        &self,
        py: Python,
        position: &PyArray<f64>,
        direction: Option<&PyArray<f64>>
    ) -> PyResult<PyObject> {
        let size = position.size();
        if let Some(direction) = direction {
            if direction.size() != size {
                let why = format!(
                    "expected a size {} array, found a size {} array",
                    size,
                    direction.size(),
                );
                let err = Error::new(ValueError)
                    .what("direction")
                    .why(&why);
                return Err(err.to_err())
            }
        }

        fn check_shape(array: &PyArray<f64>, what: &str) -> PyResult<()> {
            let shape = array.shape();
            if shape.len() == 1 && shape[0] == 3 {
                return Ok(())
            }
            if (shape.len() < 2) || (shape[shape.len() - 1] != 3) {
                let why = format!(
                    "expected a shape [.., 3] array, found a shape {:?} array",
                    shape,
                );
                let err = Error::new(ValueError)
                    .what(what)
                    .why(&why);
                Err(err.to_err())
            } else {
                Ok(())
            }
        }

        check_shape(position, "position")?;
        if let Some(direction) = direction {
            check_shape(direction, "direction")?
        }

        let (mut coordinates, mut geodetics) = match direction {
            Some(_) => {
                let coordinates: Vec<Coordinates> = Vec::new();
                (Some(coordinates), None)
            },
            None => {
                let geodetics: Vec<GeodeticCoordinates> = Vec::new();
                (None, Some(geodetics))
            },
        };

        for i in 0..(size / 3) {
            let geodetic = {
                let position = [
                    position.get(3 * i)?,
                    position.get(3 * i + 1)?,
                    position.get(3 * i + 2)?,
                ];
                GeodeticCoordinates::from_ecef(&position, self.geoid.into())
            };
            match direction {
                Some(direction) => {
                    let direction = [
                        direction.get(3 * i)?,
                        direction.get(3 * i + 1)?,
                        direction.get(3 * i + 2)?,
                    ];
                    let horizontal = HorizontalCoordinates::from_ecef(
                        &direction, self.geoid.into(), &geodetic
                    );
                    let coordinate = Coordinates {
                        position: geodetic,
                        direction: horizontal,
                    };
                    coordinates.as_mut().map(|v| v.push(coordinate));
                },
                None => {
                    geodetics.as_mut().map(|v| v.push(geodetic));
                },
            }
        }

        let result = match direction {
            Some(_) => Export::export::<CoordinatesExport>(py, coordinates.unwrap())?,
            None => Export::export::<GeodeticsExport>(py, geodetics.unwrap())?,
        };
        Ok(result)
    }

    /// Convert geodetic coordinates to ECEF ones.
    #[pyo3(signature=(elements=None, **kwargs))]
    fn to_ecef<'py>(&self,
        py: Python<'py>,
        elements: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let (position, direction) = match elements {
            None => match kwargs {
                None => return Ok(py.None()),
                Some(kwargs) => {
                    self.to_ecef_from_kwargs(kwargs)?
                },
            },
            Some(elements) => {
                if kwargs.is_some() {
                    let err = Error::new(NotImplementedError)
                        .what("arguments")
                        .why("cannot mix positional and keyword arguments");
                    return Err(err.to_err())
                }
                self.to_ecef_from_elements(elements)?
            }
        };

        let result = match direction {
            None => position,
            Some(direction) => {
                static RESULT: NamedTuple<2> = NamedTuple::new(
                    "Result", ["position", "direction"]);
                RESULT.instance(py, (position, direction))?.unbind()
            }
        };
        Ok(result)
    }
}

impl Geometry {
    pub fn apply(&mut self) -> PyResult<()> {
        let previous = CURRENT.swap(self.instance, Ordering::SeqCst);
        if (previous != self.instance) || self.modified {
            let geoid: &str = self.geoid.into();
            let geoid = CString::new(geoid)?;
            let topography = self.topography.as_ref().map(|topography| match topography {
                Topography::Dem(path) => CString::new(path.as_str()),
                Topography::Flat(z) => CString::new(format!("flat://{:.6}", z)),
            })
                .transpose()?;
            let topography = match topography {
                None => null(),
                Some(cstr) => cstr.as_ptr(),
            };
            let mut ocean = if self.ocean { 1 } else { 0 };
            to_result(
                unsafe {
                    danton::earth_model(
                        geoid.as_ptr(),
                        topography,
                        self.density,
                        &mut ocean,
                    )
                },
                None
            )?;
            self.modified = false;
        }
        Ok(())
    }

    fn to_ecef_from_elements<'py>(
        &self,
        elements: &Bound<'py, PyAny>
    ) -> PyResult<(PyObject, Option<PyObject>)> {
        let py = elements.py();
        let position = Position::new(elements)?
            .ok_or_else(|| {
                Error::new(ValueError)
                    .what("position")
                    .why("missing latitude, longitude and altitude")
                    .to_err()
            })?;
        let direction = Direction::new(elements)?;

        let size = match direction.as_ref() {
            None => position.size(),
            Some(direction) => {
                let size = position.size();
                if direction.size() != size {
                    let err = Error::new(ValueError)
                        .what("direction")
                        .why("differing arrays sizes")
                        .to_err();
                    return Err(err);
                }
                size
            },
        };

        let ecef_position = PyArray::<f64>::empty(py, &position.shape3())?;
        let ecef_direction: Option<&PyArray<f64>> = match direction.as_ref() {
            None => None,
            Some(direction) => Some(PyArray::<f64>::empty(py, &direction.shape3())?),
        };

        for i in 0..size {
            let geodetic = GeodeticCoordinates {
                latitude: position.latitude.get(i)?,
                longitude: position.longitude.get(i)?,
                altitude: position.altitude.get(i)?,
            };
            let r = geodetic.to_ecef(self.geoid.into());
            ecef_position.set(3 * i, r[0])?;
            ecef_position.set(3 * i + 1, r[1])?;
            ecef_position.set(3 * i + 2, r[2])?;
            if let Some(Direction { azimuth, elevation }) = direction {
                let horizontal = HorizontalCoordinates {
                    azimuth: azimuth.get(i)?,
                    elevation: elevation.get(i)?,
                };
                let u = horizontal.to_ecef(self.geoid.into(), &geodetic);
                let ecef_direction = ecef_direction.unwrap();
                ecef_direction.set(3 * i, u[0])?;
                ecef_direction.set(3 * i + 1, u[1])?;
                ecef_direction.set(3 * i + 2, u[2])?;
            }
        }

        let ecef_position: &PyAny = ecef_position;
        let ecef_position = ecef_position.into_py(py);

        let result = {
            match ecef_direction {
                None => (ecef_position, None),
                Some(ecef_direction) => {
                    let ecef_direction: &PyAny = ecef_direction;
                    let ecef_direction = ecef_direction.into_py(py);
                    (ecef_position, Some(ecef_direction))
                },
            }
        };

        Ok(result)
    }

    fn to_ecef_from_kwargs<'py>(
        &self,
        kwargs: &Bound<'py, PyDict>
    ) -> PyResult<(PyObject, Option<PyObject>)> {
        let py = kwargs.py();
        let mut position: Option<GeodeticCoordinates> = None;
        let mut direction: Option<HorizontalCoordinates> = None;

        let default_position = || GeodeticCoordinates {
            latitude: 0.0, longitude: 0.0, altitude: 0.0
        };
        let default_direction = || HorizontalCoordinates {
            azimuth: 0.0, elevation: 0.0
        };

        for (key, value) in kwargs.iter() {
            let key: String = key.extract()?;
            let value: f64 = value.extract()?;
            match key.as_str() {
                    "latitude" => {
                        position
                            .get_or_insert_with(default_position)
                            .latitude = value;
                    },
                    "longitude" => {
                        position
                            .get_or_insert_with(default_position)
                            .longitude = value;
                    },
                    "altitude" => {
                        position
                            .get_or_insert_with(default_position)
                            .altitude = value;
                    },
                    "azimuth" => {
                        direction
                            .get_or_insert_with(default_direction)
                            .azimuth = value;
                    },
                    "elevation" => {
                        direction
                            .get_or_insert_with(default_direction)
                            .elevation = value;
                    },
                    _ => {
                        let why = format!("unexpected keyword argument '{}'", key);
                        let err = Error::new(TypeError)
                            .what("argument")
                            .why(&why);
                        return Err(err.to_err())
                    }
            }
        }

        let position = position
            .ok_or_else(|| {
                Error::new(ValueError)
                    .what("position")
                    .why("missing latitude, longitude and altitude")
                    .to_err()
            })?;
        let ecef_position = PyArray::<f64>::empty(py, &[3])?;
        let r = position.to_ecef(self.geoid.into());
        for i in 0..3 {
            ecef_position.set(i, r[i])?
        }

        let ecef_direction: Option<&PyArray<f64>> = match direction.as_ref() {
            None => None,
            Some(direction) => {
                let ecef_direction = PyArray::<f64>::empty(py, &[3])?;
                let u = direction.to_ecef(self.geoid.into(), &position);
                for i in 0..3 {
                    ecef_direction.set(i, u[i])?
                }
                Some(ecef_direction)
            },
        };

        let ecef_position: &PyAny = ecef_position;
        let ecef_position = ecef_position.into_py(py);

        let result = {
            match ecef_direction {
                None => (ecef_position, None),
                Some(ecef_direction) => {
                    let ecef_direction: &PyAny = ecef_direction;
                    let ecef_direction = ecef_direction.into_py(py);
                    (ecef_position, Some(ecef_direction))
                },
            }
        };

        Ok(result)
    }
}


// ===============================================================================================
//
// Tracer interface.
//
// ===============================================================================================

pub struct Tracer<'a> {
    geometry: &'a Geometry,
    tracer: *mut danton::Tracer,
}

impl<'a> Drop for Tracer<'a> {
    fn drop(&mut self) {
        unsafe {
            danton::destroy(&mut (self.tracer as *mut c_void));
        }
    }
}

impl<'a> Tracer<'a> {
    pub fn new(geometry: &'a mut Geometry) -> PyResult<Self> {
        geometry.apply()?;
        let tracer = unsafe { danton::tracer_create() };
        let tracer = Self { geometry, tracer };
        Ok(tracer)
    }

    pub fn medium(&self, position: &GeodeticCoordinates) -> Medium {
        let r = position.to_ecef(self.geometry.geoid.into());
        let medium = unsafe { danton::tracer_medium(self.tracer, &r as *const f64) };
        (medium, self.geometry.ocean).into()
    }

    pub fn trace(
        &self,
        position: &GeodeticCoordinates,
        direction: &HorizontalCoordinates
    ) -> (Medium, f64, Medium) {
        let r = position.to_ecef(self.geometry.geoid.into());
        let u = direction.to_ecef(self.geometry.geoid.into(), &position);
        let mut distance: f64 = 0.0;
        let mut next_medium: c_int = -1;
        let medium = unsafe {
            danton::tracer_trace(
                self.tracer,
                &r as *const f64,
                &u as *const f64,
                &mut distance,
                &mut next_medium,
            )
        };
        let medium = (medium, self.geometry.ocean).into();
        let next_medium = (next_medium, self.geometry.ocean).into();
        (medium, distance, next_medium)
    }
}
