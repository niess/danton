use crate::bindings::danton;
use crate::utils::convert::{Array, Ellipsoid, Geoid, Medium, Reference};
use crate::utils::coordinates::{Coordinates, CoordinatesExport, GeodeticCoordinates,
    GeodeticsExport, HorizontalCoordinates};
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::{NotImplementedError, ValueError};
use crate::utils::export::Export;
use crate::utils::extract::{Direction, Distance, Position, Projection, select_coordinates,
    select_position, select_projection};
use crate::utils::float::f64x3;
use crate::utils::numpy::PyArray;
use crate::utils::tuple::NamedTuple;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString};
use ::std::ffi::{c_int, CString, c_uint, c_void};
use ::std::ptr::null;
use ::std::sync::atomic::{AtomicUsize, Ordering};


static INSTANCES: AtomicUsize = AtomicUsize::new(0);
static CURRENT: AtomicUsize = AtomicUsize::new(0);

#[pyclass(module="danton")]
pub struct Geometry {
    #[pyo3(get)]
    /// Reference geoid for the sea level.
    pub geoid: Geoid,
    #[pyo3(get)]
    /// The topography elevation model.
    topography: Option<Topography>,
    #[pyo3(get)]
    /// The topography composition.
    pub material: String,
    #[pyo3(get)]
    /// The topography density.
    density: f64,
    #[pyo3(get)]
    /// Flag to enable / disable the ocean.
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
    #[pyo3(signature=(**kwargs))]
    #[new]
    pub fn new<'py>(
        py: Python<'py>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Py<Self>> {
        let geometry = Self {
            geoid: Geoid::default(),
            topography: None,
            material: "Rock".to_string(),
            density: 2.65E+03,
            ocean: true,
            instance: INSTANCES.fetch_add(1, Ordering::SeqCst),
            modified: true,
        };
        let geometry = Bound::new(py, geometry)?;

        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key: Bound<PyString> = key.extract()?;
                geometry.setattr(key, value)?
            }
        }

        Ok(geometry.unbind())
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
        position: Array<f64>,
        direction: Option<Array<f64>>
    ) -> PyResult<PyObject> {
        let position = position.resolve();
        let size = position.size();
        let direction = match direction {
            Some(direction) => {
                let direction = direction.resolve();
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
                Some(direction)
            },
            None => None,
        };

        fn check_shape(array: &PyArray<f64>, what: &str) -> PyResult<usize> {
            let shape = array.shape();
            if shape.len() == 1 && shape[0] == 3 {
                return Ok(1)
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
                Ok(shape.len())
            }
        }

        let ndim = check_shape(position, "position")?;
        if let Some(direction) = direction {
            check_shape(direction, "direction")?;
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
        let result = if ndim == 1 {
            result
                .bind(py)
                .get_item(0)?
                .unbind()
        } else {
            result
        };
        Ok(result)
    }

    /// Get undulation(s) w.r.t. the geoid.
    #[pyo3(signature=(array=None, /, **kwargs))]
    fn geoid_undulation<'py>(
        &mut self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        match select_projection(array, kwargs)? {
            None => return Ok(py.None()),
            Some(any) => {
                self.apply()?;
                self.geoid_undulation_from_any(any)
            },
        }
    }

    /// Return medium at coordinates.
    #[pyo3(signature=(array=None, /, **kwargs))]
    fn locate<'py>(
        &mut self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let position = match select_position(array, kwargs)? {
            Some(any) => Position::new(any)?,
            None => return Ok(py.None()),
        };
        let size = position.size();
        let media = PyArray::<[u8; 16]>::empty(py, &position.shape())?;

        let tracer = Tracer::new(self, Mode::Full)?;
        for i in 0..size {
            let geodetic = position.get(i)?;
            let medium = tracer.medium(&geodetic);
            media.set(i, medium.into())?;
        }

        let media: &PyAny = media;
        Ok(media.into_py(py))
    }

    /// Convert geodetic coordinates to ECEF ones.
    #[pyo3(signature=(array=None, /, **kwargs))]
    fn to_ecef<'py>(
        &self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let (position, direction) = match select_coordinates(array, kwargs)? {
            None => return Ok(py.None()),
            Some(any) => self.to_ecef_from_any(any)?,
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

    /// Get topography elevation value(s).
    #[pyo3(signature=(array=None, /, *, reference=None, **kwargs))]
    fn topography_elevation<'py>(
        &mut self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        reference: Option<Reference>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let reference = reference.unwrap_or_else(|| Reference::default());
        match select_projection(array, kwargs)? {
            None => Ok(py.None()),
            Some(any) => {
                self.apply()?;
                self.topography_elevation_from_any(any, reference)
            },
        }
    }

    /// Trace the distance to the next medium.
    #[pyo3(signature=(array=None, /, *, backward=None, full=None, **kwargs))]
    fn trace<'py>(
        &mut self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        backward: Option<bool>,
        full: Option<bool>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let backward = backward.unwrap_or(false);
        let full = full.unwrap_or(false);
        let (position, direction) = match select_coordinates(array, kwargs)? {
            Some(any) => (Position::new(any)?, Direction::new(any)?),
            None => return Ok(py.None()),
        };
        let (size, shape) = position.common(&direction)?;
        let traces = PyArray::<Trace>::empty(py, &shape)?;

        let mode = if full { Mode::Full } else { Mode:: Merge };
        let tracer = Tracer::new(self, mode)?;
        for i in 0..size {
            let geodetic = position.get(i)?;
            let mut horizontal = direction.get(i)?;
            if backward {
                horizontal.azimuth += 180.0;
                horizontal.elevation = -horizontal.elevation;
            }
            let (current, distance, next) = tracer.trace(&geodetic, &horizontal);
            let trace = Trace { distance, current: current.into(), next: next.into() };
            traces.set(i, trace)?;
        }

        let traces: &PyAny = traces;
        Ok(traces.into_py(py))
    }

    /// Translate geodetic coordinates.
    #[pyo3(signature=(distance, array=None, /, *, copy=None, **kwargs))]
    fn translate<'py>(
        &mut self,
        py: Python<'py>,
        distance: &Bound<'py, PyAny>,
        array: Option<&Bound<'py, PyAny>>,
        copy: Option<bool>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let copy = copy.unwrap_or(true);
        let distance = Distance::new(distance)?;
        let (position, direction) = match select_coordinates(array, kwargs)? {
            Some(any) => (Position::new(any)?, Direction::new(any)?),
            None => return Ok(py.None()),
        };
        if kwargs.is_some() && !copy {
            let err = Error::new(NotImplementedError)
                .what("arguments")
                .why("cannot in-place translate kwargs");
            return Err(err.to_err())
        }
        let (size, shape) = distance.common(&position, &direction)?;
        let result = if copy {
            Some(PyArray::<Coordinates>::empty(py, &shape)?)
        } else {
            None
        };

        for i in 0..size {
            let distance = distance.get(i)?;
            let geodetic = position.get(i)?;
            let horizontal = direction.get(i)?;
            let mut r: f64x3 = (&geodetic.to_ecef(self.geoid.into())).into();
            let u: f64x3 = (&horizontal.to_ecef(self.geoid.into(), &geodetic)).into();
            r += distance * u;
            let geodetic = GeodeticCoordinates::from_ecef(&r.into(), self.geoid.into());
            let horizontal = HorizontalCoordinates::from_ecef(
                &u.into(), self.geoid.into(), &geodetic
            );
            match result {
                Some(result) => {
                    let coordinates = Coordinates { position: geodetic, direction: horizontal };
                    result.set(i, coordinates)?;
                },
                None => {
                    position.set(i, &geodetic)?;
                    direction.set(i, &horizontal)?;
                },
            }
        }

        let result = match result {
            Some(result) => {
                let result: &PyAny = result;
                result.into_py(py)
            },
            None => py.None(),
        };
        Ok(result)
    }
}

#[allow(dead_code)]
#[derive(Clone, Copy)]
#[repr(C)]
pub struct Trace {
    distance: f64,
    current: [u8; 16],
    next: [u8; 16],
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
            let mut ocean = if self.ocean { 1 } else { 0 };
            to_result(
                unsafe {
                    match topography {
                        None => danton::earth_model(
                            geoid.as_ptr(),
                            null(),
                            self.density,
                            &mut ocean,
                        ),
                        Some(topography) => danton::earth_model(
                            geoid.as_ptr(),
                            topography.as_ptr(), // Ensures that topography isn't moved.
                            self.density,
                            &mut ocean,
                        ),
                    }
                },
                None
            )?;
            self.modified = false;
        }
        Ok(())
    }

    fn geoid_undulation_from_any<'py>(
        &self,
        any: &Bound<'py, PyAny>,
    ) -> PyResult<PyObject> {
        let py = any.py();
        let projection = Projection::new(any)?;
        let size = projection.size();
        let undulations = PyArray::<f64>::empty(py, &projection.shape())?;

        for i in 0..size {
            let latitude = projection.latitude.get(i)?;
            let longitude = projection.longitude.get(i)?;
            let undulation = unsafe { danton::geoid_undulation(latitude, longitude) };
            undulations.set(i, undulation)?;
        }

        Ok(undulations.unbind(py))
    }

    fn to_ecef_from_any<'py>(
        &self,
        any: &Bound<'py, PyAny>
    ) -> PyResult<(PyObject, Option<PyObject>)> {
        let py = any.py();
        let position = Position::new(any)?;
        let direction = Direction::maybe_new(any)?;

        let (size, shape3) = match direction.as_ref() {
            None => (position.size(), position.shape3()),
            Some(direction) => {
                let (size, mut shape) = position.common(direction)?;
                shape.push(3);
                (size, shape)
            },
        };

        let ecef_position = PyArray::<f64>::empty(py, &shape3)?;
        let ecef_direction: Option<&PyArray<f64>> = match direction.as_ref() {
            None => None,
            Some(_) => Some(PyArray::<f64>::empty(py, &shape3)?),
        };

        for i in 0..size {
            let geodetic = position.get(i)?;
            let r = geodetic.to_ecef(self.geoid.into());
            ecef_position.set(3 * i, r[0])?;
            ecef_position.set(3 * i + 1, r[1])?;
            ecef_position.set(3 * i + 2, r[2])?;
            if let Some(direction) = direction.as_ref() {
                let horizontal = direction.get(i)?;
                let u = horizontal.to_ecef(self.geoid.into(), &geodetic);
                let ecef_direction = ecef_direction.unwrap();
                ecef_direction.set(3 * i, u[0])?;
                ecef_direction.set(3 * i + 1, u[1])?;
                ecef_direction.set(3 * i + 2, u[2])?;
            }
        }

        let ecef_position = ecef_position.unbind(py);

        let result = {
            match ecef_direction {
                None => (ecef_position, None),
                Some(ecef_direction) => {
                    let ecef_direction = ecef_direction.unbind(py);
                    (ecef_position, Some(ecef_direction))
                },
            }
        };

        Ok(result)
    }

    fn topography_elevation_from_any<'py>(
        &self,
        any: &Bound<'py, PyAny>,
        reference: Reference,
    ) -> PyResult<PyObject> {
        let py = any.py();
        let projection = Projection::new(any)?;
        let size = projection.size();
        let elevations = PyArray::<f64>::empty(py, &projection.shape())?;

        for i in 0..size {
            let latitude = projection.latitude.get(i)?;
            let longitude = projection.longitude.get(i)?;
            let mut elevation = unsafe { danton::topography_elevation(latitude, longitude) };
            if let Reference::Ellipsoid = reference {
                elevation += unsafe { danton::geoid_undulation(latitude, longitude) };
            }
            elevations.set(i, elevation)?;
        }

        Ok(elevations.unbind(py))
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

pub enum Mode {
    Full,
    Merge,
}

impl From<Mode> for c_uint {
    fn from(value: Mode) -> c_uint {
        match value {
            Mode::Full => danton::FULL,
            Mode::Merge => danton::MERGE,
        }
    }
}

impl<'a> Drop for Tracer<'a> {
    fn drop(&mut self) {
        unsafe {
            danton::destroy(&mut (self.tracer as *mut c_void));
        }
    }
}

impl<'a> Tracer<'a> {
    pub fn new(geometry: &'a mut Geometry, mode: Mode) -> PyResult<Self> {
        geometry.apply()?;
        let tracer = unsafe { danton::tracer_create(mode.into()) };
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
