use crate::simulation::random::Random;
use crate::utils::convert::Geodesic;
use crate::utils::coordinates::{Coordinates, GeodeticCoordinates, HorizontalCoordinates};
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{TypeError, ValueError};
use crate::utils::export::Export;
use crate::utils::float::f64x3;
use crate::utils::numpy::PyArray;
use crate::utils::tuple::NamedTuple;
use derive_more::{AsMut, AsRef, From};
use pyo3::prelude::*;
use pyo3::types::PyTuple;


#[pyclass(name="Box", module="danton")]
pub struct GeoBox {
    #[pyo3(get, set)]
    pub geodesic: Geodesic,
    #[pyo3(get, set)]
    latitude: f64,
    #[pyo3(get, set)]
    longitude: f64,
    #[pyo3(get, set)]
    altitude: f64,
    pub size: [f64; 3],
    #[pyo3(get, set)]
    declination: f64,
}

#[derive(FromPyObject)]
pub enum Size {
    #[pyo3(transparent, annotation = "[f64;3]")]
    Vector3([f64; 3]),
    #[pyo3(transparent, annotation = "[f64;2]")]
    Vector2([f64; 2]),
    #[pyo3(transparent, annotation = "f64")]
    Scalar(f64),
}

impl From<Size> for [f64; 3] {
    fn from(value: Size) -> Self {
        match value {
            Size::Vector3(size) => size,
            Size::Vector2(size) => [size[0], size[0], size[1]],
            Size::Scalar(size) => [size, size, size],
        }
    }
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
pub struct CoordinatesExport (Export<Coordinates>);

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
pub struct GeodeticsExport (Export<GeodeticCoordinates>);

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
pub struct HorizontalsExport (Export<HorizontalCoordinates>);

#[pymethods]
impl GeoBox {
    #[new]
    pub fn new(
        size: Size,
        latitude: Option<f64>,
        longitude: Option<f64>,
        altitude: Option<f64>,
        declination: Option<f64>,
        geodesic: Option<Geodesic>
    ) -> Self {
        let latitude = latitude.unwrap_or(0.0);
        let longitude = longitude.unwrap_or(0.0);
        let altitude = altitude.unwrap_or(0.0);
        let geodesic = geodesic.unwrap_or(Geodesic::default());
        let size: [f64; 3] = size.into();
        let declination = declination.unwrap_or(0.0);
        Self { geodesic, latitude, longitude, altitude, size, declination }
    }

    #[getter]
    pub fn get_area(&self) -> f64 {
        2.0 * (self.size[0] * (self.size[1] + self.size[2]) + self.size[1] * self.size[2])
    }

    #[getter]
    fn get_size<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        PyTuple::new_bound(py, self.size)
    }

    #[setter]
    fn set_size<'py>(&mut self, size: Size) {
        self.size = size.into();
    }

    fn from_local(
        &self,
        py: Python,
        position: Option<&PyArray<f64>>,
        direction: Option<&PyArray<f64>>
    ) -> PyResult<PyObject> {
        enum Case {
            Both,
            Direction,
            Position,
        }

        let case = if position.is_some() {
            if direction.is_some() {
                Case::Both
            } else {
                Case::Position
            }
        } else {
            if direction.is_some() {
                Case::Direction
            } else {
                return Ok(py.None())
            }
        };

        let size = match case {
            Case::Both => {
                let size = position.unwrap().size();
                if direction.unwrap().size() != size {
                    let why = format!(
                        "expected a size {} array, found a size {} array",
                        size,
                        direction.unwrap().size(),
                    );
                    let err = Error::new(ValueError)
                        .what("direction")
                        .why(&why);
                    return Err(err.to_err())
                }
                size
            },
            Case::Direction => direction.unwrap().size(),
            Case::Position => position.unwrap().size(),
        };

        fn check_shape(array: &PyArray<f64>, what: &str) -> PyResult<()> {
            let shape = array.shape();
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

        if let Some(position) = position {
            check_shape(position, "position")?
        }
        if let Some(direction) = direction {
            check_shape(direction, "direction")?
        }

        let (mut coordinates, mut horizontals, mut geodetics) = match case {
            Case::Both => {
                let coordinates: Vec<Coordinates> = Vec::new();
                (Some(coordinates), None, None)
            },
            Case::Direction => {
                let horizontals: Vec<HorizontalCoordinates> = Vec::new();
                (None, Some(horizontals), None)
            },
            Case::Position => {
                let geodetics: Vec<GeodeticCoordinates> = Vec::new();
                (None, None, Some(geodetics))
            },
        };

        let frame = self.local_frame();
        let origin = self.origin();

        for i in 0..(size / 3) {
            let geodetic = match case {
                Case::Both | Case::Position => {
                    let position = position.unwrap();
                    let position = [
                        position.get(3 * i)?,
                        position.get(3 * i + 1)?,
                        position.get(3 * i + 2)?,
                    ];
                    frame.to_geodetic(&position)
                },
                Case::Direction => origin,
            };
            let horizontal = match case {
                Case::Both | Case::Direction => {
                    let direction = direction.unwrap();
                    let direction = [
                        direction.get(3 * i)?,
                        direction.get(3 * i + 1)?,
                        direction.get(3 * i + 2)?,
                    ];
                    let horizontal = frame.to_horizontal(&direction, &geodetic);
                    Some(horizontal)
                },
                _ => None,
            };
            match case {
                Case::Both => {
                    let coordinate = Coordinates {
                        position: geodetic,
                        direction: horizontal.unwrap(),
                    };
                    coordinates.as_mut().map(|v| v.push(coordinate));
                },
                Case::Direction => {
                    let horizontal = horizontal.unwrap();
                    horizontals.as_mut().map(|v| v.push(horizontal));
                },
                Case::Position => {
                    geodetics.as_mut().map(|v| v.push(geodetic));
                },
            }
        }

        let result = match case {
            Case::Both => Export::export::<CoordinatesExport>(py, coordinates.unwrap())?,
            Case::Direction => Export::export::<HorizontalsExport>(py, horizontals.unwrap())?,
            Case::Position => Export::export::<GeodeticsExport>(py, geodetics.unwrap())?,
        };
        Ok(result)
    }

    fn inside<'py>(&self, elements: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let py = elements.py();
        if elements.is_none() {
            let inside: &PyAny = PyArray::<bool>::empty(py, &[0])?;
            return Ok(inside.into_py(py));
        }
        let position = Position::new(elements)?
            .ok_or_else(|| {
                let why = "missing 'latitude', 'longitude' and 'altitude'";
                let err = Error::new(ValueError)
                    .what("position")
                    .why(why);
                err.to_err()
            })?;

        let frame = self.local_frame();
        let inside = PyArray::<bool>::empty(py, &position.shape())?;

        for i in 0..position.size() {
            let geodetic = GeodeticCoordinates {
                latitude: position.latitude.get(i)?,
                longitude: position.longitude.get(i)?,
                altitude: position.altitude.get(i)?,
            };
            let r = frame.from_geodetic(&geodetic);
            const EPS: f64 = ::std::f32::EPSILON as f64;
            let b =
                (r[0].abs() <= 0.5 * self.size[0] + EPS) &&
                (r[1].abs() <= 0.5 * self.size[1] + EPS) &&
                (r[2].abs() <= 0.5 * self.size[2] + EPS);
            inside.set(i, b)?;
        }

        let inside: &PyAny = inside;
        Ok(inside.into_py(py))
    }

    fn to_local<'py>(&self, elements: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        let py = elements.py();
        let position = Position::new(elements)?;
        let direction = Direction::new(elements)?;

        let size = match position.as_ref() {
            None => match direction.as_ref() {
                None => return Ok(py.None()),
                Some(direction) => direction.size(),
            },
            Some(position) => match direction.as_ref() {
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
            },
        };

        let local_position: Option<&PyArray<f64>> = match position.as_ref() {
            None => None,
            Some(position) => Some(PyArray::<f64>::empty(py, &position.shape3())?),
        };
        let local_direction: Option<&PyArray<f64>> = match direction.as_ref() {
            None => None,
            Some(direction) => Some(PyArray::<f64>::empty(py, &direction.shape3())?),
        };

        let frame = self.local_frame();
        let origin = self.origin();

        for i in 0..size {
            let geodetic = match position {
                None => origin,
                Some(Position { latitude, longitude, altitude }) => {
                    let geodetic = GeodeticCoordinates {
                        latitude: latitude.get(i)?,
                        longitude: longitude.get(i)?,
                        altitude: altitude.get(i)?,
                    };
                    let r = frame.from_geodetic(&geodetic);
                    let local_position = local_position.unwrap();
                    local_position.set(3 * i, r[0])?;
                    local_position.set(3 * i + 1, r[1])?;
                    local_position.set(3 * i + 2, r[2])?;
                    geodetic
                },
            };
            if let Some(Direction { azimuth, elevation }) = direction {
                let horizontal = HorizontalCoordinates {
                    azimuth: azimuth.get(i)?,
                    elevation: elevation.get(i)?,
                };
                let u = frame.from_horizontal(&horizontal, &geodetic);
                let local_direction = local_direction.unwrap();
                local_direction.set(3 * i, u[0])?;
                local_direction.set(3 * i + 1, u[1])?;
                local_direction.set(3 * i + 2, u[2])?;
            }
        }

        let result = match local_position {
            None => {
                let direction: &PyAny = local_direction.unwrap();
                direction.into_py(py)
            },
            Some(position) => {
                let position: &PyAny = position;
                match local_direction {
                    None => position.into_py(py),
                    Some(direction) => {
                        static RESULT: NamedTuple<2> = NamedTuple::new(
                            "Result", ["position", "direction"]);
                        let direction: &PyAny = direction;
                        RESULT.instance(py, (position, direction))?.unbind()
                    }
                }
            },
        };
        Ok(result)
    }
}

struct Direction<'a> {
    azimuth: &'a PyArray<f64>,
    elevation: &'a PyArray<f64>,
}

impl<'a> Direction<'a> {
    fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
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

    fn shape(&self) -> Vec<usize> {
        self.azimuth.shape()
    }

    fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    fn size(&self) -> usize {
        self.azimuth.size()
    }
}

struct Position<'a> {
    latitude: &'a PyArray<f64>,
    longitude: &'a PyArray<f64>,
    altitude: &'a PyArray<f64>,
}

impl<'a> Position<'a> {
    fn new<'py: 'a>(elements: &'a Bound<'py, PyAny>) -> PyResult<Option<Self>> {
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

    fn shape(&self) -> Vec<usize> {
        self.latitude.shape()
    }

    fn shape3(&self) -> Vec<usize> {
        let mut shape = self.shape();
        shape.push(3);
        shape
    }

    fn size(&self) -> usize {
        self.latitude.size()
    }
}

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

pub struct LocalFrame {
    geodesic: Geodesic,
    origin: f64x3,
    ux: f64x3,
    uy: f64x3,
    uz: f64x3,
}

impl GeoBox {
    pub fn local_frame(&self) -> LocalFrame {
        let geodetic = self.origin();
        let origin: f64x3 = (&geodetic.to_ecef(self.geodesic)).into();
        let ux = HorizontalCoordinates { azimuth: 90.0 + self.declination, elevation: 0.0 };
        let ux: f64x3 = (&ux.to_ecef(self.geodesic, &geodetic)).into();
        let uy = HorizontalCoordinates { azimuth: self.declination, elevation: 0.0 };
        let uy: f64x3 = (&uy.to_ecef(self.geodesic, &geodetic)).into();
        let uz = HorizontalCoordinates { azimuth: 0.0, elevation: 90.0 };
        let uz: f64x3 = (&uz.to_ecef(self.geodesic, &geodetic)).into();
        LocalFrame { geodesic: self.geodesic, origin, ux, uy, uz }
    }

    pub fn origin(&self) -> GeodeticCoordinates {
        GeodeticCoordinates {
            latitude: self.latitude,
            longitude: self.longitude,
            altitude: self.altitude,
        }
    }
}

impl LocalFrame {
    pub fn from_geodetic(&self, coordinates: &GeodeticCoordinates) -> [f64; 3] {
        let r_ecef = coordinates.to_ecef(self.geodesic);
        let r_ecef: f64x3 = (&r_ecef).into();
        let r = r_ecef - self.origin;
        [
            r.dot(self.ux),
            r.dot(self.uy),
            r.dot(self.uz),
        ]
    }

    pub fn from_horizontal(
        &self,
        coordinates: &HorizontalCoordinates,
        origin: &GeodeticCoordinates
    ) -> [f64; 3] {
        let u_ecef = coordinates.to_ecef(self.geodesic, origin);
        let u_ecef: f64x3 = (&u_ecef).into();
        [
            u_ecef.dot(self.ux),
            u_ecef.dot(self.uy),
            u_ecef.dot(self.uz),
        ]
    }

    pub fn to_geodetic(&self, r: &[f64; 3]) -> GeodeticCoordinates {
        let r_ecef = r[0] * self.ux + r[1] * self.uy + r[2] * self.uz + self.origin;
        let r_ecef: [f64; 3] = r_ecef.into();
        GeodeticCoordinates::from_ecef(&r_ecef, self.geodesic)
    }

    pub fn to_horizontal(
        &self,
        u: &[f64; 3],
        origin: &GeodeticCoordinates
    ) -> HorizontalCoordinates {
        let u_ecef = u[0] * self.ux + u[1] * self.uy + u[2] * self.uz;
        let u_ecef: [f64; 3] = u_ecef.into();
        HorizontalCoordinates::from_ecef(&u_ecef, self.geodesic, origin)
    }
}


// ===============================================================================================
//
// Generator interface.
//
// ===============================================================================================

pub struct BoxGenerator {
    frame: LocalFrame,
    size: [f64; 3],
    cumulated_surface: [f64; 6],
    volume: f64,
}

impl BoxGenerator {
    pub fn new(geobox: &GeoBox) -> Self {
        let frame = geobox.local_frame();
        let size = geobox.size;
        let cumulated_surface = {
            let mut sides = [
                size[1] * size[2],
                size[1] * size[2],
                size[0] * size[2],
                size[0] * size[2],
                size[0] * size[1],
                size[0] * size[1],
            ];
            for i in 1..5 {
                sides[i] += sides[i - 1];
            }
            sides
        };
        let volume = size[0] * size[1] * size[2];
        Self { frame, size, cumulated_surface, volume }
    }

    #[inline]
    pub fn frame(&self) -> &LocalFrame {
        &self.frame
    }

    pub fn generate_inside(&self, random: &mut Random) -> GeodeticCoordinates {
        let r = [
            self.size[0] * (random.open01() - 0.5),
            self.size[1] * (random.open01() - 0.5),
            self.size[2] * (random.open01() - 0.5),
        ];
        self.frame.to_geodetic(&r)
    }

    pub fn generate_onto(
        &self,
        random: &mut Random
    ) -> ([f64; 3], [f64;3], GeodeticCoordinates, HorizontalCoordinates) {
        let (mut r, mut u) = {
            // select box side.
            let s = random.open01() * self.cumulated_surface[5];
            let mut side = 0;
            for (j, surface) in self.cumulated_surface.iter().enumerate() {
                if s <= *surface {
                    side = j;
                    break
                }
            }

            // Randomise local coordinates.
            let u = random.open01() - 0.5;
            let v = random.open01() - 0.5;
            match side {
                0 => ([  0.5, u, v ], f64x3::new(-1.0, 0.0, 0.0)),
                1 => ([ -0.5, u, v ], f64x3::new( 1.0, 0.0, 0.0)),
                2 => ([ u,  0.5, v ], f64x3::new(0.0, -1.0, 0.0)),
                3 => ([ u, -0.5, v ], f64x3::new(0.0,  1.0, 0.0)),
                4 => ([ u, v,  0.5 ], f64x3::new(0.0, 0.0, -1.0)),
                5 => ([ u, v, -0.5 ], f64x3::new(0.0, 0.0,  1.0)),
                _ => unreachable!(),
            }
        };
        for j in 0..3 {
            r[j] *= self.size[j];
        }
        let cos_theta = random.open01().sqrt();
        let phi = 2.0 * ::std::f64::consts::PI * random.open01();
        u.rotate(cos_theta, phi);
        let u: [f64; 3] = u.into();

        // Transform to geodetic.
        let geodetic = self.frame.to_geodetic(&r);
        let horizontal = self.frame.to_horizontal(&u, &geodetic);
        (r, u, geodetic, horizontal)
    }

    #[inline]
    pub fn size(&self) -> &[f64; 3] {
        &self.size
    }

    #[inline]
    pub fn surface(&self) -> f64 {
        self.cumulated_surface[5]
    }

    #[inline]
    pub fn volume(&self) -> f64 {
        self.volume
    }
}
