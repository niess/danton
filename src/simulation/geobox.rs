use crate::simulation::random::Random;
use crate::utils::convert::{Array, Ellipsoid};
use crate::utils::coordinates::{Coordinates, CoordinatesExport, GeodeticCoordinates,
    GeodeticsExport, HorizontalCoordinates, HorizontalsExport};
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::ValueError;
use crate::utils::export::Export;
use crate::utils::extract::{Direction, Position};
use crate::utils::float::f64x3;
use crate::utils::numpy::PyArray;
use crate::utils::tuple::NamedTuple;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};


#[pyclass(name="Box", module="danton")]
pub struct GeoBox {
    /// Reference ellipsoid for geodetic coordinates.
    #[pyo3(get, set)]
    pub ellipsoid: Ellipsoid,
    /// Latitude coordinate of box centre, in deg.
    #[pyo3(get, set)]
    latitude: f64,
    /// Longitude coordinate of box centre, in deg.
    #[pyo3(get, set)]
    longitude: f64,
    /// Altitude coordinate of box centre, in m (w.r.t. the ellipsoid).
    #[pyo3(get, set)]
    altitude: f64,
    pub size: [f64; 3],
    /// Box inclination, in deg, w.r.t the North.
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

#[pymethods]
impl GeoBox {
    #[new]
    pub fn new(
        size: Size,
        latitude: Option<f64>,
        longitude: Option<f64>,
        altitude: Option<f64>,
        declination: Option<f64>,
        ellipsoid: Option<Ellipsoid>
    ) -> Self {
        let latitude = latitude.unwrap_or(0.0);
        let longitude = longitude.unwrap_or(0.0);
        let altitude = altitude.unwrap_or(0.0);
        let ellipsoid = ellipsoid.unwrap_or(Ellipsoid::default());
        let size: [f64; 3] = size.into();
        let declination = declination.unwrap_or(0.0);
        Self { ellipsoid, latitude, longitude, altitude, size, declination }
    }

    /// Box size along the x, y and z-axis, in m.
    #[getter]
    fn get_size<'py>(&self, py: Python<'py>) -> Bound<'py, PyTuple> {
        PyTuple::new_bound(py, self.size)
    }

    #[setter]
    fn set_size<'py>(&mut self, size: Size) {
        self.size = size.into();
    }

    /// Box surface area, in square-metres.
    #[getter]
    pub fn get_surface(&self) -> f64 {
        2.0 * (self.size[0] * (self.size[1] + self.size[2]) + self.size[1] * self.size[2])
    }

    /// Box cubic volume, in cubic-metres.
    #[getter]
    pub fn get_volume(&self) -> f64 {
        self.size[0] * self.size[1] * self.size[2]
    }

    /// Convert local cartesian coordinates to geodetic ones.
    fn from_local(
        &self,
        py: Python,
        position: Option<Array<f64>>,
        direction: Option<Array<f64>>
    ) -> PyResult<PyObject> {
        let position = position
            .map(|position| position.resolve());
        let direction = direction
            .map(|direction| direction.resolve());

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

        let ndim = {
            let mut ndim = 0;
            if let Some(position) = position {
                ndim = ndim.max(check_shape(position, "position")?)
            }
            if let Some(direction) = direction {
                ndim = ndim.max(check_shape(direction, "direction")?)
            }
            ndim
        };

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

    /// Check if elements are inside the box.
    #[pyo3(signature=(array=None, /, **kwargs))]
    fn inside<'py>(
        &self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let position = match array {
            Some(array) => Position::new(array)?,
            None => match kwargs {
                Some(kwargs) => Position::new(kwargs.as_any())?,
                None => {
                    let inside: &PyAny = PyArray::<bool>::empty(py, &[0])?;
                    return Ok(inside.into_py(py));
                },
            }
        };
        let frame = self.local_frame();
        let inside = PyArray::<bool>::empty(py, &position.shape())?;

        for i in 0..position.size() {
            let geodetic = position.get(i)?;
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

    /// Convert geodetic coordinates to local cartesian ones.
    #[pyo3(signature=(array=None, /, **kwargs))]
    fn to_local<'py>(
        &self,
        py: Python<'py>,
        array: Option<&Bound<'py, PyAny>>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<PyObject> {
        let (position, direction) = match array {
            Some(array) => (Position::maybe_new(array)?, Direction::maybe_new(array)?),
            None => match kwargs {
                Some(kwargs) => {
                    let any = kwargs.as_any();
                    (Position::maybe_new(any)?, Direction::maybe_new(any)?)
                },
                None => (None, None),
            }
        };

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
            let geodetic = match position.as_ref() {
                None => origin,
                Some(position) => {
                    let geodetic = position.get(i)?;
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

    /// Compute the projected surface area, in square-metres.
    fn projected_surface(&self, azimuth: f64, elevation: f64) -> f64 {
        let direction = HorizontalCoordinates { azimuth, elevation };
        let projection = ProjectedBox::new(&self, &direction);
        projection.surface()
    }
}

pub struct LocalFrame {
    ellipsoid: Ellipsoid,
    origin: f64x3,
    ux: f64x3,
    uy: f64x3,
    uz: f64x3,
}

impl GeoBox {
    pub fn local_frame(&self) -> LocalFrame {
        let geodetic = self.origin();
        let origin: f64x3 = (&geodetic.to_ecef(self.ellipsoid)).into();
        let ux = HorizontalCoordinates { azimuth: 90.0 + self.declination, elevation: 0.0 };
        let ux: f64x3 = (&ux.to_ecef(self.ellipsoid, &geodetic)).into();
        let uy = HorizontalCoordinates { azimuth: self.declination, elevation: 0.0 };
        let uy: f64x3 = (&uy.to_ecef(self.ellipsoid, &geodetic)).into();
        let uz = HorizontalCoordinates { azimuth: 0.0, elevation: 90.0 };
        let uz: f64x3 = (&uz.to_ecef(self.ellipsoid, &geodetic)).into();
        LocalFrame { ellipsoid: self.ellipsoid, origin, ux, uy, uz }
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
        let r_ecef = coordinates.to_ecef(self.ellipsoid);
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
        let u_ecef = coordinates.to_ecef(self.ellipsoid, origin);
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
        GeodeticCoordinates::from_ecef(&r_ecef, self.ellipsoid)
    }

    pub fn to_horizontal(
        &self,
        u: &[f64; 3],
        origin: &GeodeticCoordinates
    ) -> HorizontalCoordinates {
        let u_ecef = u[0] * self.ux + u[1] * self.uy + u[2] * self.uz;
        let u_ecef: [f64; 3] = u_ecef.into();
        HorizontalCoordinates::from_ecef(&u_ecef, self.ellipsoid, origin)
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


// ===============================================================================================
//
// Projected box.
//
// ===============================================================================================

pub struct ProjectedBox {
    frame: LocalFrame,
    u: f64x3,
    v: f64x3,
    w: f64x3,
    o: [f64; 2],
    k0: [f64; 2],
    k1: [f64; 2],
    k2: [f64; 2],
    cumulated_surface: [f64; 3],
}

impl ProjectedBox {
    pub fn new(geobox: &GeoBox, direction: &HorizontalCoordinates) -> Self {
        let frame = geobox.local_frame();
        let size = &geobox.size;
        let u: f64x3 = (&frame.from_horizontal(direction, &geobox.origin()))
            .into();
        let (v, w) = u.covectors();
        let e = [
            [v.x(), w.x()],
            [v.y(), w.y()],
            [v.z(), w.z()],
        ];

        #[derive(Copy, Clone)]
        enum Sign {
            Negative,
            Positive,
        }

        impl Sign {
            fn new(v: f64) -> Self {
                if v >= 0.0 { Sign::Positive } else { Sign::Negative }
            }

            fn flip(&mut self) {
                match self {
                    Sign::Negative => *self = Sign::Positive,
                    Sign::Positive => *self = Sign::Negative,
                }
            }
        }

        let node = |signs: &[Sign; 3]| -> [f64; 2] {
            let mut x = 0.0;
            let mut y = 0.0;
            for i in 0..3 {
                let di = match signs[i] {
                    Sign::Negative => -0.5 * size[i],
                    Sign::Positive => 0.5 * size[i],
                };
                x += di * e[i][0];
                y += di * e[i][1];
            }
            [x, y]
        };

        let signs = [
            Sign::new(u.x()),
            Sign::new(u.y()),
            Sign::new(u.z())
        ];
        let o = node(&signs);

        let vector = |i: usize| -> [f64; 2] {
            let mut signs = signs.clone();
            signs[i].flip();
            let k = node(&signs);
            [k[0] - o[0], k[1] - o[1]]
        };

        let k0 = vector(0);
        let k1 = vector(1);
        let k2 = vector(2);

        let area = |a: &[f64; 2], b: &[f64; 2]| -> f64 {
            (a[0] * b[1] - a[1] * b[0]).abs()
        };

        let mut cumulated_surface = [
            area(&k0, &k1),
            area(&k1, &k2),
            area(&k2, &k0),
        ];
        for i in 1..3 {
            cumulated_surface[i] += cumulated_surface[i - 1];
        }

        Self {
            frame, u, v, w, o, k0, k1, k2, cumulated_surface
        }
    }

    #[inline]
    pub fn surface(&self) -> f64 {
        self.cumulated_surface[2]
    }

    pub fn generate_inside(
        &self,
        random: &mut Random
    ) -> (GeodeticCoordinates, HorizontalCoordinates) {
        // select box side.
        let s = random.open01() * self.cumulated_surface[2];
        let mut side = 0;
        for (j, surface) in self.cumulated_surface.iter().enumerate() {
            if s <= *surface {
                side = j;
                break
            }
        }

        // Randomise local coordinates.
        let x = random.open01();
        let y = random.open01();
        let (kx, ky) = match side {
            0 => (&self.k0, &self.k1),
            1 => (&self.k1, &self.k2),
            2 => (&self.k2, &self.k0),
            _ => unreachable!(),
        };
        let (x, y) = (
            self.o[0] + x * kx[0] + y * ky[0],
            self.o[1] + x * kx[1] + y * ky[1],
        );

        let r: [f64; 3] = (x * self.v + y * self.w).into();
        let position = self.frame.to_geodetic(&r);
        let direction: [f64; 3] = self.u.into();
        let direction = self.frame.to_horizontal(&direction, &position);
        (position, direction)
    }
}
