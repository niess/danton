use crate::bindings::danton;
use crate::simulation::geometry::Geometry;
use crate::simulation::random::Random;
use crate::utils::convert::Mode;
use crate::utils::coordinates::{GeodeticCoordinates, HorizontalCoordinates};
use crate::utils::error::{ctrlc_catched, Error};
use crate::utils::error::ErrorKind::{KeyboardInterrupt, KeyError, NotImplementedError, ValueError};
use crate::utils::float::f64x3;
use crate::utils::numpy::{Dtype, PyArray, ShapeArg};
use crate::utils::tuple::NamedTuple;
use pyo3::prelude::*;
use pyo3::types::PyDict;
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
    pub weight: f64,
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
#[pyo3(signature=(shape, **kwargs))]
pub fn particles(
    py: Python,
    shape: ShapeArg,
    kwargs: Option<&Bound<PyDict>>
) -> PyResult<PyObject> {
    let shape: Vec<usize> = shape.into();
    let array: &PyAny = PyArray::<Particle>::zeros(py, &shape)?;
    let mut has_pid = false;
    let mut has_energy = false;
    let mut has_weight = false;
    if let Some(kwargs) = kwargs {
        for (key, value) in kwargs.iter() {
            {
                let key: String = key.extract()?;
                match key.as_str() {
                    "pid" => { has_pid = true; },
                    "energy" => { has_energy = true; },
                    "weight" => { has_weight = true; },
                    _ => {},
                }
            }
            array.set_item(key, value)?;
        }
    }
    if !has_pid {
        array.set_item("pid", 15)?;
    }
    if !has_energy {
        array.set_item("energy", 1E+09)?;
    }
    if !has_weight {
        array.set_item("weight", 1.0)?;
    }
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
    weight: Option<&'a PyArray<f64>>,
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
        let weight = extract(elements, "weight").ok();
        let size = pid.size();
        let others = [
            energy.size(), latitude.size(), longitude.size(), altitude.size(), azimuth.size(),
            elevation.size(), weight.map(|a| a.size()).unwrap_or(size),
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
            pid, energy, latitude, longitude, altitude, azimuth, elevation, weight, size, index
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
        let weight = None;
        let iter = Self {
            pid, energy, latitude, longitude, altitude, azimuth, elevation, weight, size, index
        };
        Ok(iter)
    }

    fn get(&self, index: usize) -> PyResult<Particle> {
        let (pid, energy) = match self.pid {
            None => (0, 0.0),
            Some(pid) => (pid.get(index)?, self.energy.unwrap().get(index)?),
        };
        let weight = match self.weight {
            None => 1.0,
            Some(weight) => weight.get(index)?,
        };
        let particle = Particle {
            pid,
            energy,
            latitude: self.latitude.get(index)?,
            longitude: self.longitude.get(index)?,
            altitude: self.altitude.get(index)?,
            azimuth: self.azimuth.get(index)?,
            elevation: self.elevation.get(index)?,
            weight,
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


// ===============================================================================================
//
// Particles generator.
//
// ===============================================================================================

#[pyclass(module="danton")]

pub struct ParticlesGenerator {
    particles: PyObject,
    random: Py<Random>,
    geometry: Py<Geometry>,
    weight: bool,
    mode_hint: Option<Mode>,
    // Status flags.
    is_pid: bool,
    is_energy: bool,
    is_position: bool,
    is_direction: bool,
    // Sample size.
    size: Option<usize>,
}

#[pymethods]
impl ParticlesGenerator {
    const PI: f64 = ::std::f64::consts::PI;

    #[new]
    pub fn new<'py>(
        py: Python<'py>,
        shape: ShapeArg,
        geometry: Option<Bound<'py, Geometry>>,
        random: Option<Bound<'py, Random>>,
        weight: Option<bool>,
    ) -> PyResult<Self> {
        let shape: Vec<usize> = shape.into();
        let particles = {
            let particles = PyArray::<Particle>::zeros(py, &shape)?;
            for pi in unsafe { particles.slice_mut()? } {
                pi.weight = 1.0;
            }
            let particles: &PyAny = particles;
            let particles: PyObject = particles.into();
            particles
        };
        let random = match random {
            None => Py::new(py, Random::new(None)?)?,
            Some(random) => random.unbind(),
        };
        let geometry = match geometry {
            None => Py::new(py, Geometry::new())?,
            Some(geometry) => geometry.unbind(),
        };
        let weight = weight.unwrap_or(true);
        let generator = Self {
            particles, random, geometry, weight,
            mode_hint: None,
            is_pid: false,
            is_energy: false,
            is_position: false,
            is_direction: false,
            size: None,
        };
        Ok(generator)
    }

    fn direction<'py>(
        slf: Bound<'py, Self>,
        azimuth: &Bound<'py, PyAny>,
        elevation: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let mut generator = slf.borrow_mut();
        if generator.is_direction {
            let err = Error::new(ValueError)
                .what("direction")
                .why("already defined");
            return Err(err.to_err())
        }
        generator.particles
            .bind(py)
            .set_item("azimuth", azimuth)?;
        generator.particles
            .bind(py)
            .set_item("elevation", elevation)?;
        generator.is_direction = true;
        Ok(slf)
    }

    fn energy<'py>(
        slf: Bound<'py, Self>,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let mut generator = slf.borrow_mut();
        if generator.is_energy {
            let err = Error::new(ValueError)
                .what("energy")
                .why("already defined");
            return Err(err.to_err())
        }
        generator.particles
            .bind(py)
            .set_item("energy", value)?;
        generator.is_energy = true;
        Ok(slf)
    }

    fn generate<'py>(&mut self, py: Python<'py>) -> PyResult<Bound<'py, PyAny>> {
        if !self.is_pid {
            let mode = self.mode_hint.unwrap_or(Mode::Backward);
            let pid: i32 = match mode {
                Mode::Backward => 15,
                _ => 16,
            };
            self.particles.bind(py).set_item("pid", pid)?;
            self.is_pid = true;
        }

        if !self.is_energy {
            self.particles.bind(py).set_item("energy", 1E+09)?; // XXX Power law as default.
            self.is_energy = true;
        }

        if !self.is_position {
            // Positions are already initialised to 0.
            self.is_position = true;
        }

        if !self.is_direction {
            self.generate_solid_angle(py, None, Some([-10.0, 10.0]), None)?;
        }

        let particles = self.particles.bind(py).clone().into_any();
        let result = match self.size {
            None => particles,
            Some(size) => {
                static RESULT: NamedTuple<2> = NamedTuple::new(
                    "Result", ["particles", "size"]);
                RESULT.instance(py, (particles, size))?
            },
        };
        Ok(result)
    }

    fn pid<'py>(
        slf: Bound<'py, Self>,
        value: &Bound<'py, PyAny>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let mut generator = slf.borrow_mut();
        if generator.is_pid {
            let err = Error::new(ValueError)
                .what("pid")
                .why("already defined");
            return Err(err.to_err())
        }
        generator.particles
            .bind(py)
            .set_item("pid", value)?;
        generator.is_pid = true;
        Ok(slf)
    }

    fn position<'py>(
        slf: Bound<'py, Self>,
        latitude: &Bound<'py, PyAny>,
        longitude: &Bound<'py, PyAny>,
        altitude: Option<&Bound<'py, PyAny>>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let mut generator = slf.borrow_mut();
        if generator.is_position {
            let err = Error::new(ValueError)
                .what("position")
                .why("already defined");
            return Err(err.to_err())
        }
        generator.particles
            .bind(py)
            .set_item("latitude", latitude)?;
        generator.particles
            .bind(py)
            .set_item("longitude", longitude)?;
        if let Some(altitude) = altitude {
            generator.particles
                .bind(py)
                .set_item("altitude", altitude)?;
        }
        generator.is_position = true;
        Ok(slf)
    }

    fn powerlaw<'py>( // XXX Document all these methods (+docstring)
        slf: Bound<'py, Self>,
        energy_min: f64,
        energy_max: f64,
        exponent: Option<f64>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let exponent = exponent.unwrap_or(-1.0);
        if energy_min >= energy_max || energy_min <= 0.0 {
            let why = "expected energy_max > energy_min > 0.0";
            let err = Error::new(ValueError).what("powerlaw").why(why);
            return Err(err.to_err());
        }

        let py = slf.py();
        let generator = slf.borrow();
        let weight = weight.unwrap_or(generator.weight);

        let particles = generator.particles
            .bind(py);
        let particles_energy: &PyArray<f64> = particles
            .get_item("energy")?
            .extract()?;
        let weights: Option<&PyArray<f64>> = if weight {
            let weights: &PyArray<f64> = particles
                .get_item("weight")?
                .extract()?;
            Some(weights)
        } else {
            None
        };

        let mut random = generator.random.bind(py).borrow_mut();
        for i in 0..particles_energy.size() {
            let (energy, energy_weight) = match exponent {
                -1.0 => {
                    let lne = (energy_max / energy_min).ln();
                    let energy = energy_min * (random.open01() * lne).exp();
                    (energy, energy * lne)
                },
                0.0 => {
                    let de = energy_max - energy_min;
                    let energy = de * random.open01() + energy_min;
                    (energy, de)
                },
                exponent => {
                    let a = exponent + 1.0;
                    let b = energy_min.powf(a);
                    let de = energy_max.powf(a) - b;
                    let energy = (de * random.open01() + b).powf(1.0 / a);
                    let weight = de / (a * energy.powf(exponent));
                    (energy, weight) // XXX Test all 3 cases.
                },
            };
            particles_energy.set(i, energy)?; // XXX clip.
            if let Some(weights) = weights {
                let weight = weights.get(i)? * energy_weight;
                weights.set(i, weight)?;
            }
        }

        drop(generator);
        let mut generator = slf.borrow_mut();
        generator.is_energy = true;

        Ok(slf)
    }

    fn solid_angle<'py>(
        slf: Bound<'py, Self>,
        azimuth: Option<[f64; 2]>,
        elevation: Option<[f64; 2]>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let mut generator = slf.borrow_mut();
        generator.generate_solid_angle(py, azimuth, elevation, weight)?;
        Ok(slf)
    }

    /// Target a box volume.
    fn target<'py>(
        slf: Bound<'py, Self>,
        latitude: f64,
        longitude: f64,
        altitude: f64,
        size: [f64; 3],
        elevation: Option<[f64; 2]>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let py = slf.py();
        let generator = slf.borrow();
        if generator.is_position {
            let err = Error::new(ValueError)
                .what("target")
                .why("position already defined");
            return Err(err.to_err())
        }
        if generator.is_direction {
            let err = Error::new(ValueError)
                .what("target")
                .why("direction already defined");
            return Err(err.to_err())
        }

        // Set local frame.
        let geodesic = generator.geometry.bind(py).borrow().geodesic;
        let geodetic = GeodeticCoordinates { latitude, longitude, altitude };
        let origin: f64x3 = (&geodetic.to_ecef(geodesic)).into();
        let ux = HorizontalCoordinates { azimuth: 90.0, elevation: 0.0 };
        let ux: f64x3 = (&ux.to_ecef(geodesic, &geodetic)).into();
        let uy = HorizontalCoordinates { azimuth: 0.0, elevation: 0.0 };
        let uy: f64x3 = (&uy.to_ecef(geodesic, &geodetic)).into();
        let uz = HorizontalCoordinates { azimuth: 0.0, elevation: 90.0 };
        let uz: f64x3 = (&uz.to_ecef(geodesic, &geodetic)).into();

        // Sides surfaces.
        let sides = [
            size[1] * size[2],
            size[1] * size[2],
            size[0] * size[2],
            size[0] * size[2],
            size[0] * size[1],
            size[0] * size[1],
        ];
        let surface: f64 = sides.iter().sum();

        // Bind particles properties.
        let particles = generator.particles.bind(py);
        let latitudes: &PyArray<f64> = particles
            .get_item("latitude")?
            .extract()?;
        let longitudes: &PyArray<f64> = particles
            .get_item("longitude")?
            .extract()?;
        let altitudes: &PyArray<f64> = particles
            .get_item("altitude")?
            .extract()?;
        let azimuths: &PyArray<f64> = particles
            .get_item("azimuth")?
            .extract()?;
        let elevations: &PyArray<f64> = particles
            .get_item("elevation")?
            .extract()?;
        let weights = if weight.unwrap_or(generator.weight) {
            let weight: &PyArray<f64> = particles
                .get_item("weight")?
                .extract()?;
            Some(weight)
        } else {
            None
        };

        // Bind the PRNG.
        let mut random = generator.random.bind(py).borrow_mut();

        // Loop over particles.
        let mut trials: usize = 0;
        let n = latitudes.size();
        for i in 0..n {
            loop {
                let (mut r, mut u) = {
                    // select box side.
                    let s = random.open01() * surface;
                    let mut side = 0;
                    let mut sum = 0.0;
                    for (j, surface) in sides.iter().enumerate() {
                        sum += *surface;
                        if s <= sum {
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
                    r[j] *= size[j];
                }
                let cos_theta = random.open01().sqrt();
                let phi = 2.0 * Self::PI * random.open01();
                u.rotate(cos_theta, phi);

                // Transform to ECEF.
                let r_ecef = r[0] * ux + r[1] * uy + r[2] * uz + origin;
                let r_ecef: [f64; 3] = r_ecef.into();
                let u_ecef = u.x() * ux + u.y() * uy + u.z() * uz;
                let u_ecef: [f64; 3] = u_ecef.into();

                // Transform to geodetic.
                let geodetic = GeodeticCoordinates::from_ecef(&r_ecef, geodesic);
                let horizontal = HorizontalCoordinates::from_ecef(&u_ecef, geodesic, &geodetic);

                trials += 1;
                if ((trials % 1000) == 0) && ctrlc_catched() {
                    return Err(Error::new(KeyboardInterrupt).to_err())
                }
                if let Some([el_min, el_max]) = elevation {
                    // Find the second (outgoing) intersection (using the local frame).
                    // Ref: https://iquilezles.org/articles/intersectors/ (axis-aligned box).
                    let m = [ 1.0 / u.x(), 1.0 / u.y(), 1.0 / u.z() ];
                    let n = [ m[0] * r[0], m[1] * r[1], m[2] * r[2] ];
                    let k = [
                        0.5 * m[0].abs() * size[0],
                        0.5 * m[1].abs() * size[1],
                        0.5 * m[2].abs() * size[2],
                    ];
                    let t2 = [ -n[0] + k[0], -n[1] + k[1], -n[2] + k[2] ];
                    let tf = t2[0].min(t2[1]).min(t2[2]);
                    let r2 = [
                        r[0] + tf * u.x(),
                        r[1] + tf * u.y(),
                        r[2] + tf * u.z(),
                    ];

                    // Transform to ECEF, and then to geodetic coordinates.
                    let r2_ecef = r2[0] * ux + r2[1] * uy + r2[2] * uz + origin;
                    let r2_ecef: [f64; 3] = r2_ecef.into();
                    let geodetic2 = GeodeticCoordinates::from_ecef(&r2_ecef, geodesic);
                    let horizontal2 = HorizontalCoordinates::from_ecef(
                        &u_ecef, geodesic, &geodetic2
                    );

                    // Check the elevation values (inclusively).
                    let el0 = horizontal.elevation.min(horizontal2.elevation);
                    let el1 = horizontal.elevation.max(horizontal2.elevation);
                    if (el1 < el_min) || (el0 > el_max) {
                        continue;
                    }
                }

                // Update particle properties.
                latitudes.set(i, geodetic.latitude)?;
                longitudes.set(i, geodetic.longitude)?;
                altitudes.set(i, geodetic.altitude)?;
                azimuths.set(i, horizontal.azimuth)?;
                elevations.set(i, horizontal.elevation)?;
                if let Some(weights) = weights {
                    let w = weights.get(i)? * surface * Self::PI;
                    weights.set(i, w)?;
                }
                break;
            }
        }

        drop(generator);
        let mut generator = slf.borrow_mut();
        generator.is_position = true;
        generator.is_direction = true;
        if generator.mode_hint.is_none() {
            generator.mode_hint = Some(Mode::Forward);
        }
        if elevation.is_some() {
            if let Some(s) = generator.size {
                trials += s - n;
            }
            generator.size = Some(trials);
        }

        Ok(slf)
    }
}

impl ParticlesGenerator {
    const DEG: f64 = 180.0 / std::f64::consts::PI;
    const RAD: f64 = std::f64::consts::PI / 180.0;

    fn generate_solid_angle(
        &mut self,
        py: Python,
        azimuth: Option<[f64; 2]>,
        elevation: Option<[f64; 2]>,
        weight: Option<bool>,
    ) -> PyResult<()> {
        let weight = weight.unwrap_or(self.weight);
        let (sin_el0, sin_el1) = match elevation {
            None => (-1.0, 1.0),
            Some([el0, el1]) => {
                let el0 = el0 * Self::RAD; // XXX Bound to [-90.0, 90.0].
                let el1 = el1 * Self::RAD;
                (el0.sin(), el1.sin())
            },
        };

        let (az0, az1) = match azimuth {
            None => (-180.0, 180.0),
            Some([az0, az1]) => {
                let az0 = az0 % 360.0; // XXX Bound to [-180, 180].
                let az1 = az1 % 360.0;
                (az0, az1)
            },
        };

        let particles = self.particles.bind(py);
        let azimuths: &PyArray<f64> = particles
            .get_item("azimuth")?
            .extract()?;
        let elevations: &PyArray<f64> = particles
            .get_item("elevation")?
            .extract()?;
        let weights: Option<&PyArray<f64>> = if weight {
            let weights: &PyArray<f64> = particles
                .get_item("weight")?
                .extract()?;
            Some(weights)
        } else {
            None
        };
        let mut random = self.random.bind(py).borrow_mut();
        let n = azimuths.size();
        for i in 0..n {
            let sin_el = (sin_el1 - sin_el0) * random.open01() + sin_el0;
            let el = sin_el.asin() * Self::DEG;
            let az = (az1 - az0) * random.open01() + az0;
            azimuths.set(i, az)?;
            elevations.set(i, el)?;
        }
        if let Some(weights) = weights {
            let solid_angle = (sin_el1 - sin_el0).abs() * (az1 - az0).abs() * Self::RAD;
            for i in 0..n {
                let weight = weights.get(i)? * solid_angle;
                weights.set(i, weight)?;
            }
        }

        self.is_direction = true;
        Ok(())
    }
}
