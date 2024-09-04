use crate::bindings::danton;
use crate::simulation::geobox::{GeoBox, BoxGenerator, ProjectedBox};
use crate::simulation::geometry::{Geometry, Mode, Tracer};
use crate::simulation::random::Random;
use crate::utils::coordinates::HorizontalCoordinates;
use crate::utils::error::{ctrlc_catched, Error};
use crate::utils::error::ErrorKind::{KeyboardInterrupt, KeyError, NotImplementedError, ValueError};
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
#[pyo3(signature=(shape=None, /, **kwargs))]
pub fn particles(
    py: Python,
    shape: Option<ShapeArg>,
    kwargs: Option<&Bound<PyDict>>
) -> PyResult<PyObject> {
    let shape: Vec<usize> = match shape {
        Some(shape) => shape.into(),
        None => Vec::new(),
    };
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

pub fn extract<'a, 'py, T>(elements: &'a Bound<'py, PyAny>, key: &str) -> PyResult<&'a PyArray<T>>
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
    geometry: Py<Geometry>,
    random: Py<Random>,
    // Configuration.
    direction: Direction,
    energy: Energy,
    pid: Option<c_int>,
    position: Position,
    // Weight flags.
    weight_direction: bool,
    weight_energy: bool,
    weight_position: bool,
}

#[derive(Default)]
enum Direction {
    #[default]
    None,
    Point { azimuth: f64, elevation: f64 },
    SolidAngle { azimuth: [f64; 2], elevation: [f64; 2], sin_elevation: [f64; 2] },
}

#[derive(Default)]
enum Energy {
    #[default]
    None,
    Point(f64),
    PowerLaw { energy_min: f64, energy_max: f64, exponent: f64 },
}

#[derive(Default)]
enum Position {
    Inside { geobox: Py<GeoBox>, limit: Option<f64> },
    #[default]
    None,
    Point { latitude: f64, longitude: f64, altitude: f64 },
    Target { geobox: Py<GeoBox> },
}

#[derive(FromPyObject)]
enum Limit {
    Bool(bool),
    Float(f64),
}

impl Limit {
    const DEFAULT_LIMIT: f64 = 3.0;
}

impl From<Limit> for Option<f64> {
    fn from(value: Limit) -> Self {
        match value {
            Limit::Bool(b) => if b { Some(Limit::DEFAULT_LIMIT) } else { None },
            Limit::Float(v) => Some(v),
        }
    }
}

#[pymethods]
impl ParticlesGenerator {
    const PI: f64 = ::std::f64::consts::PI;

    #[new]
    pub fn new<'py>(
        py: Python<'py>,
        geometry: Option<Bound<'py, Geometry>>,
        random: Option<Bound<'py, Random>>,
        weight: Option<bool>,
    ) -> PyResult<Self> {
        let geometry = match geometry {
            None => Geometry::new(py, None)?,
            Some(geometry) => geometry.unbind(),
        };
        let random = match random {
            None => Py::new(py, Random::new(None, None)?)?,
            Some(random) => random.unbind(),
        };
        let weight = weight.unwrap_or(true);
        let generator = Self {
            geometry, random,
            direction: Direction::default(),
            energy: Energy::default(),
            pid: None,
            position: Position::default(),
            weight_direction: weight,
            weight_energy: weight,
            weight_position: weight,
        };
        Ok(generator)
    }

    /// Set Monte Carlo particles direction.
    fn direction<'py>(
        slf: Bound<'py, Self>,
        azimuth: Option<f64>,
        elevation: Option<f64>,
    ) -> PyResult<Bound<'py, Self>> {
        let azimuth = azimuth.unwrap_or(0.0);
        let elevation = elevation.unwrap_or(0.0);
        let mut generator = slf.borrow_mut();
        generator.direction = Direction::Point { azimuth, elevation };
        Ok(slf)
    }

    /// Set Monte Carlo particles energy.
    fn energy<'py>(
        slf: Bound<'py, Self>,
        value: f64,
    ) -> PyResult<Bound<'py, Self>> {
        let mut generator = slf.borrow_mut();
        generator.energy = Energy::Point(value);
        Ok(slf)
    }

    /// Generate Monte Carlo particles according to current selection(s).
    #[pyo3(signature=(shape=None, /))]
    fn generate<'py>(&self, py: Python<'py>, shape: Option<ShapeArg>) -> PyResult<PyObject> {
        // Check configuration.
        let is_rejection = match self.position {
            Position::Inside { limit, .. } => {
                if limit.is_some() {
                    if self.weight_direction != self.weight_position {
                        let err = Error::new(ValueError)
                            .what("configuration")
                            .why("'direction' weight conflicts with 'inside' mode");
                        return Err(err.to_err())
                    }
                    if self.weight_energy != self.weight_position {
                        let err = Error::new(ValueError)
                            .what("configuration")
                            .why("'energy' weight conflicts with 'inside' mode");
                        return Err(err.to_err())
                    }
                }
                true
            },
            Position::Target { .. } => {
                if self.weight_direction != self.weight_position {
                    let err = Error::new(ValueError)
                        .what("configuration")
                        .why("'direction' weight conflicts with 'target' mode");
                    return Err(err.to_err())
                }
                match self.direction {
                    Direction::None => false,
                    Direction::Point { .. } => false,
                    Direction::SolidAngle { .. } => true,
                }
            },
            _ => false,
        };

        // Create particles container.
        let shape: Vec<usize> = match shape {
            Some(shape) => shape.into(),
            None => Vec::new(),
        };
        let array = PyArray::<Particle>::zeros(py, &shape)?;
        let particles = unsafe { array.slice_mut()? };

        // Bind geometry etc.
        let mut geometry = self.geometry.bind(py).borrow_mut();
        let tracer = Tracer::new(&mut geometry, Mode::Merge)?;
        let mut random = self.random.bind(py).borrow_mut();

        // Prepare any box generator.
        let (bg, pb) = match &self.position {
            Position::Inside { geobox, .. } => {
                let bg = BoxGenerator::new(&geobox.bind(py).borrow());
                (Some(bg), None)
            },
            Position::Target { geobox } => match self.direction {
                Direction::Point { azimuth, elevation } => {
                    let direction = HorizontalCoordinates { azimuth, elevation };
                    let pb = ProjectedBox::new(&geobox.bind(py).borrow(), &direction);
                    (None, Some(pb))
                },
                _ => {
                    let bg = BoxGenerator::new(&geobox.bind(py).borrow());
                    (Some(bg), None)
                },
            },
            _ => (None, None),
        };

        // Loop over events.
        let mut trials: usize = 0;
        for particle in particles.iter_mut() {
            particle.pid = match self.pid {
                None => match self.position {
                    Position::Target { .. } => 16,
                    _ => 15,
                },
                Some(pid) => pid,
            };
            loop {
                trials += 1;
                if (trials % 1000) == 0 && ctrlc_catched() {
                    return Err(Error::new(KeyboardInterrupt).to_err())
                }

                particle.weight = 1.0;

                let (direction, energy, position) = match self.position {
                    Position::Inside { .. } => self.generate_inside(
                        bg.as_ref().unwrap(), &tracer, &mut random, particle
                    ),
                    Position::None => (false, false, true),
                    Position::Point { latitude, longitude, altitude } => {
                        particle.latitude = latitude;
                        particle.longitude = longitude;
                        particle.altitude = altitude;
                        (false, false, true)
                    },
                    Position::Target { .. } => match bg.as_ref() {
                        Some(bg) => self.generate_target(bg, &mut random, particle),
                        None => self.generate_through(pb.as_ref().unwrap(), &mut random, particle),
                    },
                };
                if !position {
                    continue;
                }

                if !direction {
                    self.generate_direction(&mut random, particle);
                }

                if !energy {
                    self.generate_energy(&mut random, particle);
                }
                break;
            }
        }

        // Return result.
        let array: &PyAny = array;
        let array: PyObject = array.into();
        let result = if is_rejection {
            static RESULT: NamedTuple<2> = NamedTuple::new(
                "Result", ["particles", "size"]);
            RESULT.instance(py, (array, trials))?.unbind()
        } else {
            array
        };
        Ok(result)
    }

    /// Select Monte Carlo particles inside a box.
    fn inside<'py>(
        slf: Bound<'py, Self>,
        r#box: &Bound<'py, GeoBox>,
        limit: Option<Limit>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let limit: Option<f64> = limit.and_then(|v| v.into());
        let mut generator = slf.borrow_mut();
        if let Some(weight) = weight {
            generator.weight_position = weight;
        }
        generator.position = Position::Inside { geobox: r#box.clone().unbind(), limit };
        Ok(slf)
    }

    /// Set Monte Carlo particles type.
    fn pid<'py>(
        slf: Bound<'py, Self>,
        value: c_int,
    ) -> PyResult<Bound<'py, Self>> {
        let mut generator = slf.borrow_mut();
        generator.pid = Some(value);
        Ok(slf)
    }

    /// Set Monte Carlo particles position.
    fn position<'py>(
        slf: Bound<'py, Self>,
        latitude: Option<f64>,
        longitude: Option<f64>,
        altitude: Option<f64>,
    ) -> PyResult<Bound<'py, Self>> {
        let latitude = latitude.unwrap_or(0.0);
        let longitude = longitude.unwrap_or(0.0);
        let altitude = altitude.unwrap_or(0.0);
        let mut generator = slf.borrow_mut();
        generator.position = Position::Point { latitude, longitude, altitude };
        Ok(slf)
    }

    /// Select particles energy according to a power-law.
    fn powerlaw<'py>(
        slf: Bound<'py, Self>,
        energy_min: f64,
        energy_max: f64,
        exponent: Option<f64>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        if energy_min >= energy_max || energy_min <= 0.0 {
            let why = "expected energy_max > energy_min > 0.0";
            let err = Error::new(ValueError).what("powerlaw").why(why);
            return Err(err.to_err());
        }
        let exponent = exponent.unwrap_or(-1.0);
        let mut generator = slf.borrow_mut();
        if let Some(weight) = weight {
            generator.weight_energy = weight;
        }
        generator.energy = Energy::PowerLaw { energy_min, energy_max, exponent };
        Ok(slf)
    }

    /// Select particles direction uniformly over a solid-angle.
    fn solid_angle<'py>(
        slf: Bound<'py, Self>,
        azimuth: Option<[f64; 2]>,
        elevation: Option<[f64; 2]>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let azimuth = azimuth.map(|[a, b]| if a <= b { [ a, b ] } else { [ b, a ] });
        let elevation = elevation.map(|[a, b]| if a <= b { [ a, b ] } else { [ b, a ] });
        let sin_elevation = match elevation {
            None => [-1.0, 1.0],
            Some([el0, el1]) => {
                check_angle(el0, -90.0, 90.0, "elevation")?;
                check_angle(el1, -90.0, 90.0, "elevation")?;
                let el0 = el0 * Self::RAD;
                let el1 = el1 * Self::RAD;
                [el0.sin(), el1.sin()]
            },
        };
        let elevation = elevation.unwrap_or([-90.0, 90.0]);
        let azimuth = match azimuth {
            None => [-180.0, 180.0],
            Some([az0, az1]) => {
                check_angle(az0, -180.0, 180.0, "azimuth")?;
                check_angle(az1, -180.0, 180.0, "azimuth")?;
                [az0, az1]
            },
        };

        let mut generator = slf.borrow_mut();
        if let Some(weight) = weight {
            generator.weight_direction = weight;
        }
        generator.direction = Direction::SolidAngle { azimuth, elevation, sin_elevation };
        Ok(slf)
    }

    /// Select Monte Carlo particles targeting a box volume.
    fn target<'py>(
        slf: Bound<'py, Self>,
        r#box: Bound<'py, GeoBox>,
        weight: Option<bool>,
    ) -> PyResult<Bound<'py, Self>> {
        let mut generator = slf.borrow_mut();
        if let Some(weight) = weight {
            generator.weight_position = weight;
        }
        generator.position = Position::Target { geobox: r#box.clone().unbind() };
        Ok(slf)
    }
}

impl ParticlesGenerator {
    const RAD: f64 = std::f64::consts::PI / 180.0;

    fn generate_direction(&self, random: &mut Random, particle: &mut Particle) {
        let (azimuth, elevation, weight) = self.direction.generate(random);
        particle.azimuth = azimuth;
        particle.elevation = elevation;
        if self.weight_direction {
            particle.weight *= weight;
        }
    }

    fn generate_energy(&self, random: &mut Random, particle: &mut Particle) {
        let (energy, weight) = self.energy.generate(random);
        particle.energy = energy;
        if self.weight_energy {
            particle.weight *= weight;
        }
    }

    fn generate_inside(
        &self,
        bg: &BoxGenerator,
        tracer: &Tracer,
        random: &mut Random,
        particle: &mut Particle
    ) -> (bool, bool, bool) {
        let geodetic = bg.generate_inside(random);

        // Check that the location is in the air.
        let medium = tracer.medium(&geodetic);
        if !medium.is_atmosphere() {
            return (false, false, false)
        }

        let Position::Inside { limit, .. } = self.position else { unreachable!() };

        let mut is_direction = false;
        let mut is_energy = false;
        if let Some(limit) = limit {
            let (azimuth, elevation, d_weight) = self.direction.generate(random);
            let (energy, e_weight) = self.energy.generate(random);

            // Check the path length along the track origin.
            let horizontal = HorizontalCoordinates {
                azimuth: azimuth + 180.0,
                elevation: -elevation
            };
            let decay_length = {
                const TAU_CTAU: f64 = 8.718E-05;
                const TAU_MASS: f64 = 1.77682;
                energy * TAU_CTAU / TAU_MASS
            };
            let distance_limit = Some(limit * decay_length);
            let (_, distance, next_medium) = tracer.trace(&geodetic, &horizontal, distance_limit);
            let p = if next_medium.is_atmosphere() {
                (-limit).exp()
            } else {
                (-distance / decay_length).exp()
            };
            if random.open01() > p {
                return (true, true, false);
            } else {
                particle.energy = energy;
                particle.azimuth = azimuth;
                particle.elevation = elevation;
                if self.weight_position {
                    particle.weight *= d_weight * e_weight / p;
                }
                is_direction = true;
                is_energy = true;
            }
        }

        particle.latitude = geodetic.latitude;
        particle.longitude = geodetic.longitude;
        particle.altitude = geodetic.altitude;

        if self.weight_position {
            particle.weight *= bg.volume();
        }
        (is_direction, is_energy, true)
    }

    fn generate_target(
        &self,
        bg: &BoxGenerator,
        random: &mut Random,
        particle: &mut Particle
    ) -> (bool, bool, bool) {
        let (r, u, geodetic, horizontal) = bg.generate_onto(random);

        match self.direction {
            Direction::SolidAngle { azimuth, elevation, .. } => {
                // Find the second (outgoing) intersection (using the local frame).
                // Ref: https://iquilezles.org/articles/intersectors/ (axis-aligned box).
                let m = [ 1.0 / u[0], 1.0 / u[1], 1.0 / u[2] ];
                let n = [ m[0] * r[0], m[1] * r[1], m[2] * r[2] ];
                let k = [
                    0.5 * m[0].abs() * bg.size()[0],
                    0.5 * m[1].abs() * bg.size()[1],
                    0.5 * m[2].abs() * bg.size()[2],
                ];
                let t2 = [ -n[0] + k[0], -n[1] + k[1], -n[2] + k[2] ];
                let tf = t2[0].min(t2[1]).min(t2[2]);
                let r2 = [
                    r[0] + tf * u[0],
                    r[1] + tf * u[1],
                    r[2] + tf * u[2],
                ];

                let geodetic2 = bg.frame().to_geodetic(&r2);
                let horizontal2 = bg.frame().to_horizontal(&u, &geodetic2);

                // Check the azimuth and elevation values (inclusively).
                let az0 = horizontal.azimuth.min(horizontal2.azimuth);
                let az1 = horizontal.azimuth.max(horizontal2.azimuth);
                let el0 = horizontal.elevation.min(horizontal2.elevation);
                let el1 = horizontal.elevation.max(horizontal2.elevation);

                let [az_min, az_max] = azimuth;
                let [el_min, el_max] = elevation;
                if (az1 < az_min) || (az0 > az_max) || (el1 < el_min) || (el0 > el_max) {
                    return (false, false, false)
                }
            },
            Direction::None => (),
            Direction::Point { .. } => unreachable!(),
        }

        particle.latitude = geodetic.latitude;
        particle.longitude = geodetic.longitude;
        particle.altitude = geodetic.altitude;
        particle.azimuth = horizontal.azimuth;
        particle.elevation = horizontal.elevation;

        if self.weight_position {
            particle.weight *= bg.surface() * Self::PI;
        }
        (true, false, true)
    }

    fn generate_through(
        &self,
        pb: &ProjectedBox,
        random: &mut Random,
        particle: &mut Particle
    ) -> (bool, bool, bool) {
        let (geodetic, horizontal) = pb.generate_inside(random);

        particle.latitude = geodetic.latitude;
        particle.longitude = geodetic.longitude;
        particle.altitude = geodetic.altitude;
        particle.azimuth = horizontal.azimuth;
        particle.elevation = horizontal.elevation;

        if self.weight_position {
            particle.weight *= pb.surface();
        }
        (true, false, true)
    }
}

fn check_angle(value: f64, min: f64, max: f64, what: &str) -> PyResult<()> {
    if (value < min) || (value > max) {
        let why = format!(
            "expected a value in [{}, {}], found {}",
            min,
            max,
            value,
        );
        let err = Error::new(ValueError)
            .what(what)
            .why(&why);
        Err(err.to_err())
    } else {
        Ok(())
    }
}

impl Direction {
    const DEG: f64 = 180.0 / std::f64::consts::PI;
    const RAD: f64 = std::f64::consts::PI / 180.0;

    fn generate(&self, random: &mut Random) -> (f64, f64, f64) {
        match self {
            Self::None => (0.0, 0.0, 1.0),
            Self::Point { azimuth, elevation } => (*azimuth, *elevation, 1.0),
            Self::SolidAngle { .. } => self.generate_solid_angle(random),
        }
    }

    fn generate_solid_angle(&self, random: &mut Random) -> (f64, f64, f64) {
        let Direction::SolidAngle { azimuth, sin_elevation, .. } = *self else { unreachable!() };
        let [az0, az1] = azimuth;
        let [sin_el0, sin_el1] = sin_elevation;
        let sin_el = (sin_el1 - sin_el0) * random.open01() + sin_el0;
        let el = sin_el.asin() * Self::DEG;
        let az = (az1 - az0) * random.open01() + az0;
        let solid_angle = (sin_el1 - sin_el0).abs() * (az1 - az0).abs() * Self::RAD;
        (az, el, solid_angle)
    }
}

impl Energy {
    fn generate(&self, random: &mut Random) -> (f64, f64) {
        match self {
            Self::None => (1E+09, 1.0),
            Self::Point(value) => (*value, 1.0),
            Self::PowerLaw { .. } => self.generate_powerlaw(random),
        }
    }

    fn generate_powerlaw(&self, random: &mut Random) -> (f64, f64) {
        let Energy::PowerLaw { energy_min, energy_max, exponent } = *self else { unreachable!() };
        let (energy, weight) = if exponent == -1.0 {
            let lne = (energy_max / energy_min).ln();
            let energy = energy_min * (random.open01() * lne).exp();
            (energy, energy * lne)
        } else if exponent == 0.0 {
            let de = energy_max - energy_min;
            let energy = de * random.open01() + energy_min;
            (energy, de)
        } else {
            let a = exponent + 1.0;
            let b = energy_min.powf(a);
            let de = energy_max.powf(a) - b;
            let energy = (de * random.open01() + b).powf(1.0 / a);
            let weight = de / (a * energy.powf(exponent));
            (energy, weight)
        };
        let energy = energy.clamp(energy_min, energy_max);
        (energy, weight)
    }
}
