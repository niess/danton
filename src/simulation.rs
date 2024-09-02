use crate::bindings::danton;
use crate::utils::convert::Mode;
use crate::utils::error::{ctrlc_catched, Error, to_result};
use crate::utils::error::ErrorKind::{KeyboardInterrupt, NotImplementedError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString};
use ::std::pin::Pin;
use ::std::ptr::null_mut;

pub mod geobox;
pub mod geometry;
pub mod particles;
pub mod physics;
mod primary;
pub mod random;
pub mod recorder;
mod sampler;
pub mod stepper;


#[pyclass(module="danton")]
pub struct Simulation {
    /// The Monte Carlo geometry.
    #[pyo3(get, set)]
    geometry: Py<geometry::Geometry>,
    /// The Monte Carlo physics models.
    #[pyo3(get)]
    physics: Py<physics::Physics>,
    /// The Monte Carlo random engine.
    #[pyo3(get, set)]
    random: Py<random::Random>,
    context: *mut danton::Context,
    recorder: Pin<Box<recorder::Recorder>>,
    primaries: [Pin<Box<danton::Primary>>; 6],
    stepper: Option<Pin<Box<stepper::Stepper>>>,
}

unsafe impl Send for Simulation {}

#[pymethods]
impl Simulation {
    #[pyo3(signature=(**kwargs))]
    #[new]
    pub fn new<'py>(
        py: Python<'py>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Py<Self>> {
        let (geometry_kwargs, physics_kwargs) = match kwargs {
            None => (None, None),
            Some(kwargs) => {
                let extract = |fields: &[&str]| -> PyResult<Option<Bound<'py, PyDict>>> {
                    let mut result = None;
                    for field in fields {
                        if let Some(value) = kwargs.get_item(field)? {
                            result
                                .get_or_insert_with(|| {
                                    PyDict::new_bound(py)
                                })
                                .set_item(field, value)?;
                            kwargs.del_item(field)?;
                        }
                    }
                    Ok(result)
                };
                let geometry_kwargs = extract(
                    &["density", "geoid", "materials", "ocean", "topography"]
                )?;
                let physics_kwargs = extract(
                    &["bremsstrahlung", "dis", "pair_production", "pdf", "photonuclear"]
                )?;
                (geometry_kwargs, physics_kwargs)
            },
        };

        let geometry = geometry::Geometry::new(py, geometry_kwargs.as_ref())?;
        let physics = physics::Physics::new(py, physics_kwargs.as_ref())?;
        let random = Py::new(py, random::Random::new(None)?)?;
        let mut recorder = recorder::Recorder::new();
        let mut primaries = core::array::from_fn(|_| danton::Primary::new());
        let context = unsafe {
            let context = danton::context_create();
            (*context).sampler = danton::sampler_create();
            let recorder = recorder::Recorder::base(&mut recorder);
            (*context).recorder = recorder;
            for i in 0..primaries.len() {
                (*context).primary[i] = &mut *primaries[i];
            }
            context
        };
        let simulation = Self {
            geometry, physics, random, context, recorder, primaries, stepper: None
        };
        let simulation = Bound::new(py, simulation)?;

        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key: Bound<PyString> = key.extract()?;
                simulation.setattr(key, value)?
            }
        }

        Ok(simulation.unbind())
    }

    /// Flag to control the sampling of tau decays.
    #[getter]
    fn get_tau_decays(&self, py: Python) -> bool {
        self.context(py).decay != 0
    }

    #[setter]
    fn set_tau_decays(&mut self, py: Python, value: bool) {
        self.context_mut(py).decay = if value { 1 } else { 0 };
    }

    /// Flag to control the longitudinal approximation.
    #[getter]
    fn get_longitudinal(&self, py: Python) -> bool {
        self.context(py).longitudinal != 0
    }

    #[setter]
    fn set_longitudinal(&mut self, py: Python, value: bool) {
        self.context_mut(py).longitudinal = if value { 1 } else { 0 };
    }

    /// The Monte Carlo simulation mode.
    #[getter]
    fn get_mode(&self, py: Python) -> Mode {
        self.context(py).mode.into()
    }

    #[setter]
    fn set_mode(&mut self, py: Python, value: Mode) {
        if let Mode::Backward = value {
            for primary in self.primaries.iter_mut() {
                primary.configure(None)
            }
        }
        self.context_mut(py).mode = value.into();
    }

    #[setter]
    fn set_physics<'py>(&mut self, value: Bound<'py, physics::Physics>) {
        unsafe { danton::context_reset(self.context) };
        self.physics = value.clone().unbind();
    }

    /// Flag to control the recording of Monte Carlo steps.
    #[getter]
    fn get_record_steps(&self) -> bool {
        self.stepper.is_some()
    }

    #[setter]
    fn set_record_steps(&mut self, py: Python, value: bool) {
        if value != self.stepper.is_some() {
            match value {
                true => {
                    let mut stepper = stepper::Stepper::new();
                    self.context_mut(py).run_action = stepper::Stepper::base(&mut stepper);
                    self.stepper = Some(stepper);
                },
                false => {
                    self.context_mut(py).run_action = null_mut();
                    self.stepper = None;
                },
            }
        }
    }

    /// Create a local box.
    #[pyo3(name="r#box", signature=(size=None, *, latitude=None, longitude=None, altitude=None, declination=None))]
    fn geobox(
        &self,
        py: Python,
        size: Option<geobox::Size>,
        latitude: Option<f64>,
        longitude: Option<f64>,
        altitude: Option<f64>,
        declination: Option<f64>,
    ) -> geobox::GeoBox {
        let ellipsoid = Some(self.geometry.bind(py).borrow().geoid.into());
        geobox::GeoBox::new(size, latitude, longitude, altitude, declination, ellipsoid)
    }

    /// Generate Monte Carlo particles.
    fn particles(
        &self,
        py: Python,
        weight: Option<bool>,
    ) -> PyResult<particles::ParticlesGenerator> {
        let geometry = Some(self.geometry.bind(py).clone());
        let random = Some(self.random.bind(py).clone());
        particles::ParticlesGenerator::new(py, geometry, random, weight)
    }

    /// Run a Danton Monte Carlo simulation.
    fn run<'py>(&mut self, particles: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        // Configure physics, geometry, samplers etc.
        let py = particles.py();
        let mode = self.get_mode(py);
        let decay = self.get_tau_decays(py);
        let (geoid, ocean) = {
            let mut geometry = self.geometry.bind(py).borrow_mut();
            geometry.apply()?;
            let mut physics = self.physics.bind(py).borrow_mut();
            if physics.modified {
                unsafe { danton::context_reset(self.context) };
            }
            physics.apply(py, geometry.material.as_str())?;
            (geometry.geoid, geometry.ocean)
        };
        self.recorder.ellipsoid = geoid.into();
        self.recorder.mode = mode;
        self.recorder.decay = decay;
        if let Some(stepper) = self.stepper.as_mut() {
            stepper.mode = mode;
            stepper.ellipsoid = geoid.into();
            stepper.ocean = ocean;
        }

        // Set random context.
        let mut random = self.random.bind(py).borrow_mut();
        let mut random_context = random.new_context(self);
        unsafe { danton::context_random_set(self.context, &mut random_context) };

        // Loop over events.
        let particles = match mode {
            Mode::Backward | Mode::Forward => particles::ParticlesIterator::new(particles)?,
            Mode::Grammage => particles::ParticlesIterator::coordinates(particles)?,
        };
        for (i, particle) in particles.enumerate() {
            let particle = particle?;

            self.recorder.event = i;
            self.recorder.weight = particle.weight;
            self.recorder.random_index = random_context.index();
            if let Some(stepper) = self.stepper.as_mut() {
                stepper.event = i;
            }

            // Configure the particles sampler.
            self.sampler_mut(py).set(mode, decay, &particle)?;

            // Configure primaries.
            if let Mode::Forward = mode {
                let particle_index = particle.index()? as usize;
                if particle_index >= 6 {
                    let err = Error::new(NotImplementedError)
                        .what("pid")
                        .why("primary is not a neutrino");
                    return Err(err.into())
                }
                for i in 0..self.primaries.len() {
                    self.context_mut(py).primary[i] = if i == particle_index {
                        self.primaries[i].configure(Some(particle.energy));
                        &mut *self.primaries[i]
                    } else {
                        null_mut()
                    };
                }
            }

            if ctrlc_catched() {
                self.recorder.clear();
                self.stepper.as_mut()
                    .map(|stepper| stepper.clear());
                return Err(Error::new(KeyboardInterrupt).to_err())
            }

            // Run the Monte Carlo.
            to_result(
                unsafe { danton::context_run(self.context, 1, 0) },
                Some(self.context),
            )?;
        }

        // Export the collected results.
        let steps = self.stepper.as_mut()
            .map(|stepper| stepper.export(py))
            .transpose()?;
        let result = self.recorder.export(py, steps)?;
        Ok(result)
    }
}

impl Simulation {
    #[inline]
    fn context(&self, _py: Python) -> &danton::Context {
        unsafe { &*self.context }
    }

    #[inline]
    fn context_mut(&mut self, _py: Python) -> &mut danton::Context {
        unsafe { &mut *self.context }
    }

    #[inline]
    fn sampler_mut(&mut self, py: Python) -> &mut danton::Sampler {
        let context = self.context_mut(py);
        unsafe { &mut *context.sampler }
    }
}

impl Drop for Simulation {
    fn drop(&mut self) {
        unsafe {
            let context = &mut *self.context;
            danton::Sampler::destroy(&mut context.sampler);
            danton::context_destroy(&mut self.context);
        }
    }
}
