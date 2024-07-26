use crate::bindings::danton;
use crate::utils::convert::Mode;
use crate::utils::error::{Error, to_result};
use crate::utils::error::ErrorKind::NotImplementedError;
use pyo3::prelude::*;
use ::std::pin::Pin;
use ::std::ptr::null_mut;

pub mod geometry;
pub mod particles;
pub mod physics;
mod primary;
pub mod recorder;
mod sampler;
pub mod stepper;


#[pyclass(module="danton")]
pub struct Simulation {
    #[pyo3(get, set)]
    geometry: Py<geometry::Geometry>,
    #[pyo3(get, set)]
    physics: Py<physics::Physics>,
    context: *mut danton::Context,
    recorder: Pin<Box<recorder::Recorder>>,
    primaries: [Pin<Box<danton::Primary>>; 6],
    stepper: Option<Pin<Box<stepper::Stepper>>>,
}

unsafe impl Send for Simulation {}

#[pymethods]
impl Simulation {
    #[new]
    fn new(py: Python) -> PyResult<Self> {
        let geometry = Py::new(py, geometry::Geometry::new())?;
        let physics = Py::new(py, physics::Physics::new())?;
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
        let simulation = Self { geometry, physics, context, recorder, primaries, stepper: None };
        Ok(simulation)
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

    /// Run a Danton Monte Carlo simulation.
    fn run<'py>(&mut self, particles: &Bound<'py, PyAny>) -> PyResult<PyObject> {
        // Configure samplers etc.
        let py = particles.py();
        let mode = self.get_mode(py);
        let decay = self.get_tau_decays(py);
        self.physics.bind(py).borrow_mut().apply()?;
        let (geodesic, ocean) = {
            let mut geometry = self.geometry.bind(py).borrow_mut();
            geometry.apply()?;
            (geometry.geodesic, geometry.ocean)
        };
        self.recorder.geodesic = geodesic;
        self.recorder.mode = mode;
        self.recorder.decay = decay;
        if let Some(stepper) = self.stepper.as_mut() {
            stepper.mode = mode;
            stepper.geodesic = geodesic;
            stepper.ocean = ocean;
        }

        // Loop over events.
        let particles = match mode {
            Mode::Backward | Mode::Forward => particles::ParticlesIterator::new(particles)?,
            Mode::Grammage => particles::ParticlesIterator::coordinates(particles)?,
        };
        for (i, particle) in particles.enumerate() {
            let particle = particle?;

            self.recorder.event = i;
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
