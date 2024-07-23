use crate::bindings::danton;
use crate::utils::convert::Mode;
use crate::utils::error::to_result;
use crate::utils::numpy::PyArray;
use pyo3::prelude::*;
use ::std::pin::Pin;

pub mod geometry;
pub mod particles;
mod primary;
pub mod recorder;
mod sampler;


#[pyclass(module="danton")]
pub struct Simulation {
    #[pyo3(get, set)]
    geometry: Py<geometry::Geometry>,
    context: *mut danton::Context,
    recorder: Pin<Box<recorder::Recorder>>,
    primaries: [Pin<Box<danton::Primary>>; 6],
}

unsafe impl Send for Simulation {}

#[pymethods]
impl Simulation {
    #[new]
    fn new(py: Python) -> PyResult<Self> {
        let geometry = Py::new(py, geometry::Geometry::new())?;
        let mut recorder = recorder::Recorder::new();
        let mut primaries = core::array::from_fn(|_| danton::Primary::new(1E+03, 1E+12));
        let context = unsafe {
            let context = danton::context_create();
            (*context).sampler = danton::sampler_create();
            let recorder: &mut danton::Recorder = recorder::Recorder::base(&mut recorder);
            (*context).recorder = recorder;
            for i in 0..primaries.len() {
                (*context).primary[i] = &mut *primaries[i];
            }
            context
        };
        let simulation = Self { geometry, context, recorder, primaries };
        Ok(simulation)
    }

    /// Flag controlling the decay of tau final states.
    #[getter]
    fn get_decay(&self, py: Python) -> bool {
        self.context(py).decay != 0
    }

    #[setter]
    fn set_decay(&mut self, py: Python, value: bool) {
        self.context_mut(py).decay = if value { 1 } else { 0 };
    }

    /// Flag for controlling the tranverse transport.
    #[getter]
    fn get_longitudinal(&self, py: Python) -> bool {
        self.context(py).longitudinal != 0
    }

    #[setter]
    fn set_longitudinal(&mut self, py: Python, value: bool) {
        self.context_mut(py).longitudinal = if value { 1 } else { 0 };
    }

    /// The run mode.
    #[getter]
    fn get_mode(&self, py: Python) -> Mode {
        self.context(py).mode.into()
    }

    #[setter]
    fn set_mode(&mut self, py: Python, value: Mode) {
        self.context_mut(py).mode = value.into();
    }

    fn run(
        &mut self,
        py: Python,
        particles: &PyArray<particles::Particle>,
    ) -> PyResult<PyObject> {
        self.recorder.geodesic = {
            let mut geometry = self.geometry.bind(py).borrow_mut();
            geometry.apply()?;
            geometry.geodesic
        };
        let n = particles.size();
        for i in 0..n {
            self.recorder.event = i;
            let particle = particles.get(i)?;
            self.sampler_mut(py).set(&particle)?;
            to_result(
                unsafe { danton::context_run(self.context, 1, 0) },
                Some(self.context),
            )?;
        }
        self.recorder.export(py)
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
