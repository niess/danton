use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::convert::{Geodesic, Mode};
use crate::utils::export::Export;
use derive_more::{AsMut, AsRef, From};
use pyo3::prelude::*;
use ::std::ffi::{c_int, CString, c_uint};
use ::std::pin::Pin;


#[repr(C)]
pub struct Stepper {
    base: danton::RunAction,
    pub event: usize,
    pub mode: Mode,
    pub geodesic: Geodesic,
    pub ocean: bool,
    steps: Option<Vec<Step>>,
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct StepsExport (Export<Step>);

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Step {
    event: usize,
    pid: c_int,
    energy: f64,
    latitude: f64,
    longitude: f64,
    altitude: f64,
    azimuth: f64,
    elevation: f64,
    medium: [u8; 16],
}

impl From<Particle> for Step {
    fn from(value: Particle) -> Self {
        Self {
            event: 0,
            pid: value.pid,
            energy: value.energy,
            latitude: value.latitude,
            longitude: value.longitude,
            altitude: value.altitude,
            azimuth: value.azimuth,
            elevation: value.elevation,
            medium: [0; 16],
        }
    }
}

impl Stepper {
    pub fn new() -> Pin<Box<Self>> {
        let base = danton::RunAction {
            call: Some(Self::call),
        };
        let stepper = Self {
            base,
            event: 0,
            mode: Mode::Backward,
            geodesic: Geodesic::Prem,
            ocean: true,
            steps: None,
        };
        Box::pin(stepper)
    }

    pub fn base(stepper: &mut Pin<Box<Self>>) -> &mut danton::RunAction {
        &mut stepper.base
    }

    pub fn clear(&mut self) {
        self.steps = None;
    }

    pub fn export(&mut self, py: Python) -> PyResult<PyObject> {
        let steps = match self.steps.take() {
            None => py.None(),
            Some(steps) => Export::export::<StepsExport>(py, steps)?,
        };
        return Ok(steps)
    }

    unsafe extern "C" fn call(
        _context: *mut danton::Context,
        run_action: *mut danton::RunAction,
        event: c_uint,
        medium: c_int,
        state: *mut danton::State,
    ) -> c_int {
        if event != danton::STEP {
            return danton::SUCCESS;
        }
        let stepper = unsafe { &mut *(run_action as *mut Self) };
        if !state.is_null() {
            let state = unsafe { &*state };
            let step = {
                let particle: Particle = (state, stepper.geodesic).into();
                let mut step: Step = particle.into();
                step.event = stepper.event;
                if let Mode::Grammage = stepper.mode {
                    step.pid = 0;
                    step.energy = 0.0;
                    step.elevation = -step.elevation;
                    if step.elevation.abs() != 90.0 {
                        if step.azimuth > 0.0 {
                            step.azimuth -= 180.0;
                        } else {
                            step.azimuth += 180.0;
                        }
                    }
                }
                let ground = "UpperCrust";
                let medium = match medium {
                    0 => "InnerCore",
                    1 => "OuterCore",
                    2 => "LowerMantle",
                    3 => "Mantle2",
                    4 => "Mantle1",
                    5 => "Mantle0",
                    6 => "UpperMantle",
                    7 => "LowerCrust",
                    8 => ground,
                    9 => if stepper.ocean { "Ocean" } else { ground },
                    10 => "Troposphere0",
                    11 => "Troposhpere1",
                    12 => "Stratosphere",
                    13 => "Mesosphere",
                    14 => "Exosphere",
                    -1 => "Exit",
                    100 => "Topography",
                    _ => "Unknown",
                };
                let medium = CString::new(medium).unwrap();
                let bytes = medium.as_bytes();
                for (i, bi) in bytes.iter().enumerate() {
                    if i >= step.medium.len() - 1 {
                        break;
                    }
                    step.medium[i] = *bi;
                }
                step
            };
            match stepper.steps.as_mut() {
                None => stepper.steps = Some(vec![step]),
                Some(steps) => steps.push(step),
            };
        }

        danton::SUCCESS
    }
}
