use crate::bindings::danton;
use crate::simulation::Simulation;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::SystemError;
use crate::utils::numpy::{PyArray, ShapeArg};
use getrandom::getrandom;
use pyo3::prelude::*;
use rand::Rng;
use rand::distributions::Open01;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;
use ::std::ptr::null_mut;


// ===============================================================================================
//
// Generator interface.
//
// ===============================================================================================

/// A Pseudo-Random Numbers Generator (PRNG).
#[pyclass(module = "danton")]
pub struct Random {
    rng: Pcg64Mcg,
    /// Prng stream index.
    #[pyo3(get)]
    index: u128,
    /// Prng initial seed.
    #[pyo3(get)]
    seed: u128,
}

#[pymethods]
impl Random {
    #[new]
    pub fn new(seed: Option<u128>) -> PyResult<Self> {
        let rng = Pcg64Mcg::new(0xCAFEF00DD15EA5E5);
        let mut random = Self { rng, seed: 0, index: 0 };
        random.initialise(seed)?;
        Ok(random)
    }

    #[setter]
    fn set_index(&mut self, index: Option<u128>) -> PyResult<()> {
        match index {
            None => self.initialise(Some(self.seed))?,
            Some(index) => {
                let delta: u128 = index.wrapping_sub(self.index);
                self.rng.advance(delta);
                self.index = index;
            },
        }
        Ok(())
    }

    #[setter]
    fn set_seed(&mut self, seed: Option<u128>) -> PyResult<()> {
        self.initialise(seed)
    }

    /// Generate pseudo-random number(s) uniformly distributed over (0,1).
    fn uniform01(
        &mut self,
        py: Python,
        shape: Option<ShapeArg>,
    ) -> PyResult<PyObject> {
        match shape {
            None => {
                let value = self.open01();
                Ok(value.into_py(py))
            },
            Some(shape) => {
                let shape: Vec<usize> = shape.into();
                let n = shape.iter().product();
                let iter = (0..n).map(|_| self.open01());
                let array: &PyAny = PyArray::<f64>::from_iter(py, &shape, iter)?;
                Ok(array.into())
            },
        }
    }
}

impl Random {
    fn initialise(&mut self, seed: Option<u128>) -> PyResult<()> {
        match seed {
            None => {
                let mut seed = [0_u8; 16];
                getrandom(&mut seed)
                    .map_err(|_| Error::new(SystemError)
                        .what("random")
                        .why("could not seed random engine")
                        .to_err()
                )?;
                self.rng = Pcg64Mcg::from_seed(seed);
                self.seed = u128::from_ne_bytes(seed);
            },
            Some(seed) => {
                self.seed = seed;
                let seed = u128::to_ne_bytes(seed);
                self.rng = Pcg64Mcg::from_seed(seed);
            },
        }
        self.index = 0;
        Ok(())
    }

    pub fn new_context<'a>(&'a mut self, simulation: &Simulation) -> RandomContext<'a> {
        RandomContext { context: simulation.context, random: self }
    }

    fn open01(&mut self) -> f64 {
        self.index += 1;
        self.rng.sample::<f64, Open01>(Open01)
    }
}


// ===============================================================================================
//
// C interface.
//
// ===============================================================================================

#[repr(C)]
pub struct RandomContext<'a> {
    context: *mut danton::Context,
    random: &'a mut Random,
}

#[no_mangle]
pub extern "C" fn danton_context_random_open01(context: *mut RandomContext) -> f64 {
    let context = unsafe { &mut *context };
    context.random.open01()
}

impl<'a> Drop for RandomContext<'a> {
    fn drop(&mut self) {
        unsafe { danton::context_random_set(self.context, null_mut()) };
    }
}
