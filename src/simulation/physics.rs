use crate::bindings::danton;
use crate::utils::convert::{Bremsstrahlung, Dis, PairProduction, Pdf, Photonuclear};
use crate::utils::error::to_result;
use pyo3::prelude::*;
use ::std::ffi::CString;
use ::std::sync::atomic::{AtomicUsize, Ordering};


// XXX Reset / modify physics.

static INSTANCES: AtomicUsize = AtomicUsize::new(0);
static CURRENT: AtomicUsize = AtomicUsize::new(0);

#[derive(Default)]
#[pyclass(module="danton")]
pub struct Physics {
    #[pyo3(get, set)]
    bremsstrahlung: Bremsstrahlung,
    #[pyo3(get, set)]
    pair_production: PairProduction,
    #[pyo3(get, set)]
    photonuclear: Photonuclear,
    #[pyo3(get, set)]
    dis: Dis,
    #[pyo3(get, set)]
    pdf: Option<Pdf>,

    instance: usize,
    modified: bool,
}

#[pymethods]
impl Physics {
    #[new]
    pub fn new() -> Self {
        let mut physics = Self::default();
        physics.instance = INSTANCES.fetch_add(1, Ordering::SeqCst);
        physics.modified = true;
        physics
    }
}

impl Physics {
    pub fn apply(&mut self) -> PyResult<()> {
        let previous = CURRENT.swap(self.instance, Ordering::SeqCst);
        if (previous != self.instance) || self.modified {
            let set = |key: &str, value: &str| -> PyResult<()> {
                let key = CString::new(key)?;
                let value = CString::new(value)?;
                to_result(unsafe {
                    danton::physics_set(key.as_ptr(), value.as_ptr())
                }, None)
            };
            set("bremsstrahlung", self.bremsstrahlung.into())?;
            set("pair-production", self.pair_production.into())?;
            set("photonuclear", self.photonuclear.into())?;
            set("DIS", self.dis.as_ref())?;
            if let Some(pdf) = self.pdf.as_ref() {
                set("PDF", pdf.as_ref())?;
            }
            self.modified = false;
        }
        Ok(())
    }
}
