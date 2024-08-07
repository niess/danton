use crate::bindings::danton;
use crate::utils::convert::{Geodesic, Medium};
use crate::utils::coordinates::{GeodeticCoordinates, HorizontalCoordinates};
use crate::utils::error::to_result;
use pyo3::prelude::*;
use ::std::ffi::{c_int, CString, c_void};
use ::std::ptr::null;
use ::std::sync::atomic::{AtomicUsize, Ordering};


static INSTANCES: AtomicUsize = AtomicUsize::new(0);
static CURRENT: AtomicUsize = AtomicUsize::new(0);

#[pyclass(module="danton")]
pub struct Geometry {
    #[pyo3(get)]
    /// Reference geodesic for the sea level.
    pub geodesic: Geodesic,
    #[pyo3(get)]
    /// Topography elevation data.
    topography: Option<Topography>,
    #[pyo3(get)]
    /// The topography composition.
    pub material: String,
    #[pyo3(get)]
    /// The topography density.
    density: f64,
    #[pyo3(get)]
    /// Flag to enable/disable the ocean.
    pub ocean: bool,

    instance: usize,
    modified: bool,
}

#[derive(Clone, FromPyObject, PartialEq)]
enum Topography {
    #[pyo3(transparent, annotation = "str")]
    Dem(String),
    #[pyo3(transparent, annotation = "float")]
    Flat(f64),
}

impl IntoPy<PyObject> for Topography {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::Dem(path) => path.into_py(py),
            Self::Flat(z) => z.into_py(py),
        }
    }
}

#[pymethods]
impl Geometry {
    #[new]
    pub fn new() -> Self {
        Self {
            geodesic: Geodesic::default(),
            topography: None,
            material: "Rock".to_string(),
            density: 2.65E+03,
            ocean: true,
            instance: INSTANCES.fetch_add(1, Ordering::SeqCst),
            modified: true,
        }
    }

    #[setter]
    fn set_density(&mut self, value: f64) {
        if value != self.density {
            self.modified = true;
            self.density = value;
        }
    }

    #[setter]
    fn set_geodesic(&mut self, value: Geodesic) {
        if value != self.geodesic {
            self.modified = true;
            self.geodesic = value;
        }
    }

    #[setter]
    fn set_material(&mut self, value: String) {
        if value != self.material {
            self.modified = true;
            self.material = value;
        }
    }

    #[setter]
    fn set_ocean(&mut self, value: bool) {
        if value != self.ocean {
            self.modified = true;
            self.ocean = value;
        }
    }

    #[setter]
    fn set_topography(&mut self, value: Option<Topography>) {
        if value != self.topography {
            self.modified = true;
            self.topography = value;
        }
    }
}

impl Geometry {
    pub fn apply(&mut self) -> PyResult<()> {
        let previous = CURRENT.swap(self.instance, Ordering::SeqCst);
        if (previous != self.instance) || self.modified {
            let geodesic: &str = self.geodesic.into();
            let geodesic = CString::new(geodesic)?;
            let topography = self.topography.as_ref().map(|topography| match topography {
                Topography::Dem(path) => CString::new(path.as_str()),
                Topography::Flat(z) => CString::new(format!("flat://{:.6}", z)),
            })
                .transpose()?;
            let topography = match topography {
                None => null(),
                Some(cstr) => cstr.as_ptr(),
            };
            let mut ocean = if self.ocean { 1 } else { 0 };
            to_result(
                unsafe {
                    danton::earth_model(
                        geodesic.as_ptr(),
                        topography,
                        self.density,
                        &mut ocean,
                    )
                },
                None
            )?;
            self.modified = false;
        }
        Ok(())
    }
}


// ===============================================================================================
//
// Tracer interface.
//
// ===============================================================================================

pub struct Tracer<'a> {
    geometry: &'a Geometry,
    tracer: *mut danton::Tracer,
}

impl<'a> Drop for Tracer<'a> {
    fn drop(&mut self) {
        unsafe {
            danton::destroy(&mut (self.tracer as *mut c_void));
        }
    }
}

impl<'a> Tracer<'a> {
    pub fn new(geometry: &'a mut Geometry) -> PyResult<Self> {
        geometry.apply()?;
        let tracer = unsafe { danton::tracer_create() };
        let tracer = Self { geometry, tracer };
        Ok(tracer)
    }

    pub fn medium(&self, position: &GeodeticCoordinates) -> Medium {
        let r = position.to_ecef(self.geometry.geodesic);
        let medium = unsafe { danton::tracer_medium(self.tracer, &r as *const f64) };
        (medium, self.geometry.ocean).into()
    }

    pub fn trace(
        &self,
        position: &GeodeticCoordinates,
        direction: &HorizontalCoordinates
    ) -> (Medium, f64, Medium) {
        let r = position.to_ecef(self.geometry.geodesic);
        let u = direction.to_ecef(self.geometry.geodesic, &position);
        let mut distance: f64 = 0.0;
        let mut next_medium: c_int = -1;
        let medium = unsafe {
            danton::tracer_trace(
                self.tracer,
                &r as *const f64,
                &u as *const f64,
                &mut distance,
                &mut next_medium,
            )
        };
        let medium = (medium, self.geometry.ocean).into();
        let next_medium = (next_medium, self.geometry.ocean).into();
        (medium, distance, next_medium)
    }
}
