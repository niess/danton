use console::style;
use crate::bindings::{alouette, danton, ent, pumas};
use crate::simulation::materials::{Materials, MaterialsData};
use crate::utils::cache;
use crate::utils::convert::{Bremsstrahlung, Dis, Mdf, PairProduction, Pdf, Photonuclear};
use crate::utils::error::{ctrlc_catched, Error};
use crate::utils::error::ErrorKind::{CLibraryException, KeyboardInterrupt, ValueError};
use indicatif::{ProgressBar, ProgressStyle};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyString, PyTuple};
use temp_dir::TempDir;
use ::std::borrow::Cow;
use ::std::ffi::{c_char, c_int, CStr, CString, c_uint};
use ::std::fs::File;
use ::std::os::fd::AsRawFd;
use ::std::ptr::{null, null_mut};


#[pyclass(module="danton")]
pub struct Physics {
    /// The Bremsstrahlung model for tau energy losses.
    #[pyo3(get)]
    bremsstrahlung: Bremsstrahlung,
    /// The e+e- pair-production model for tau energy losses.
    #[pyo3(get)]
    pair_production: PairProduction,
    /// The photonuclear model for tau energy losses.
    #[pyo3(get)]
    photonuclear: Photonuclear,
    /// The Deep Inelastic Scattering (DIS) model for neutrino interactions.
    #[pyo3(get)]
    dis: Dis,
    /// The Parton Distribution Functions (PDF) for neutrino interactions.
    #[pyo3(get)]
    pdf: Option<Pdf>,

    physics: danton::Physics,
    material_index: danton::MaterialIndex,
    materials_instance: Option<usize>,
    pub modified: bool,
}

unsafe impl Send for Physics {}

#[pymethods]
impl Physics {
    #[pyo3(signature=(*, **kwargs))]
    #[new]
    pub fn new<'py>(
        py: Python<'py>,
        kwargs: Option<&Bound<'py, PyDict>>,
    ) -> PyResult<Py<Self>> {
        let bremsstrahlung = Bremsstrahlung::default();
        let pair_production = PairProduction::default();
        let photonuclear = Photonuclear::default();
        let dis = Dis::default();
        let pdf = None;
        let physics = danton::Physics { ent: null_mut(), pumas: null_mut() };
        let material_index = danton::MaterialIndex::default();
        let materials_instance = None;
        let modified = false;

        let physics = Self {
            bremsstrahlung, pair_production, photonuclear, dis, pdf, physics,
            material_index, materials_instance, modified,
        };
        let physics = Bound::new(py, physics)?;

        if let Some(kwargs) = kwargs {
            for (key, value) in kwargs.iter() {
                let key: Bound<PyString> = key.extract()?;
                physics.setattr(key, value)?
            }
        }

        Ok(physics.unbind())
    }

    #[setter]
    fn set_bremsstrahlung(&mut self, value: Bremsstrahlung) {
        if value != self.bremsstrahlung {
            self.bremsstrahlung = value;
            self.modified = true;
        }
    }

    #[setter]
    fn set_pair_production(&mut self, value: PairProduction) {
        if value != self.pair_production {
            self.pair_production = value;
            self.modified = true;
        }
    }

    #[setter]
    fn set_photonuclear(&mut self, value: Photonuclear) {
        if value != self.photonuclear {
            self.photonuclear = value;
            self.modified = true;
        }
    }

    #[setter]
    fn set_dis(&mut self, value: Dis) {
        if value != self.dis {
            self.dis = value;
            self.modified = true;
        }
    }

    #[setter]
    fn set_pdf(&mut self, value: Option<Pdf>) {
        if value != self.pdf {
            self.pdf = value;
            self.modified = true;
        }
    }
}

impl Drop for Physics {
    fn drop(&mut self) {
        self.destroy_physics();
    }
}

impl Physics {
    pub fn apply(
        &mut self,
        py: Python,
        materials: &Materials,
        topography_material: &str
    ) -> PyResult<bool> {
        let modified = if self.physics.ent.is_null() || self.modified ||
            (self.materials_instance != Some(materials.instance)) {

            self.destroy_physics();
            self.create_physics(py, &materials.tag)?;
            self.materials_instance = Some(materials.instance);
            self.modified = false;
            true
        } else {
            false
        };

        self.material_index.topography = match topography_material {
            "Air" => self.material_index.air,
            "Rock" => self.material_index.rock,
            "Water" => self.material_index.water,
            _ => unsafe {
                let mut index: c_int = 0;
                let rc = pumas::physics_material_index(
                    self.physics.pumas,
                    CString::new(topography_material)?.as_ptr(),
                    &mut index,
                );
                if rc != pumas::SUCCESS {
                    let why = format!("unknown material '{}'", topography_material);
                    let err = Error::new(ValueError)
                        .what("topography material")
                        .why(&why)
                        .to_err();
                    return Err(err);
                }
                index
            },
        };

        let physics = unsafe { &mut *danton::physics };
        *physics = self.physics;
        let material_index = unsafe { &mut *danton::material_index };
        *material_index = self.material_index;

        // Apply materials.
        unsafe { danton::materials_set() };

        Ok(modified)
    }

    fn check_alouette(rc: c_uint) -> PyResult<()> {
        if rc == ent::SUCCESS {
            Ok(())
        } else {
            let why = format!("could not initialise alouette (rc = {})", rc);
            let err = Error::new(CLibraryException)
                .what("physics")
                .why(&why)
                .to_err();
            return Err(err);
        }
    }

    fn check_ent(rc: c_uint) -> PyResult<()> {
        if rc == ent::SUCCESS {
            Ok(())
        } else {
            let why = format!("could not create ent physics (rc = {})", rc);
            let err = Error::new(CLibraryException)
                .what("physics")
                .why(&why)
                .to_err();
            return Err(err);
        }
    }

    fn check_pumas(rc: c_uint) -> PyResult<()> {
        if rc == pumas::SUCCESS {
            Ok(())
        } else {
            let why = format!("could not create pumas physics (rc = {})", rc);
            let err = Error::new(CLibraryException)
                .what("physics")
                .why(&why)
                .to_err();
            return Err(err);
        }
    }

    fn create_ent(&self, py: Python) -> PyResult<*mut ent::Physics> {
        let pdf = Pdf::into_c_string(py, self.pdf.as_ref(), &self.dis)?;
        let dis = self.dis.to_c_string(py)?;
        let c_dis = match dis.as_ref() {
            None => null(),
            Some(dis) => dis.as_c_str().as_ptr(),
        };
        let physics = unsafe {
            let mut physics: *mut ent::Physics = null_mut();
            let rc = ent::physics_create(
                &mut physics,
                pdf.as_c_str().as_ptr(),
                c_dis,
            );
            Self::check_ent(rc)?;
            physics
        };
        Ok(physics)
    }

    fn create_physics(&mut self, py: Python, materials: &str) -> PyResult<()> {
        // Load or create Pumas physics.
        let pumas = match self.load_pumas(materials) {
            None => self.create_pumas(py, materials)?,
            Some(pumas) => pumas,
        };
        let mut material_index = danton::MaterialIndex::default();
        unsafe {
            if *danton::tau_mass <= 0.0 {
                Self::check_pumas(
                    pumas::physics_particle(
                        pumas,
                        null_mut(),
                        danton::tau_ctau0,
                        danton::tau_mass
                    )
                )?;
            }
            Self::check_pumas(
                pumas::physics_material_index(
                    pumas,
                    CString::new("Air")?.as_ptr(),
                    &mut material_index.air,
                )
            )?;
            Self::check_pumas(
                pumas::physics_material_index(
                    pumas,
                    CString::new("Rock")?.as_ptr(),
                    &mut material_index.rock,
                )
            )?;
            material_index.topography = material_index.rock;
            Self::check_pumas(
                pumas::physics_material_index(
                    pumas,
                    CString::new("Water")?.as_ptr(),
                    &mut material_index.water,
                )
            )?;
        }

        let ent = self.create_ent(py)?;
        self.physics = danton::Physics { ent, pumas };
        Ok(())
    }

    fn create_pumas(&self, py: Python, materials: &str) -> PyResult<*mut pumas::Physics> {
        let tag = self.pumas_physics_tag();
        let dump_path = cache::path()?
            .join("materials");
        let description = dump_path.join(format!("{}.toml", materials));
        let description = MaterialsData::from_file(py, &description)?;
        std::fs::create_dir_all(&dump_path)?;
        let dump_path = dump_path
            .join(format!("{}-{}.pumas", materials, tag));
        let dedx_path = TempDir::new()?;
        let mdf_path = dedx_path.path().join("materials.xml");
        Mdf::new(py, &description)
            .dump(&mdf_path)?;

        let c_bremsstrahlung = CString::new::<&str>(self.bremsstrahlung.into())?;
        let c_pair_production = CString::new::<&str>(self.pair_production.into())?;
        let c_photonuclear = CString::new::<&str>(self.photonuclear.into())?;
        let mut settings = pumas::PhysicsSettings {
            cutoff: 0.0,
            elastic_ratio: 0.0,
            bremsstrahlung: c_bremsstrahlung.as_c_str().as_ptr(),
            pair_production: c_pair_production.as_c_str().as_ptr(),
            photonuclear: c_photonuclear.as_c_str().as_ptr(),
            n_energies: 0,
            energy: null_mut(),
            update: 0,
            dry: 0,
        };

        let physics = unsafe {
            let mut physics: *mut pumas::Physics = null_mut();
            let mdf_path = CString::new(mdf_path.to_string_lossy().as_ref())?;
            let dedx_path = CString::new(dedx_path.path().to_string_lossy().as_ref())?;
            let mut notifier = Notifier::new(materials);
            let rc = pumas::physics_create(
                &mut physics,
                pumas::TAU,
                mdf_path.as_c_str().as_ptr(),
                dedx_path.as_c_str().as_ptr(),
                &mut settings,
                &mut notifier as *mut Notifier as *mut pumas::PhysicsNotifier,
            );
            if rc == pumas::INTERRUPT {
                // Ctrl-C has been catched.
                let err = Error::new(KeyboardInterrupt)
                    .why("while computing materials tables");
                return Err(err.to_err())
            } else {
                Self::check_pumas(rc)?;
                physics
            }
        };

        // Cache physics data for subsequent usage.
        if let Ok(file) = File::create(dump_path) {
            unsafe {
                let stream = libc::fdopen(
                    file.as_raw_fd(),
                    CStr::from_bytes_with_nul_unchecked(b"wb\0").as_ptr(),
                );
                pumas::physics_dump(physics, stream);
                libc::fclose(stream);
            }
        };

        Ok(physics)
    }

    fn destroy_physics(&mut self) {
        unsafe {
            ent::physics_destroy(&mut self.physics.ent);
            pumas::physics_destroy(&mut self.physics.pumas);
        }
    }

    pub fn initialise() -> PyResult<()> {
        let rc = unsafe { alouette::initialise(null_mut()) };
        Self::check_alouette(rc)
    }

    fn load_pumas(&self, materials: &str) -> Option<*mut pumas::Physics> {
        let tag = self.pumas_physics_tag();
        let path = cache::path().ok()?
            .join(format!("materials/{}-{}.pumas", materials, tag));
        let file = File::open(path).ok()?;
        let mut physics = null_mut();
        let rc = unsafe {
            let stream = libc::fdopen(
                file.as_raw_fd(),
                CStr::from_bytes_with_nul_unchecked(b"rb\0").as_ptr(),
            );
            let rc = pumas::physics_load(&mut physics, stream);
            libc::fclose(stream);
            rc
        };
        if rc != pumas::SUCCESS {
            return None;
        }

        let check_particle = || -> bool {
            let mut particle: c_uint = pumas::MUON;
            let rc = unsafe {
                pumas::physics_particle(physics, &mut particle, null_mut(), null_mut())
            };
            (rc == pumas::SUCCESS) && (particle == pumas::TAU)
        };

        let check_process = |process: c_uint, expected: &str| -> bool {
            unsafe {
                let mut name: *const c_char = null();
                let rc = pumas::physics_dcs(
                    physics,
                    process,
                    &mut name,
                    null_mut(),
                );
                if rc == pumas::SUCCESS {
                    CStr::from_ptr(name)
                        .to_str()
                        .ok()
                        .map(|name| name == expected)
                        .unwrap_or(false)
                } else {
                    false
                }
            }
        };
        if  check_particle() &&
            check_process(pumas::BREMSSTRAHLUNG, self.bremsstrahlung.into()) &&
            check_process(pumas::PAIR_PRODUCTION, self.pair_production.into()) &&
            check_process(pumas::PHOTONUCLEAR, self.photonuclear.into()) {
            Some(physics)
        } else {
            unsafe { pumas::physics_destroy(&mut physics) };
            None
        }
    }

    fn pumas_physics_tag(&self)-> String {
        let bremsstrahlung: &str = self.bremsstrahlung.into();
        let pair_production: &str = self.pair_production.into();
        let photonuclear: &str = self.photonuclear.into();
        format!(
            "{}-{}-{}",
            bremsstrahlung,
            pair_production,
            photonuclear,
        )
    }
}


// ===============================================================================================
//
// Notifier for physics computations.
//
// ===============================================================================================

#[repr(C)]
struct Notifier {
    interface: pumas::PhysicsNotifier,
    bar: Option<ProgressBar>,
    client: String,
    section: usize,
}

impl Notifier {
    const SECTIONS: usize = 4;

    fn new(client: &str) -> Self {
        let interface = pumas::PhysicsNotifier {
            configure: Some(pumas_physics_notifier_configure),
            notify: Some(pumas_physics_notifier_notify),
        };
        Self { interface, bar: None, client: client.to_string(), section: 0 }
    }

    fn configure(&mut self, title: Option<&str>, steps: c_int) {
        self.bar = match self.bar {
            None => {
                self.section += 1;
                let title = title.unwrap();
                let title = if title.starts_with("multiple") {
                    Cow::Borrowed(title)
                } else {
                    Cow::Owned(format!("{}s", title))
                };
                let bar = ProgressBar::new(steps as u64);
                let bar_style = ProgressStyle::with_template(
                    "{msg} [{wide_bar:.dim}] {percent}%, {elapsed})"
                )
                    .unwrap()
                    .progress_chars("=> ");
                bar.set_style(bar_style);
                let section = style(format!("[{}/{}]", self.section, Self::SECTIONS)).dim();
                bar.set_message(format!("({} {} Computing {}", self.client, section, title));
                bar.set_position(0);
                Some(bar)
            },
            Some(_) => None,
        }
    }

    fn notify(&self) {
        self.bar.as_ref().unwrap().inc(1)
    }
}

impl Drop for Notifier {
    fn drop(&mut self) {
        if let Some(bar) = self.bar.as_ref() {
            bar.finish_and_clear()
        }
    }
}

#[no_mangle]
extern "C" fn pumas_physics_notifier_configure(
    slf: *mut pumas::PhysicsNotifier,
    title: *const c_char,
    steps: c_int
) -> c_uint {
    if ctrlc_catched() {
        pumas::INTERRUPT
    } else {
        let notifier = unsafe { &mut *(slf as *mut Notifier) };
        let title = if title.is_null() {
            None
        } else {
            let title = unsafe { CStr::from_ptr(title) };
            Some(title.to_str().unwrap())
        };
        notifier.configure(title, steps);
        pumas::SUCCESS
    }
}

#[no_mangle]
extern "C" fn pumas_physics_notifier_notify(slf: *mut pumas::PhysicsNotifier) -> c_uint {
    if ctrlc_catched() {
        pumas::INTERRUPT
    } else {
        let notifier = unsafe { &*(slf as *mut Notifier) };
        notifier.notify();
        pumas::SUCCESS
    }
}


// ===============================================================================================
//
// Pre-computation interface.
//
// ===============================================================================================

/// Compute materials tables.
#[pyfunction]
#[pyo3(signature=(*args, **kwargs))]
pub fn compute(
    py: Python,
    args: &Bound<PyTuple>,
    kwargs: Option<&Bound<PyDict>>,
) -> PyResult<()> {
    let mut physics = Physics::new(py, kwargs)?
        .bind(py)
        .borrow_mut();
    if args.is_empty() {
        let materials = Materials::new(py, None)?;
        physics
            .create_physics(py, materials.tag.as_str())?;
    } else {
        for arg in args.iter() {
            let arg: String = arg.extract()?;
            let materials = Materials::new(py, Some(arg.as_str()))?;
            physics
                .create_physics(py, materials.tag.as_str())?;
        }
    }
    Ok(())
}
