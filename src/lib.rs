use process_path::get_dylib_path;
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use pyo3::exceptions::PySystemError;
use ::std::ffi::CString;
use ::std::path::Path;

mod bindings;
mod simulation;
mod utils;

// XXX Document the geometry and the physics.
// XXX Upgrade PRNG.


static PREFIX: GILOnceCell<String> = GILOnceCell::new();

fn set_prefix(py: Python) -> PyResult<()> {
    let filename = match get_dylib_path() {
        Some(path) => path
                        .to_string_lossy()
                        .to_string(),
        None => return Err(PySystemError::new_err("could not resolve module path")),
    };
    let prefix = match Path::new(&filename).parent() {
        None => ".",
        Some(path) => path.to_str().unwrap(),
    };
    PREFIX
        .set(py, prefix.to_string()).unwrap();
    Ok(())
}

#[pyfunction]
fn finalise_danton() {
    unsafe { bindings::danton::finalise() }
}

/// DecAyiNg Taus frOm Neutrinos (Danton)
#[pymodule]
fn danton(module: &Bound<PyModule>) -> PyResult<()> {
    let py = module.py();

    // Initialise numpy interface.
    utils::numpy::initialise(py)?;

    // Set package prefix.
    set_prefix(py)?;

    // Initialise the C library.
    let prefix = CString::new(PREFIX.get(py).unwrap().as_str()).unwrap();
    unsafe { bindings::danton::initialise(prefix.as_ptr(), None, None); }

    // Register C library finalisation.
    let finalise = wrap_pyfunction!(finalise_danton, module)?;
    py.import_bound("atexit")?
      .call_method1("register", (finalise,))?;

    // Initialise the physics.
    simulation::physics::Physics::initialise()?;

    // Register class object(s).
    module.add_class::<simulation::Simulation>()?;
    module.add_class::<simulation::geometry::Geometry>()?;
    module.add_class::<simulation::physics::Physics>()?;
    module.add_class::<simulation::random::Random>()?;

    // Register function(s).
    module.add_function(wrap_pyfunction!(simulation::particles::particles, module)?)?;

    Ok(())
}
