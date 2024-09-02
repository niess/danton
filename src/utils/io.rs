use pyo3::prelude::*;
use pyo3::exceptions::PyFileNotFoundError;
use pyo3::types::PyDict;
use std::path::Path;


// ===============================================================================================
//
// Dict loader(s).
//
// ===============================================================================================

pub trait ConfigFormat {
    fn import_module<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyModule>>;

    fn load_dict<'py>(py: Python<'py>, path: &Path) -> PyResult<Bound<'py, PyDict>> {
        let content = std::fs::read_to_string(path)
            .map_err(|err| match err.kind() {
                std::io::ErrorKind::NotFound => {
                    let path = format!("No such file or directory '{}'", path.display());
                    PyFileNotFoundError::new_err(path)
                },
                _ => err.into(),
            })?;
        let module = Self::import_module(py)?;
        let loads = module.getattr("loads")?;
        let content = loads.call1((content,))?;
        let dict: Bound<PyDict> = content.extract()?;
        Ok(dict)
    }
}

pub struct Toml;

impl ConfigFormat for Toml {
    fn import_module<'py>(py: Python<'py>) -> PyResult<Bound<'py, PyModule>> {
        py.import_bound("tomllib")
            .or_else(|_| py.import_bound("tomli"))
    }
}
