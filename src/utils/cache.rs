use crate::utils::error::Error;
use crate::utils::error::ErrorKind::ValueError;
use pyo3::prelude::*;
use ::std::env;
use ::std::path::{Path, PathBuf};


pub fn default_path() -> PyResult<PathBuf> {
    let home = env::var("HOME")
        .map_err(|_| Error::new(ValueError)
            .what("cache")
            .why("could not resolve $HOME")
        )?;
    let cache = Path::new(&home)
        .join(".cache/danton");
    Ok(cache)
}


pub fn path() -> PyResult<PathBuf> {
    let cache = match env::var("DANTON_CACHE") {
        Ok(cache) => Path::new(&cache).to_path_buf(),
        Err(_) => default_path()?,
    };
    Ok(cache)
}
