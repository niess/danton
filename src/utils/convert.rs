use crate::utils::error::{Error, variant_explain};
use crate::utils::error::ErrorKind::ValueError;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;

mod array;
mod geodesy;
mod materials;
mod mdf;
mod medium;
mod mode;
mod physics;
mod state;
mod toml;

pub use array::Array;
pub use geodesy::{Ellipsoid, Geoid, Reference};
pub use mdf::Mdf;
pub use medium::Medium;
pub use mode::Mode;
pub use physics::{Bremsstrahlung, Dis, DisModel, PairProduction, Pdf, Photonuclear};
pub use toml::ToToml;


trait Convert {
    fn what() -> &'static str;

    #[inline]
    fn from_any<'py>(any: &Bound<'py, PyAny>) -> PyResult<Self>
    where
        Self: EnumVariantsStrings,
    {
        let name: String = any.extract()?;
        let value = Self::from_str(&name)
            .map_err(|options| {
                let why = variant_explain(&name, options);
                Error::new(ValueError).what(Self::what()).why(&why).to_err()
            })?;
        Ok(value)
    }

    #[inline]
    fn into_any(self, py: Python) -> PyObject
    where
        Self: EnumVariantsStrings,
    {
        self.to_str().into_py(py)
    }
}
