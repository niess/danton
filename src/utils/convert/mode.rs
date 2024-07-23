use crate::bindings::danton;
use crate::utils::error::{Error, variant_explain};
use crate::utils::error::ErrorKind::ValueError;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::c_uint;


#[derive(Clone, Copy, EnumVariantsStrings)]
#[enum_variants_strings_transform(transform="lower_case")]
pub enum Mode {
    Backward,
    Forward,
    Grammage,
}

impl<'py> FromPyObject<'py> for Mode {
    fn extract_bound(mode: &Bound<'py, PyAny>) -> PyResult<Self> {
        let mode: String = mode.extract()?;
        let mode = Mode::from_str(&mode)
            .map_err(|options| {
                let why = variant_explain(&mode, options);
                Error::new(ValueError).what("mode").why(&why).to_err()
            })?;
        Ok(mode)
    }
}

impl IntoPy<PyObject> for Mode {
    fn into_py(self, py: Python) -> PyObject {
        self.to_str().into_py(py)
    }
}

impl From<c_uint> for Mode {
    fn from(value: c_uint) -> Self {
        match value {
            danton::BACKWARD => Mode::Backward,
            danton::FORWARD => Mode::Forward,
            danton::GRAMMAGE => Mode::Grammage,
            _ => unreachable!(),
        }
    }
}

impl From<Mode> for c_uint {
    fn from(value: Mode) -> Self {
        match value {
            Mode::Backward => danton::BACKWARD,
            Mode::Forward => danton::FORWARD,
            Mode::Grammage => danton::GRAMMAGE,
        }
    }
}
