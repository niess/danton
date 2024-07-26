use crate::bindings::danton;
use crate::utils::convert::Convert;
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

impl Convert for Mode {
    #[inline]
    fn what() -> &'static str {
        "mode"
    }
}

impl<'py> FromPyObject<'py> for Mode {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Mode {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Mode> for &'static str {
    fn from(value: Mode) -> Self {
        value.to_str()
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
