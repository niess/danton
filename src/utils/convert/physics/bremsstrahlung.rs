use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::CString;


#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum Bremsstrahlung {
    ABB94,
    KKP95,
    #[default]
    SSR19,
}

impl Bremsstrahlung {
    pub fn as_pumas(&self) -> &str {
        let value = self.to_str();
        &value[0..value.len()-2]
    }
}

impl Convert for Bremsstrahlung {
    #[inline]
    fn what() -> &'static str {
        "bremsstrahlung model"
    }
}

impl<'py> FromPyObject<'py> for Bremsstrahlung {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Bremsstrahlung {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Bremsstrahlung> for &'static str {
    fn from(value: Bremsstrahlung) -> Self {
        value.to_str()
    }
}

impl From<Bremsstrahlung> for CString {
    fn from(value: Bremsstrahlung) -> Self {
        CString::new(value.as_pumas()).unwrap()
    }
}
