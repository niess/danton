use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::CString;


#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum PairProduction {
    KKP68,
    #[default]
    SSR19,
}

impl PairProduction {
    pub fn as_pumas(&self) -> &str {
        let value = self.to_str();
        &value[0..value.len()-2]
    }
}

impl Convert for PairProduction {
    #[inline]
    fn what() -> &'static str {
        "pair-production model"
    }
}

impl<'py> FromPyObject<'py> for PairProduction {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for PairProduction {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<PairProduction> for &'static str {
    fn from(value: PairProduction) -> Self {
        value.to_str()
    }
}

impl From<PairProduction> for CString {
    fn from(value: PairProduction) -> Self {
        CString::new(value.as_pumas()).unwrap()
    }
}
