use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum Bremsstrahlung {
    ABB,
    KKP,
    #[default]
    SSR,
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
