use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum PairProduction {
    KKP,
    #[default]
    SSR,
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
