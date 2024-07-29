use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum Photonuclear {
    BBKS,
    BM,
    #[default]
    DRSS,
}

impl Convert for Photonuclear {
    #[inline]
    fn what() -> &'static str {
        "photonuclear model"
    }
}

impl<'py> FromPyObject<'py> for Photonuclear {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Photonuclear {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Photonuclear> for &'static str {
    fn from(value: Photonuclear) -> Self {
        value.to_str()
    }
}
