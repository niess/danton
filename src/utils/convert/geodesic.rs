use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


#[derive(Clone, Copy, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="upper_case")]
pub enum Geodesic {
    Egm96,
    Prem,
    Wgs84,
}

impl Convert for Geodesic {
    #[inline]
    fn what() -> &'static str {
        "geodesic"
    }
}

impl<'py> FromPyObject<'py> for Geodesic {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Geodesic {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Geodesic> for &'static str {
    fn from(value: Geodesic) -> Self {
        value.to_str()
    }
}
