use crate::utils::error::{Error, variant_explain};
use crate::utils::error::ErrorKind::ValueError;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


#[derive(Clone, Copy, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="upper_case")]
pub enum Geodesic {
    Egm96,
    Prem,
    Wgs84,
}

impl<'py> FromPyObject<'py> for Geodesic {
    fn extract_bound(geodesic: &Bound<'py, PyAny>) -> PyResult<Self> {
        let geodesic: String = geodesic.extract()?;
        let geodesic = Geodesic::from_str(&geodesic)
            .map_err(|options| {
                let why = variant_explain(&geodesic, options);
                Error::new(ValueError).what("geodesic").why(&why).to_err()
            })?;
        Ok(geodesic)
    }
}

impl IntoPy<PyObject> for Geodesic {
    fn into_py(self, py: Python) -> PyObject {
        self.to_str().into_py(py)
    }
}

impl From<Geodesic> for &'static str {
    fn from(value: Geodesic) -> Self {
        value.to_str()
    }
}
