use crate::utils::convert::Convert;
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;


// ===============================================================================================
//
// Geoid models.
//
// ===============================================================================================

#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="upper_case")]
pub enum Geoid {
    Egm96,
    #[default]
    Prem,
    Wgs84,
}

impl Convert for Geoid {
    #[inline]
    fn what() -> &'static str {
        "geoid"
    }
}

impl<'py> FromPyObject<'py> for Geoid {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Geoid {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Geoid> for &'static str {
    fn from(value: Geoid) -> Self {
        value.to_str()
    }
}


// ===============================================================================================
//
// Ellipsoid models.
//
// ===============================================================================================

#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="upper_case")]
pub enum Ellipsoid {
    #[default]
    Prem,
    Wgs84,
}

impl Convert for Ellipsoid {
    #[inline]
    fn what() -> &'static str {
        "geoid"
    }
}

impl<'py> FromPyObject<'py> for Ellipsoid {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Ellipsoid {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Ellipsoid> for &'static str {
    fn from(value: Ellipsoid) -> Self {
        value.to_str()
    }
}

impl From<Geoid> for Ellipsoid {
    fn from(value: Geoid) -> Self {
        match value {
            Geoid::Prem => Ellipsoid::Prem,
            Geoid::Egm96 | Geoid::Wgs84 => Ellipsoid::Wgs84,
        }
    }
}

impl std::fmt::Display for Ellipsoid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.to_str())
    }
}


// ===============================================================================================
//
// Elevation reference.
//
// ===============================================================================================

#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="lower_case")]
pub enum Reference {
    #[default]
    Ellipsoid,
    Geoid,
}

impl Convert for Reference {
    #[inline]
    fn what() -> &'static str {
        "reference"
    }
}

impl<'py> FromPyObject<'py> for Reference {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        Self::from_any(any)
    }
}

impl IntoPy<PyObject> for Reference {
    fn into_py(self, py: Python) -> PyObject {
        self.into_any(py)
    }
}

impl From<Reference> for &'static str {
    fn from(value: Reference) -> Self {
        value.to_str()
    }
}
