use crate::utils::convert::Convert;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{NotImplementedError, TypeError};
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::OsStr;
use ::std::path::Path;


#[derive(Clone)]
pub enum Pdf {
    Model(PdfModel),
    Path(String),
}

#[derive(Clone, Copy, EnumVariantsStrings)]
#[enum_variants_strings_transform(transform="none")]
pub enum PdfModel {
    CT14nlo,
    HERAPDF15NLO,
    NNPDF31sx,
}

impl Convert for PdfModel {
    #[inline]
    fn what() -> &'static str {
        "PDF model"
    }
}

impl<'py> FromPyObject<'py> for Pdf {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        let value: String = any.extract()
            .map_err(|_| {
                let err = Error::new(TypeError)
                    .what(PdfModel::what())
                    .why("expected a 'str'");
                err.to_err()
            })?;
        let path = Path::new(&value);
        let result = match path.extension().and_then(OsStr::to_str) {
            None => {
                let model = PdfModel::from_any(any)?;
                Self::Model(model)
            },
            Some("txt") =>  Self::Path(value),
            Some(other) => {
                let why = format!("cannot load PDF from a '{}' file", other);
                let err = Error::new(NotImplementedError)
                    .what(PdfModel::what())
                    .why(&why);
                return Err(err.to_err())
            }
        };
        Ok(result)
    }
}

impl IntoPy<PyObject> for Pdf {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::Model(model) => model.into_any(py),
            Self::Path(path) => path.into_py(py),
        }
    }
}

impl<'a> AsRef<str> for Pdf {
    fn as_ref(&self) -> &str {
        match self {
            Pdf::Model(model) => model.to_str(),
            Pdf::Path(path) => path.as_str(),
        }
    }
}
