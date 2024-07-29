use crate::utils::convert::{Convert, Dis, DisModel};
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{NotImplementedError, TypeError};
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::{CString, OsStr};
use ::std::path::Path;


#[derive(Clone, PartialEq)]
pub enum Pdf {
    Model(PdfModel),
    Path(String),
}

#[derive(Clone, Copy, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum PdfModel {
    CT14nlo,
    HERAPDF15NLO,
    NNPDF31sx,
}

impl Pdf {
    pub fn to_c_string(&self, py: Python) -> PyResult<CString> {
        let c_string = match self {
            Self::Model(model) => {
                let prefix = Path::new(crate::PREFIX.get(py).unwrap());
                let filename = match model {
                    PdfModel::CT14nlo => "CT14nlo_0000.dat",
                    PdfModel::HERAPDF15NLO => "HERAPDF15NLO_EIG_0000.dat",
                    PdfModel::NNPDF31sx => "NNPDF31sx_nlo_as_0118_LHCb_nf_6_0000.dat",
                };
                let path = {
                    let mut path = prefix.join("data/pdf");
                    path.push(filename);
                    path
                };
                CString::new(path.to_string_lossy().as_ref())
            },
            Self::Path(path) => CString::new(path.as_str())
        }?;
        Ok(c_string)
    }

    pub fn into_c_string(py: Python, value: Option<&Pdf>, dis: &Dis) -> PyResult<CString> {
        match value {
            None => match dis {
                Dis::Model(model) => match model {
                    DisModel::BGR18 => Pdf::Model(PdfModel::NNPDF31sx).to_c_string(py),
                    DisModel::CSMS => Pdf::Model(PdfModel::HERAPDF15NLO).to_c_string(py),
                    DisModel::LO => Pdf::Model(PdfModel::CT14nlo).to_c_string(py),
                },
                Dis::Path(_) => Pdf::Model(PdfModel::CT14nlo).to_c_string(py),
            },
            Some(pdf) => pdf.to_c_string(py),
        }
    }
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
