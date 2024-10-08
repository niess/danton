use crate::utils::convert::Convert;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{NotImplementedError, TypeError};
use enum_variants_strings::EnumVariantsStrings;
use pyo3::prelude::*;
use ::std::ffi::{CString, OsStr};
use ::std::path::Path;


#[derive(Clone, PartialEq)]
pub enum Dis {
    Model(DisModel),
    Path(String),
}

#[derive(Clone, Copy, Default, EnumVariantsStrings, PartialEq)]
#[enum_variants_strings_transform(transform="none")]
pub enum DisModel {
    BGR18,
    #[default]
    CMS11,
    LO,
}

impl Dis {
    pub fn to_c_string(&self, py: Python) -> PyResult<Option<CString>> {
        let dis = match self {
            Dis::Model(model) => match model {
                DisModel::BGR18 => {
                    let prefix = Path::new(crate::PREFIX.get(py).unwrap());
                    Some(prefix.join(format!("data/cs/BGR18.txt")))
                },
                DisModel::CMS11 => {
                    let prefix = Path::new(crate::PREFIX.get(py).unwrap());
                    Some(prefix.join(format!("data/cs/CMS11.txt")))
                },
                DisModel::LO => None,
            },
            Dis::Path(path) => Some(Path::new(path).to_path_buf()),
        };
        let dis = dis
            .map(|dis| CString::new(dis.to_string_lossy().as_ref()))
            .transpose()?;
        Ok(dis)
    }
}

impl Default for Dis {
    fn default() -> Self {
        Self::Model(DisModel::default())
    }
}

impl Convert for DisModel {
    #[inline]
    fn what() -> &'static str {
        "DIS model"
    }
}

impl<'py> FromPyObject<'py> for Dis {
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        let value: String = any.extract()
            .map_err(|_| {
                let err = Error::new(TypeError)
                    .what(DisModel::what())
                    .why("expected a 'str'");
                err.to_err()
            })?;
        let path = Path::new(&value);
        let result = match path.extension().and_then(OsStr::to_str) {
            None => {
                let model = DisModel::from_any(any)?;
                Self::Model(model)
            },
            Some("txt") =>  Self::Path(value),
            Some(other) => {
                let why = format!("cannot load DIS cross-sections from '{}' files", other);
                let err = Error::new(NotImplementedError)
                    .what(DisModel::what())
                    .why(&why);
                return Err(err.to_err())
            }
        };
        Ok(result)
    }
}

impl IntoPy<PyObject> for Dis {
    fn into_py(self, py: Python) -> PyObject {
        match self {
            Self::Model(model) => model.into_any(py),
            Self::Path(path) => path.into_py(py),
        }
    }
}

impl<'a> AsRef<str> for Dis {
    fn as_ref(&self) -> &str {
        match self {
            Dis::Model(model) => model.to_str(),
            Dis::Path(path) => path.as_str(),
        }
    }
}
