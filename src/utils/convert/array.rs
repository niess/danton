use crate::utils::numpy::{Dtype, PyArray};
use pyo3::prelude::*;


pub enum Array<'py, T> {
    Copy(Bound<'py, PyAny>),
    Ref(&'py PyArray<T>)
}

impl<'py, T> FromPyObject<'py> for Array<'py, T>
where
    T: Dtype,
{
    fn extract_bound(any: &Bound<'py, PyAny>) -> PyResult<Self> {
        let array: PyResult<&PyArray<T>> = any.extract();
        let array = match array {
            Ok(array) => Self::Ref(array),
            Err(_) => {
                let py = any.py();
                let numpy = PyModule::import_bound(py, "numpy")?;
                let array = numpy
                    .getattr("array")?
                    .call1((any, T::dtype(py)?))?;
                Self::Copy(array)
            },
        };
        Ok(array)
    }
}

impl<'py, T> Array<'py, T>
where
    T: Dtype,
{
    pub fn resolve(&self) -> &'py PyArray<T> {
        match self {
            Self::Copy(any) => any.extract().unwrap(),
            Self::Ref(array) => array,
        }
    }
}
