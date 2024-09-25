use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};
use pyo3::sync::GILOnceCell;


pub struct NamedTuple<const N: usize> {
    name: &'static str,
    fields: [&'static str; N],
    object: GILOnceCell<PyObject>,
}

impl<const N: usize> NamedTuple<N> {
    pub const fn new(name: &'static str, fields: [&'static str; N]) -> Self {
        let object = GILOnceCell::new();
        Self { name, fields, object }
    }

    pub fn instance<'py>(
        &self,
        py: Python<'py>,
        args: impl IntoPy<Py<PyTuple>>,
    ) -> PyResult<Bound<'py, PyAny>> {
        let kwargs = PyTuple::new_bound(py, [("module", "danton")]);
        let kwargs = PyDict::from_sequence_bound(kwargs.as_any())?;
        let tuple = self.object.get_or_try_init(py, || {
            let tag = self.fields.join("_");
            let name = format!("_{}_{}", self.name, tag);
            let tp = py.import_bound("collections")
                .and_then(|m| m.getattr("namedtuple"))
                .and_then(|m| m.call((name.as_str(), self.fields), Some(&kwargs)))?;
            py.import_bound("danton")
                .and_then(|m| m.setattr(name.as_str(), tp.clone()))?;
            Ok::<_, PyErr>(tp.unbind())
        })?.bind(py);
        let tuple = tuple.call1(args)?;
        Ok(tuple)
    }
}
