use crate::simulation::materials::{Component, Material, MaterialsData};
use crate::utils::error::{Error, ErrorKind};
use crate::utils::error::ErrorKind::{KeyError, TypeError, ValueError};
use crate::utils::namespace::Namespace;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyTuple};


// ===============================================================================================
//
// Materials data extraction.
//
// ===============================================================================================

impl<'py> TryFrom<Bound<'py, PyDict>> for MaterialsData {
    type Error = PyErr;

    fn try_from(value: Bound<'py, PyDict>) -> PyResult<Self> {
        let to_err = |expected: &str, found: &Bound<PyAny>| {
            let tp = found.get_type();
            let why = format!("expected a '{}', found a '{:?}'", expected, tp);
            Error::new(TypeError)
                .what("materials")
                .why(&why)
                .to_err()
        };

        let mut materials = Self::new();
        for (k, v) in value.iter() {
            let k: String = k.extract()
                .map_err(|_| to_err("string", &k))?;
            let v: Bound<PyDict> = v.extract()
                .map_err(|_| to_err("dict", &v))?;
            let material: Material = (k.as_str(), &v, &materials).try_into()?;
            materials.0.insert(k, material);
        }

        Ok(materials)
    }
}

// ===============================================================================================
//
// Material extraction.
//
// ===============================================================================================

type Context<'a, 'py> = (&'a str, &'a Bound<'py, PyDict>, &'a MaterialsData);

impl<'a, 'py> TryFrom<Context<'a, 'py>> for Material {
    type Error = PyErr;

    fn try_from(value: Context) -> PyResult<Self> {
        let (name, data, others) = value;
        let py = data.py();
        let to_err = |kind: ErrorKind, why: &str| -> PyErr {
            let what = format!("'{}' material", name);
            Error::new(kind)
                .what(&what)
                .why(why)
                .to_err()
        };

        let mut density: Option<f64> = None;
        let mut mee: Option<f64> = None;
        let mut composition: Option<Bound<PyAny>> = None;
        for (k, v) in data.iter() {
            let k: String = k.extract()
                .map_err(|_| to_err(TypeError, "key is not a string"))?;
            match k.as_str() {
                "density" => {
                    let v: f64 = v.extract()
                        .map_err(|_| to_err(ValueError, "'density' is not a float"))?;
                    density = Some(v);
                },
                "composition" => composition = Some(v),
                "I" => {
                    let v: f64 = v.extract()
                        .map_err(|_| to_err(ValueError, "'I' is not a float"))?;
                    mee = Some(v);
                },
                _ => {
                    return Err(to_err(KeyError, &format!("invalid property '{}'", k)));
                },
            }
        }
        let density = density
            .ok_or_else(|| to_err(KeyError, "missing 'density'"))?;
        let composition = composition
            .ok_or_else(|| to_err(KeyError, "missing 'composition'"))?;

        let formula: Option<String> = composition.extract().ok();
        let material = match formula {
            Some(formula) => Material::from_formula(py, density, formula.as_str(), mee)
                .map_err(|(kind, why)| to_err(kind, &why))?,
            None => {
                let composition: Bound<PyDict> = composition.extract()
                    .map_err(|_| {
                        let tp = composition.get_type();
                        let why = format!(
                            "expected a 'dict' or a 'string' for 'composition', found a '{:?}'",
                            tp,
                        );
                        to_err(TypeError, &why)
                    })?;
                let mut components = Vec::<Component>::new();
                for (k, v) in composition {
                    let name: String = k.extract()
                        .map_err(|_| to_err(TypeError, "key is not a string"))?;
                    let weight: f64 = v.extract()
                        .map_err(|_| {
                            let why = format!("weight for '{}' is not a float", k);
                            to_err(TypeError, &why)
                        })?;
                    if weight > 0.0 {
                        components.push(Component::new(name, weight))
                    }
                }
                Material::from_composition(py, density, &components, others, mee)
                    .map_err(|(kind, why)| to_err(kind, &why))?
            },
        };
        Ok(material)
    }
}


// ===============================================================================================
//
// Conversion to a Namespace.
//
// ===============================================================================================

impl ToPyObject for Material {
    fn to_object(&self, py: Python) -> PyObject {
        let composition = PyTuple::new_bound(py, &self.composition);
        match self.I {
            Some(mee) => {
                Namespace::new(py, &[
                    ("density", self.density.to_object(py)),
                    ("I", mee.to_object(py)),
                    ("composition", composition.into_any().unbind()),
                ]).unwrap().unbind()
            },
            None => {
                Namespace::new(py, &[
                    ("density", self.density.to_object(py)),
                    ("composition", composition.into_any().unbind()),
                ]).unwrap().unbind()
            },
        }
    }
}

impl ToPyObject for Component {
    fn to_object(&self, py: Python) -> PyObject {
        Namespace::new(py, &[
            ("name", self.name.to_object(py)),
            ("weight", self.weight.to_object(py)),
        ]).unwrap().unbind()
    }
}
