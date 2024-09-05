use crate::simulation::materials::{AtomicElement, Component, ELEMENTS, Material, MaterialsData};
use pyo3::prelude::*;
use ::std::path::Path;


// ===============================================================================================
//
// Materials Description File (MDF) for Pumas.
//
// ===============================================================================================

pub struct Mdf (String);

impl Mdf {
    pub fn new(py: Python, materials: &MaterialsData) -> Self {
        let mut lines = Vec::<String>::new();
        lines.push("<pumas>".to_string());

        let table = ELEMENTS.get(py).unwrap();
        let mut elements = Vec::<&str>::new();
        for material in materials.0.values() {
            for Component { name, .. } in material.composition.iter() {
                elements.push(name)
            }
        }
        elements.sort();
        elements.dedup();

        for key in elements {
            let element = table.0.get(key).unwrap();
            let element = element.to_xml(key);
            lines.push(element);
        }

        let mut keys: Vec<_> = materials.0.keys().collect();
        keys.sort();
        for key in keys.drain(..) {
            let material = &materials.0[key];
            let material = material.to_xml(key);
            lines.push(material);
        }

        lines.push("</pumas>".to_string());
        let mdf = lines.join("\n");
        Self (mdf)
    }

    pub fn dump<P: AsRef<Path>>(&self, destination: P) -> PyResult<()> {
        std::fs::write(destination, self.0.as_str())?;
        Ok(())
    }
}


// ===============================================================================================
//
// Xml conversions.
//
// ===============================================================================================

trait ToXml {
    fn to_xml(&self, key: &str) -> String;
}

impl ToXml for AtomicElement {
    fn to_xml(&self, key: &str) -> String {
        format!(
            "<element name=\"{}\" Z=\"{}\" A=\"{}\" I=\"{}\" />",
            key,
            self.Z,
            self.A, // g/mol
            self.I, // eV.
        )
    }
}

impl ToXml for Material {
    fn to_xml(&self, key: &str) -> String {
        let mut lines = Vec::<String>::new();
        let header = match self.I {
            Some(mee) => format!(
                "<material name=\"{}\" density=\"{}\" I=\"{}\">",
                key,
                self.density * 1E-03, // g/cm3.
                mee * 1E+09, // eV
            ),
            None => format!(
                "<material name=\"{}\" density=\"{}\">",
                key,
                self.density * 1E-03, // g/cm3.
            ),
        };
        lines.push(header);
        let mut composition: Vec<_> = self.composition.iter().collect();
        composition.sort_by(|a, b| a.name.cmp(&b.name));
        for Component { name, weight } in composition.drain(..) {
            let line = format!(
                "    <component name=\"{}\" fraction=\"{}\" />",
                name,
                weight,
            );
            lines.push(line);
        }
        lines.push("</material>".to_string());
        lines.join("\n")
    }
}
