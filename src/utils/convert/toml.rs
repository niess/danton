use crate::simulation::materials::{Component, Material, MaterialsData};


// ===============================================================================================
//
// Toml writer, for materials.
//
// ===============================================================================================

pub trait ToToml {
    fn to_toml(&self) -> String;
}

impl ToToml for MaterialsData {
    fn to_toml(&self) -> String {
        let mut keys: Vec<_> = self.0.keys().collect();
        keys.sort();
        let mut lines = Vec::<String>::new();
        let n = keys.len();
        for (i, key) in keys.iter().enumerate() {
            lines.push(format!("[{}]", key));
            lines.push(self.0[key.as_str()].to_toml());
            if i < n - 1 {
                lines.push("".to_string());
            }
        }
        lines.join("\n")
    }
}

impl ToToml for Material {
    fn to_toml(&self) -> String {
        let mut lines = Vec::<String>::new();
        lines.push(format!("density = {}", self.density));
        if let Some(mee) = self.I {
            lines.push(format!("I = {}", mee));
        }
        let components: Vec<_> = self.composition.iter()
            .map(|component| component.to_toml())
            .collect();
        let composition = components.join(", ");
        lines.push(format!("composition = {{ {} }}", composition));
        lines.join("\n")
    }
}

impl<'a> ToToml for Component {
    fn to_toml(&self) -> String {
        format!("{} = {}", self.name, self.weight)
    }
}
