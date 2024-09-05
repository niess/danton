use crate::utils::cache;
use crate::utils::convert::ToToml;
use crate::utils::error::Error;
use crate::utils::error::ErrorKind::{self, KeyError, ValueError};
use crate::utils::io::{ConfigFormat, Toml};
use pyo3::prelude::*;
use pyo3::sync::GILOnceCell;
use regex::Regex;
use ::std::borrow::Cow;
use ::std::collections::HashMap;
use ::std::ffi::OsStr;
use ::std::path::Path;
use ::std::sync::atomic::{AtomicUsize, Ordering};


// ===============================================================================================
//
// Materials interface.
//
// ===============================================================================================

pub const DEFAULT_MATERIALS: &'static str = "default";

pub fn initialise(py: Python) -> PyResult<()> {
    // Initialise atomic elements data.
    ElementsTable::initialise(py);

    // Synchronize default materials.
    let path = Path::new(crate::PREFIX.get(py).unwrap())
        .join(format!("data/materials/{}.toml", DEFAULT_MATERIALS));
    let data = MaterialsData::from_file(py, path)?;
    let _unused = data.sync(py, DEFAULT_MATERIALS)?;

    Ok(())
}

static INSTANCES: AtomicUsize = AtomicUsize::new(0);

#[pyclass(frozen, module="danton")]
pub struct Materials {
    #[pyo3(get)]
    /// A name identifying the materials set.
    pub tag: String,

    pub data: MaterialsData,
    pub instance: usize,
}

#[pymethods]
impl Materials {
    #[new]
    pub fn new(py: Python, arg: Option<&str>) -> PyResult<Self> {
        let arg = arg.unwrap_or(DEFAULT_MATERIALS);
        let path = Path::new(arg);
        let get_tag = || -> PyResult<&str> {
            path
                .file_stem()
                .and_then(OsStr::to_str)
                .ok_or_else(|| {
                    let stem = path.file_stem()
                        .and_then(|stem| Some(stem.to_string_lossy()))
                        .unwrap_or(path.to_string_lossy());
                    let why = format!("invalid tag '{}'", stem);
                    Error::new(ValueError)
                        .what("materials")
                        .why(&why)
                        .to_err()
                })
        };

        let (tag, data) = match path.extension().and_then(OsStr::to_str) {
            None => match path.parent().filter(|parent| *parent != Path::new("")) {
                // Value should be a valid materials tag.
                Some(_) => {
                    let why = format!("invalid tag '{}'", arg);
                    let err = Error::new(ValueError)
                        .what("materials")
                        .why(&why);
                    return Err(err.to_err())
                },
                None => {
                    let tag = get_tag()?;
                    let cached = cache::path()?
                        .join(format!("materials/{}.toml", tag));
                    if cached.try_exists().unwrap_or(false) {
                        let data = MaterialsData::from_file(py, cached)?;
                        (tag.to_string(), data)
                    } else {
                        let why = format!("unknown tag '{}'", tag);
                        let err = Error::new(ValueError)
                            .what("materials")
                            .why(&why);
                        return Err(err.to_err())
                    }
                },
            },
            Some("toml") => { // Value should be a materials definition file.
                let tag = get_tag()?;
                if tag == DEFAULT_MATERIALS {
                    let why = format!("cannot modify '{}'", DEFAULT_MATERIALS);
                    let err = Error::new(ValueError)
                        .what("materials")
                        .why(&why);
                    return Err(err.to_err())
                }

                let data = {
                    let mut data = MaterialsData::from_file(py, path)?;
                    let mut default_data: Option<MaterialsData> = None;
                    for material in ["Air", "Rock", "Water"] {
                        data.0.entry(material.to_string()).or_insert_with(|| {
                            let default_data = default_data.get_or_insert_with(|| {
                                let path = Path::new(crate::PREFIX.get(py).unwrap())
                                    .join(format!(
                                        "data/materials/{}.toml",
                                        DEFAULT_MATERIALS
                                    ));
                                MaterialsData::from_file(py, path).unwrap()
                            });
                            default_data.0[material].clone()
                        });
                    }
                    data
                };
                data.sync(py, tag)?;
                (tag.to_string(), data)
            },
            _ => {
                let why = format!("invalid file format '{}'", path.display());
                let err = Error::new(ValueError)
                    .what("materials")
                    .why(&why);
                return Err(err.to_err())
            },
        };

        let instance = INSTANCES.fetch_add(1, Ordering::SeqCst);
        let materials = Self { tag, data, instance };
        Ok(materials)
    }

    fn __getitem__(&self, py: Python, name: &str) -> PyResult<PyObject> {
        let material = self.data.0.get(name)
            .ok_or_else(|| {
                let why = format!("unknown material '{}'", name);
                Error::new(KeyError)
                    .what("material")
                    .why(&why)
                    .to_err()
            })?;
        Ok(material.to_object(py))
    }
}

// ===============================================================================================
//
// Atomic elements.
//
// ===============================================================================================

#[allow(non_snake_case)]
pub struct AtomicElement {
    pub Z: u32,
    pub A: f64,
    pub I: f64,
}

#[repr(transparent)]
pub struct ElementsTable (pub HashMap<&'static str, AtomicElement>);

pub static ELEMENTS: GILOnceCell<ElementsTable> = GILOnceCell::new();

impl ElementsTable {
    fn initialise(py: Python) {
        // Data from https://pdg.lbl.gov/2024/AtomicNuclearProperties/index.html.
        let table = HashMap::from([
            ("H",  AtomicElement { Z: 1,   A: 1.008,   I: 19.2   }),
            ("D",  AtomicElement { Z: 1,   A: 2.0141,  I: 19.2   }),
            ("He", AtomicElement { Z: 2,   A: 4.0026,  I: 41.8   }),
            ("Li", AtomicElement { Z: 3,   A: 6.94,    I: 40.0   }),
            ("Be", AtomicElement { Z: 4,   A: 9.01218, I: 63.7   }),
            ("B",  AtomicElement { Z: 5,   A: 10.81,   I: 76.0   }),
            ("C",  AtomicElement { Z: 6,   A: 12.0107, I: 78.0   }),
            ("N",  AtomicElement { Z: 7,   A: 14.007,  I: 82.0   }),
            ("O",  AtomicElement { Z: 8,   A: 15.999,  I: 95.0   }),
            ("F",  AtomicElement { Z: 9,   A: 18.9984, I: 115.0  }),
            ("Ne", AtomicElement { Z: 10,  A: 20.1797, I: 137.0  }),
            ("Rk", AtomicElement { Z: 11,  A: 22.0,    I: 136.4  }), // Fictitious Rockium.
            ("Na", AtomicElement { Z: 11,  A: 22.9898, I: 149.0  }),
            ("Mg", AtomicElement { Z: 12,  A: 24.305,  I: 156.0  }),
            ("Al", AtomicElement { Z: 13,  A: 26.9815, I: 166.0  }),
            ("Si", AtomicElement { Z: 14,  A: 28.0855, I: 173.0  }),
            ("P",  AtomicElement { Z: 15,  A: 30.9738, I: 173.0  }),
            ("S",  AtomicElement { Z: 16,  A: 32.065,  I: 180.0  }),
            ("Cl", AtomicElement { Z: 17,  A: 35.453,  I: 174.0  }),
            ("Ar", AtomicElement { Z: 18,  A: 39.948,  I: 188.0  }),
            ("K",  AtomicElement { Z: 19,  A: 39.0983, I: 190.0  }),
            ("Ca", AtomicElement { Z: 20,  A: 40.078,  I: 191.0  }),
            ("Sc", AtomicElement { Z: 21,  A: 44.9559, I: 216.0  }),
            ("Ti", AtomicElement { Z: 22,  A: 47.867,  I: 233.0  }),
            ("V",  AtomicElement { Z: 23,  A: 50.9415, I: 245.0  }),
            ("Cr", AtomicElement { Z: 24,  A: 51.9961, I: 257.0  }),
            ("Mn", AtomicElement { Z: 25,  A: 54.938,  I: 272.0  }),
            ("Fe", AtomicElement { Z: 26,  A: 55.845,  I: 286.0  }),
            ("Co", AtomicElement { Z: 27,  A: 58.9332, I: 297.0  }),
            ("Ni", AtomicElement { Z: 28,  A: 58.6934, I: 311.0  }),
            ("Cu", AtomicElement { Z: 29,  A: 63.546,  I: 322.0  }),
            ("Zn", AtomicElement { Z: 30,  A: 65.38,   I: 330.0  }),
            ("Ga", AtomicElement { Z: 31,  A: 69.723,  I: 334.0  }),
            ("Ge", AtomicElement { Z: 32,  A: 72.63,   I: 350.0  }),
            ("As", AtomicElement { Z: 33,  A: 74.9216, I: 347.0  }),
            ("Se", AtomicElement { Z: 34,  A: 78.971,  I: 348.0  }),
            ("Br", AtomicElement { Z: 35,  A: 79.904,  I: 357.0  }),
            ("Kr", AtomicElement { Z: 36,  A: 83.798,  I: 352.0  }),
            ("Rb", AtomicElement { Z: 37,  A: 85.4678, I: 363.0  }),
            ("Sr", AtomicElement { Z: 38,  A: 87.62,   I: 366.0  }),
            ("Y",  AtomicElement { Z: 39,  A: 88.9058, I: 379.0  }),
            ("Zr", AtomicElement { Z: 40,  A: 91.224,  I: 393.0  }),
            ("Nb", AtomicElement { Z: 41,  A: 92.9064, I: 417.0  }),
            ("Mo", AtomicElement { Z: 42,  A: 95.95,   I: 424.0  }),
            ("Tc", AtomicElement { Z: 43,  A: 97.9072, I: 428.0  }),
            ("Ru", AtomicElement { Z: 44,  A: 101.07,  I: 441.0  }),
            ("Rh", AtomicElement { Z: 45,  A: 102.906, I: 449.0  }),
            ("Pd", AtomicElement { Z: 46,  A: 106.42,  I: 470.0  }),
            ("Ag", AtomicElement { Z: 47,  A: 107.868, I: 470.0  }),
            ("Cd", AtomicElement { Z: 48,  A: 112.414, I: 469.0  }),
            ("In", AtomicElement { Z: 49,  A: 114.818, I: 488.0  }),
            ("Sn", AtomicElement { Z: 50,  A: 118.71,  I: 488.0  }),
            ("Sb", AtomicElement { Z: 51,  A: 121.76,  I: 487.0  }),
            ("Te", AtomicElement { Z: 52,  A: 127.6,   I: 485.0  }),
            ("I",  AtomicElement { Z: 53,  A: 126.904, I: 491.0  }),
            ("Xe", AtomicElement { Z: 54,  A: 131.293, I: 482.0  }),
            ("Cs", AtomicElement { Z: 55,  A: 132.905, I: 488.0  }),
            ("Ba", AtomicElement { Z: 56,  A: 137.327, I: 491.0  }),
            ("La", AtomicElement { Z: 57,  A: 138.905, I: 501.0  }),
            ("Ce", AtomicElement { Z: 58,  A: 140.116, I: 523.0  }),
            ("Pr", AtomicElement { Z: 59,  A: 140.908, I: 535.0  }),
            ("Nd", AtomicElement { Z: 60,  A: 144.242, I: 546.0  }),
            ("Pm", AtomicElement { Z: 61,  A: 144.913, I: 560.0  }),
            ("Sm", AtomicElement { Z: 62,  A: 150.36,  I: 574.0  }),
            ("Eu", AtomicElement { Z: 63,  A: 151.964, I: 580.0  }),
            ("Gd", AtomicElement { Z: 64,  A: 157.25,  I: 591.0  }),
            ("Tb", AtomicElement { Z: 65,  A: 158.925, I: 614.0  }),
            ("Dy", AtomicElement { Z: 66,  A: 162.5,   I: 628.0  }),
            ("Ho", AtomicElement { Z: 67,  A: 164.93,  I: 650.0  }),
            ("Er", AtomicElement { Z: 68,  A: 167.259, I: 658.0  }),
            ("Tm", AtomicElement { Z: 69,  A: 168.934, I: 674.0  }),
            ("Yb", AtomicElement { Z: 70,  A: 173.054, I: 684.0  }),
            ("Lu", AtomicElement { Z: 71,  A: 174.967, I: 694.0  }),
            ("Hf", AtomicElement { Z: 72,  A: 178.49,  I: 705.0  }),
            ("Ta", AtomicElement { Z: 73,  A: 180.948, I: 718.0  }),
            ("W",  AtomicElement { Z: 74,  A: 183.84,  I: 727.0  }),
            ("Re", AtomicElement { Z: 75,  A: 186.207, I: 736.0  }),
            ("Os", AtomicElement { Z: 76,  A: 190.23,  I: 746.0  }),
            ("Ir", AtomicElement { Z: 77,  A: 192.217, I: 757.0  }),
            ("Pt", AtomicElement { Z: 78,  A: 195.084, I: 790.0  }),
            ("Au", AtomicElement { Z: 79,  A: 196.967, I: 790.0  }),
            ("Hg", AtomicElement { Z: 80,  A: 200.592, I: 800.0  }),
            ("Tl", AtomicElement { Z: 81,  A: 204.38,  I: 810.0  }),
            ("Pb", AtomicElement { Z: 82,  A: 207.2,   I: 823.0  }),
            ("Bi", AtomicElement { Z: 83,  A: 208.98,  I: 823.0  }),
            ("Po", AtomicElement { Z: 84,  A: 208.982, I: 830.0  }),
            ("At", AtomicElement { Z: 85,  A: 209.987, I: 825.0  }),
            ("Rn", AtomicElement { Z: 86,  A: 222.018, I: 794.0  }),
            ("Fr", AtomicElement { Z: 87,  A: 223.02,  I: 827.0  }),
            ("Ra", AtomicElement { Z: 88,  A: 226.025, I: 826.0  }),
            ("Ac", AtomicElement { Z: 89,  A: 227.028, I: 841.0  }),
            ("Th", AtomicElement { Z: 90,  A: 232.038, I: 847.0  }),
            ("Pa", AtomicElement { Z: 91,  A: 231.036, I: 878.0  }),
            ("U",  AtomicElement { Z: 92,  A: 238.029, I: 890.0  }),
            ("Np", AtomicElement { Z: 93,  A: 237.048, I: 902.0  }),
            ("Pu", AtomicElement { Z: 94,  A: 244.064, I: 921.0  }),
            ("Am", AtomicElement { Z: 95,  A: 243.061, I: 934.0  }),
            ("Cm", AtomicElement { Z: 96,  A: 247.07,  I: 939.0  }),
            ("Bk", AtomicElement { Z: 97,  A: 247.07,  I: 952.0  }),
            ("Cf", AtomicElement { Z: 98,  A: 251.08,  I: 966.0  }),
            ("Es", AtomicElement { Z: 99,  A: 252.083, I: 980.0  }),
            ("Fm", AtomicElement { Z: 100, A: 257.095, I: 994.0  }),
            ("Md", AtomicElement { Z: 101, A: 258.098, I: 1007.0 }),
            ("No", AtomicElement { Z: 102, A: 259.101, I: 1020.0 }),
            ("Lr", AtomicElement { Z: 103, A: 262.11,  I: 1034.0 }),
            ("Rf", AtomicElement { Z: 104, A: 267.122, I: 1047.0 }),
            ("Db", AtomicElement { Z: 105, A: 268.126, I: 1061.0 }),
            ("Sg", AtomicElement { Z: 106, A: 269.129, I: 1074.0 }),
            ("Bh", AtomicElement { Z: 107, A: 270.133, I: 1087.0 }),
            ("Hs", AtomicElement { Z: 108, A: 269.134, I: 1102.0 }),
            ("Mt", AtomicElement { Z: 109, A: 278.156, I: 1115.0 }),
            ("Ds", AtomicElement { Z: 110, A: 281.164, I: 1129.0 }),
            ("Rg", AtomicElement { Z: 111, A: 282.169, I: 1143.0 }),
            ("Cn", AtomicElement { Z: 112, A: 285.177, I: 1156.0 }),
            ("Nh", AtomicElement { Z: 113, A: 286.182, I: 1171.0 }),
            ("Fl", AtomicElement { Z: 114, A: 289.19,  I: 1185.0 }),
            ("Mc", AtomicElement { Z: 115, A: 289.194, I: 1199.0 }),
            ("Lv", AtomicElement { Z: 116, A: 293.204, I: 1213.0 }),
            ("Ts", AtomicElement { Z: 117, A: 294.211, I: 1227.0 }),
            ("Og", AtomicElement { Z: 118, A: 294.214, I: 1242.0 }),
        ]);
        let _unused = ELEMENTS
            .set(py, Self (table));
    }
}


// ===============================================================================================
//
// Material data.
//
// ===============================================================================================

#[allow(non_snake_case)]
#[derive(Clone, Debug)]
pub struct Material {
    pub density: f64,
    mass: f64,
    pub composition: Vec<Component>, // Beware: mass fractions.
    pub I: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct Component {
    pub name: String,
    pub weight: f64,
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        if self.density != other.density {
            return false;
        }
        if (self.mass - other.mass).abs() > 1E-09 { // Prevent rounding errors.
            return false;
        }
        if self.I != other.I {
            return false;
        }
        self.composition.eq(&other.composition)
    }
}

impl PartialEq for Component {
    fn eq(&self, other: &Self) -> bool {
        if self.name.eq(&other.name) {
            (self.weight - other.weight).abs() <= 1E-09 // Prevent rounding errors.
        } else {
            false
        }
    }
}

type ErrorData = (ErrorKind, String);

impl Material {
    #[allow(non_snake_case)]
    pub fn new(density: f64, mass:f64, mut composition: Vec<Component>, I: Option<f64>) -> Self {
        composition.sort_by(|a,b| a.name.cmp(&b.name));
        Self { density, mass, composition, I }
    }

    pub fn from_composition(
        py: Python,
        density: f64,
        composition: &[Component], // Beware: mass fractions.
        others: &MaterialsData,
        #[allow(non_snake_case)]
        I: Option<f64>,
    ) -> Result<Self, ErrorData> {
        let table = ELEMENTS.get(py).unwrap();
        let mut weights = HashMap::<&'static str, f64>::new();
        let mut sum = 0.0;
        for Component { name, weight: wi } in composition.iter() {
            let xi = match table.0.get(name.as_str()) {
                Some(element) => {
                    let xi = wi / element.A;
                    weights
                        .entry(name.as_ref())
                        .and_modify(|x| *x += xi)
                        .or_insert(xi);
                    xi
                },
                None => {
                    let material: Option<Cow<Material>> = others.0.get(name.as_str())
                        .map(|material| Cow::Borrowed(material))
                        .or_else(||
                            Self::from_formula(py, 0.0, name.as_str(), None)
                                .ok()
                                .map(|material| Cow::Owned(material))
                        );

                    match material {
                        Some(material) => {
                            let xi = wi / material.mass;
                            for cj in material.composition.iter() {
                                let Component { name: symbol, weight: wj } = cj;
                                let (symbol, element) = table.0
                                    .get_key_value(symbol.as_str()).unwrap();
                                let xj = wj / element.A * material.mass;
                                let xij = xi * xj;
                                weights
                                    .entry(symbol)
                                    .and_modify(|x| *x += xij)
                                    .or_insert(xij);
                            }
                            xi
                        },
                        None => {
                            let why = format!(
                                "unknown element, material or molecule '{}'",
                                name.as_str(),
                            );
                            return Err((KeyError, why))
                        },
                    }
                },
            };
            sum += xi;
        }
        let n = composition.len();
        let mut composition = Vec::<Component>::with_capacity(n);
        for (symbol, weight) in weights.iter() {
            composition.push(Component { name: symbol.to_string(), weight: weight / sum })
        }
        Self::from_elements(py, density, &composition, I)
    }

    pub fn from_elements(
        py: Python,
        density: f64,
        composition: &[Component], // Beware: mole fractions.
        #[allow(non_snake_case)]
        I: Option<f64>,
    ) -> Result<Self, ErrorData> {
        let table = ELEMENTS.get(py).unwrap();
        let n = composition.len();
        let mut mass = 0.0;
        let mut mass_composition = Vec::<Component>::with_capacity(n);
        for component in composition.iter() {
            let element = table.0.get::<str>(&component.name)
                .ok_or_else(|| {
                    let why = format!("unkown element '{}'", &component.name);
                    (KeyError, why)
                })?;
            let weight = component.weight * element.A;
            mass_composition.push(Component { name: component.name.clone(), weight });
            mass += weight;
        }
        for i in 0..n {
            mass_composition[i].weight /= mass;
        }

        let material = Self::new(
            density,
            mass,
            mass_composition,
            I,
        );
        Ok(material)
    }

    pub fn from_formula(
        py: Python,
        density: f64,
        formula: &str,
        #[allow(non_snake_case)]
        I: Option<f64>,
    ) -> Result<Self, ErrorData> {
        let table = ELEMENTS.get(py).unwrap();
        let re = Regex::new(r"([A-Z][a-z]?)([0-9]*)").unwrap();
        let mut composition = Vec::<Component>::new();
        let mut sum = 0.0;
        for captures in re.captures_iter(formula) {
            let symbol = captures.get(1).unwrap().as_str();
            if !table.0.contains_key(symbol) {
                let why = format!("unknown element '{}'", symbol);
                return Err((KeyError, why))
            }
            let weight = captures.get(2).unwrap().as_str();
            let weight: f64 = if weight.len() == 0 {
                1.0
            } else {
                weight.parse::<f64>()
                    .map_err(|_| {
                        let why = format!(
                            "could not parse weight ('{}') for '{}'",
                            weight,
                            symbol,
                        );
                        (ValueError, why)
                    })?
            };
            composition.push(Component { name: symbol.to_string(), weight });
            sum += weight;
        }
        if (sum - 1.0).abs() > 1E-09 {
            for component in composition.iter_mut() {
                component.weight /= sum;
            }
        }
        Material::from_elements(py, density, &composition, I)
    }
}

impl Component {
    pub fn new(name: String, weight: f64) -> Self {
        Self { name, weight }
    }
}


// ===============================================================================================
//
// Data relative to a set of materials.
//
// ===============================================================================================

#[derive(Debug, PartialEq)]
pub struct MaterialsData (pub HashMap<String, Material>);

impl MaterialsData {
    pub fn new() -> Self {
        Self (HashMap::new())
    }

    pub fn from_file<P: AsRef<Path>>(py: Python, path: P) -> PyResult<Self> {
        Toml::load_dict(py, path.as_ref())?
            .try_into()
    }

    pub fn sync(&self, py: Python, tag: &str) -> PyResult<bool> {
        let materials_cache = cache::path()?
            .join("materials");
        let cached = materials_cache
            .join(format!("{}.toml", tag));

        let update_materials_definition = || -> PyResult<()> {
            match std::fs::read_dir(&materials_cache) {
                Ok(content) => {
                    // Remove any cached data.
                    for entry in content {
                        if let Ok(entry) = entry {
                            if let Some(filename) = entry.file_name().to_str() {
                                if filename.starts_with(tag) &&
                                   filename.ends_with(".pumas") {
                                    std::fs::remove_file(&entry.path())?;
                                }
                            }
                        }
                    }
                },
                Err(_) => std::fs::create_dir_all(&materials_cache)?,
            }

            std::fs::write(&cached, self.to_toml())?;
            Ok(())
        };

        let updated = if cached
            .try_exists()
            .unwrap_or(false) {
            // Compare the materials definitions.
            let cached = MaterialsData::from_file(py, &cached)?;
            if cached != *self {
                update_materials_definition()?;
                true
            } else {
                false
            }
        } else {
            // Apply the materials definition.
            update_materials_definition()?;
            true
        };
        Ok(updated)
    }
}
