use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::convert::{Ellipsoid, Mode};
use crate::utils::export::Export;
use crate::utils::float::f64x3;
use crate::utils::namespace::Namespace;
use derive_more::{AsMut, AsRef, From};
use pyo3::prelude::*;
use ::std::ffi::c_int;
use ::std::pin::Pin;


#[repr(C)]
pub struct Recorder {
    base: danton::Recorder,
    pub event: usize,
    pub random_index: u128,
    pub weight: f64,
    pub mode: Mode,
    pub decay: bool,
    pub ellipsoid: Ellipsoid,
    grammages: Option<Vec<f64>>,
    primaries: Option<Vec<Primary>>,
    secondaries: Option<Vec<Secondary>>,
    products: Option<Vec<Product>>,
    vertices: Option<Vec<Vertex>>,
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct GrammagesExport (Export<f64>);

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Primary {
    pub event: usize,
    pub particle: Particle,
    pub random_index: [u64; 2],
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct PrimariesExport (Export<Primary>);

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Secondary {
    event: usize,
    particle: Particle,
    random_index: [u64; 2],
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct SecondariesExport (Export<Secondary>);

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Product {
    event: usize,
    pid: i32,
    momentum: f64,
    theta: f64,
    phi: f64,
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct ProductsExport (Export<Product>);

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Vertex {
    event: usize,
    particle: Particle,
    generation: i32,
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct VerticesExport (Export<Vertex>);

impl Recorder {
    pub fn new() -> Pin<Box<Self>> {
        let base = danton::Recorder {
            record_event: Some(Self::record_event),
            record_grammage: Some(Self::record_grammage),
        };
        let recorder = Self {
            base,
            event: 0,
            random_index: 0,
            weight: 1.0,
            mode: Mode::Backward,
            decay: true,
            ellipsoid: Ellipsoid::Prem81,
            grammages: None,
            primaries: None,
            secondaries: None,
            products: None,
            vertices: None
        };
        Box::pin(recorder)
    }

    pub fn base(recorder: &mut Pin<Box<Self>>) -> &mut danton::Recorder {
        &mut recorder.base
    }

    pub fn clear(&mut self) {
        self.grammages = None;
        self.primaries = None;
        self.secondaries = None;
        self.products = None;
        self.vertices = None;
    }

    pub fn export(&mut self, py: Python, steps: Option<PyObject>) -> PyResult<PyObject> {
        if let Mode::Grammage = self.mode {
            let grammages = match self.grammages.take() {
                None => Export::<f64>::empty(py)?,
                Some(grammages) => Export::export::<GrammagesExport>(py, grammages)?,
            };
            let result = match steps {
                None => grammages,
                Some(steps) => {
                    Namespace::new(py, &[
                        ("grammages", grammages),
                        ("steps", steps),
                    ])?.unbind()
                },
            };
            return Ok(result)
        }

        let products = match self.products.take() {
            None => Export::<Product>::empty(py)?,
            Some(products) => Export::export::<ProductsExport>(py, products)?,
        };
        let vertices = match self.vertices.take() {
            None => Export::<Vertex>::empty(py)?,
            Some(vertices) => Export::export::<VerticesExport>(py, vertices)?,
        };
        let result = match self.mode {
            Mode::Backward => {
                let primaries = match self.primaries.take() {
                    None => Export::<Primary>::empty(py)?,
                    Some(primaries) => Export::export::<PrimariesExport>(py, primaries)?,
                };
                if self.decay {
                    match steps {
                        None => {
                            Namespace::new(py, &[
                                ("primaries", primaries),
                                ("products", products),
                                ("vertices", vertices),
                            ])?.unbind()
                        },
                        Some(steps) => {
                            Namespace::new(py, &[
                                ("primaries", primaries),
                                ("products", products),
                                ("vertices", vertices),
                                ("steps", steps),
                            ])?.unbind()
                        },
                    }
                } else {
                    match steps {
                        None => {
                            Namespace::new(py, &[
                                ("primaries", primaries),
                                ("vertices", vertices),
                            ])?.unbind()
                        },
                        Some(steps) => {
                            Namespace::new(py, &[
                                ("primaries", primaries),
                                ("vertices", vertices),
                                ("steps", steps),
                            ])?.unbind()
                        },
                    }
                }
            },
            Mode::Forward => {
                let secondaries = match self.secondaries.take() {
                    None => Export::<Secondary>::empty(py)?,
                    Some(secondaries) => Export::export::<SecondariesExport>(py, secondaries)?,
                };
                if self.decay {
                    match steps {
                        None => {
                            Namespace::new(py, &[
                                ("secondaries", secondaries),
                                ("products", products),
                                ("vertices", vertices),
                            ])?.unbind()
                        },
                        Some(steps) => {
                            Namespace::new(py, &[
                                ("secondaries", secondaries),
                                ("products", products),
                                ("vertices", vertices),
                                ("steps", steps),
                            ])?.unbind()
                        },
                    }
                } else {
                    match steps {
                        None => {
                            Namespace::new(py, &[
                                ("secondaries", secondaries),
                                ("vertices", vertices),
                            ])?.unbind()
                        },
                        Some(steps) => {
                            Namespace::new(py, &[
                                ("secondaries", secondaries),
                                ("vertices", vertices),
                                ("steps", steps),
                            ])?.unbind()
                        },
                    }
                }
            },
            Mode::Grammage => unreachable!(),
        };

        Ok(result)
    }

    unsafe extern "C" fn record_event(
        _context: *mut danton::Context,
        recorder: *mut danton::Recorder,
        event: *const danton::Event,
    ) -> c_int {
        let recorder = unsafe { &mut *(recorder as *mut Self) };
        let event = unsafe { &*event };

        if let Mode::Backward = recorder.mode {
            let primary = {
                let state = unsafe { &*event.primary };
                let mut particle: Particle = (state, recorder.ellipsoid).into();
                particle.weight = event.weight * recorder.weight;
                Primary {
                    event: recorder.event,
                    particle,
                    random_index: [
                        (recorder.random_index >> 64) as u64,
                        recorder.random_index as u64
                    ],
                }
            };
            match recorder.primaries.as_mut() {
                None => recorder.primaries = Some(vec![primary]),
                Some(primaries) => primaries.push(primary),
            };
        }

        if let Mode::Forward = recorder.mode {
            let secondary = {
                let state = unsafe { &*event.secondary };
                let mut particle: Particle = (state, recorder.ellipsoid).into();
                particle.weight = event.weight * recorder.weight;
                Secondary {
                    event: recorder.event,
                    particle,
                    random_index: [
                        (recorder.random_index >> 64) as u64,
                        recorder.random_index as u64
                    ],
                }
            };
            match recorder.secondaries.as_mut() {
                None => recorder.secondaries = Some(vec![secondary]),
                Some(secondaries) => secondaries.push(secondary),
            };
        }

        if !event.vertex.is_null() {
            let vertex = {
                let state = unsafe { &*event.vertex };
                let particle: Particle = (state, recorder.ellipsoid).into();
                Vertex {
                    event: recorder.event,
                    particle,
                    generation: event.generation,
                }
            };
            match recorder.vertices.as_mut() {
                None => recorder.vertices = Some(vec![vertex]),
                Some(vertices) => vertices.push(vertex),
            };
        }

        if recorder.decay {
            let state = unsafe { &*event.secondary };
            let u: f64x3 = (&state.direction).into();
            let (v, w) = u.covectors();
            for i in 0..event.n_products {
                let product = unsafe { &*event.product.offset(i as isize) };
                let p: f64x3 = (&product.momentum).into();
                let momentum = p.norm();
                let cos_theta = (p.dot(u) / momentum).clamp(-1.0, 1.0);
                let theta = cos_theta.acos() * 180.0 / ::std::f64::consts::PI;
                let phi = if theta == 0.0 {
                    0.0
                } else {
                    let x = p.dot(v);
                    let y = p.dot(w);
                    y.atan2(x) * 180.0 / ::std::f64::consts::PI
                };
                let product = Product {
                    event: recorder.event,
                    pid: product.pid,
                    momentum,
                    theta,
                    phi,
                };
                match recorder.products.as_mut() {
                    None => recorder.products = Some(vec![product]),
                    Some(products) => products.push(product),
                };
            }
        }

        danton::SUCCESS
    }

    unsafe extern "C" fn record_grammage(
        _context: *mut danton::Context,
        recorder: *mut danton::Recorder,
        grammage: *const danton::Grammage,
    ) -> c_int {
        let recorder = unsafe { &mut *(recorder as *mut Self) };
        let grammage = unsafe { &*grammage };

        match recorder.grammages.as_mut() {
            None => recorder.grammages = Some(vec![grammage.value]),
            Some(grammages) => grammages.push(grammage.value),
        };

        danton::SUCCESS
    }
}
