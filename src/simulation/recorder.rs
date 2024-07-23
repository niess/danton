use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::convert::Geodesic;
use crate::utils::export::Export;
use crate::utils::float::f64x3;
use crate::utils::tuple::NamedTuple;
use derive_more::{AsMut, AsRef, From};
use pyo3::prelude::*;
use ::std::ffi::c_int;
use ::std::pin::Pin;


#[repr(C)]
pub struct Recorder {
    base: danton::Recorder,
    pub event: usize,
    pub geodesic: Geodesic,
    primaries: Option<Vec<Primary>>,
    products: Option<Vec<Product>>,
    vertices: Option<Vec<Vertex>>,
}

#[repr(C)]
#[derive(Clone, Copy)]
pub struct Primary {
    event: usize,
    particle: Particle,
    weight: f64,
}

#[derive(AsMut, AsRef, From)]
#[pyclass(module="danton")]
struct PrimariesExport (Export<Primary>);

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
            geodesic: Geodesic::Prem,
            primaries: None,
            products: None,
            vertices: None
        };
        Box::pin(recorder)
    }

    pub fn base(recorder: &mut Pin<Box<Self>>) -> &mut danton::Recorder {
        &mut recorder.base
    }

    pub fn export(&mut self, py: Python) -> PyResult<PyObject> {
        let primaries = match self.primaries.take() {
            None => py.None(),
            Some(primaries) => Export::export::<PrimariesExport>(py, primaries)?,
        };
        let products = match self.products.take() {
            None => py.None(),
            Some(products) => Export::export::<ProductsExport>(py, products)?,
        };
        let vertices = match self.vertices.take() {
            None => py.None(),
            Some(vertices) => Export::export::<VerticesExport>(py, vertices)?,
        };

        static RESULT: NamedTuple<3> = NamedTuple::new(
            "Result", ["primaries", "products", "vertices"]);
        let result = RESULT.instance(py, (primaries, products, vertices))?.unbind();
        Ok(result)
    }

    unsafe extern "C" fn record_event(
        _context: *mut danton::Context,
        recorder: *mut danton::Recorder,
        event: *const danton::Event,
    ) -> c_int {
        let recorder = unsafe { &mut *(recorder as *mut Self) };
        let event = unsafe { &*event };

        let primary = {
            let state = unsafe { &*event.primary };
            let particle: Particle = (state, recorder.geodesic).into();
            Primary {
                event: recorder.event,
                particle,
                weight: event.weight,
            }
        };
        match recorder.primaries.as_mut() {
            None => recorder.primaries = Some(vec![primary]),
            Some(primaries) => primaries.push(primary),
        };

        let vertex = {
            let state = unsafe { &*event.vertex };
            let particle: Particle = (state, recorder.geodesic).into();
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

        let state = unsafe { &*event.primary };
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

        danton::SUCCESS
    }

    unsafe extern "C" fn record_grammage(
        _context: *mut danton::Context,
        _recorder: *mut danton::Recorder,
        _grammage: *const danton::Grammage,
    ) -> c_int {
        unimplemented!();
    }
}
