use crate::bindings::danton;
use std::pin::Pin;


impl danton::Primary {
    pub fn new(min: f64, max: f64) -> Pin<Box<Self>> {
        let energy = [ min, max ];
        Box::pin(Self { flux: Some(Self::flux), energy })
    }

    unsafe extern "C" fn flux(_primary: *mut Self, _energy: f64) -> f64 {
        1.0
    }
}
