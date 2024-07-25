use crate::bindings::danton;
use std::pin::Pin;


impl danton::Primary {
    pub const ENERGY_MIN: f64 = 1E+02;
    pub const ENERGY_MAX: f64 = 1E+12;

    pub fn new() -> Pin<Box<Self>> {
        let energy = [ Self::ENERGY_MIN, Self::ENERGY_MAX ];
        Box::pin(Self { flux: Some(Self::flux), energy })
    }

    pub fn configure(&mut self, energy: Option<f64>) {
        match energy {
            None => self.energy = [ Self::ENERGY_MIN, Self::ENERGY_MAX ],
            Some(energy) => self.energy = [ energy, energy ],
        }
    }

    unsafe extern "C" fn flux(_primary: *mut Self, _energy: f64) -> f64 {
        1.0
    }
}
