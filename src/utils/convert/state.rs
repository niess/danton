use crate::bindings::{danton, turtle};
use crate::utils::convert::Geodesic;
use crate::simulation::particles::Particle;


impl From<(&danton::State, Geodesic)> for Particle {
    fn from(value: (&danton::State, Geodesic)) -> Self {
        let (state, geodesic) = value;
        let mut latitude: f64 = 0.0;
        let mut longitude: f64 = 0.0;
        let mut altitude: f64 = 0.0;
        let mut azimuth: f64 = 0.0;
        let mut elevation: f64 = 0.0;
        match geodesic {
            Geodesic::Prem => {
                unimplemented!() // XXX Implement this case.
            },
            _ => unsafe {
                turtle::ecef_to_geodetic(
                    &state.position as *const f64,
                    &mut latitude,
                    &mut longitude,
                    &mut altitude
                );
                turtle::ecef_to_horizontal(
                    latitude,
                    longitude,
                    &state.direction as *const f64,
                    &mut azimuth,
                    &mut elevation,
                );
            }
        }
        fn tozero(f: f64) -> f64 {
            if f.abs() <= 3.0 * ::std::f64::EPSILON {
                0.0
            } else {
                f
            }
        }
        Self {
            pid: state.pid,
            energy: state.energy,
            latitude: tozero(latitude),
            longitude: tozero(longitude),
            altitude: tozero(altitude),
            azimuth: tozero(azimuth),
            elevation: tozero(elevation),
        }
    }
}
