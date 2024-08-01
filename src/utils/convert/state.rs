use crate::bindings::danton;
use crate::simulation::particles::Particle;
use crate::utils::coordinates::{GeodeticCoordinates, HorizontalCoordinates};
use crate::utils::convert::Geodesic;


// ===============================================================================================
//
// Convert Danton Monte Carlo states (i.e. ECEF to geodetic coordinates).
//
// ===============================================================================================

impl From<(&danton::State, Geodesic)> for Particle {
    fn from(value: (&danton::State, Geodesic)) -> Self {
        let (state, geodesic) = value;
        let geodetic = GeodeticCoordinates::from_ecef(&state.position, geodesic);
        let horizontal = HorizontalCoordinates::from_ecef(&state.direction, geodesic, &geodetic);

        fn tozero(f: f64) -> f64 {
            if f.abs() <= ::std::f32::EPSILON as f64 {
                0.0
            } else {
                f
            }
        }

        Self {
            pid: state.pid,
            energy: state.energy,
            latitude: tozero(geodetic.latitude),
            longitude: tozero(geodetic.longitude),
            altitude: tozero(geodetic.altitude),
            azimuth: tozero(horizontal.azimuth),
            elevation: tozero(horizontal.elevation),
            weight: 1.0,
        }
    }
}
