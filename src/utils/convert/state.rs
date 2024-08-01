use crate::bindings::{danton, turtle};
use crate::utils::convert::Geodesic;
use crate::simulation::particles::Particle;


// ===============================================================================================
//
// Convert Danton Monte Carlo states (i.e. ECEF to geodetic coordinates).
//
// ===============================================================================================

impl From<(&danton::State, Geodesic)> for Particle {
    fn from(value: (&danton::State, Geodesic)) -> Self {
        let (state, geodesic) = value;
        let mut latitude: f64 = 0.0;
        let mut longitude: f64 = 0.0;
        let mut altitude: f64 = 0.0;
        let mut azimuth: f64 = 0.0;
        let mut elevation: f64 = 0.0;
        match geodesic {
            Geodesic::Prem => { // XXX Use GeodeticCoordinates?
                const DEG: f64 = 180.0 / ::std::f64::consts::PI;
                let r = &state.position;
                let rho2 = r[0] * r[0] + r[1] * r[1];
                let rho = rho2.sqrt();
                let theta = rho.atan2(r[2]);
                let phi = r[1].atan2(r[0]);
                latitude = 90.0 - theta * DEG;
                longitude = phi * DEG;
                altitude = (rho2 + r[2] * r[2]).sqrt() - danton::PREM_EARTH_RADIUS;

                let st = theta.sin();
                let ct = theta.cos();
                let sp = phi.sin();
                let cp = phi.cos();
                #[allow(non_snake_case)]
                let R = [[-sp, cp, 0. ], [-ct * cp, -ct * sp, st], [st * cp, st * sp, ct]];

                let u = &state.direction;
                let vx = R[0][0] * u[0] + R[0][1] * u[1] + R[0][2] * u[2];
                let vy = R[1][0] * u[0] + R[1][1] * u[1] + R[1][2] * u[2];
                let vz = R[2][0] * u[0] + R[2][1] * u[1] + R[2][2] * u[2];
                let rho = (vx * vx + vy * vy).sqrt();
                elevation = 90.0 - rho.atan2(vz) * DEG;
                azimuth = if elevation.abs() == 90.0 {
                    0.0
                } else {
                    90.0 - vy.atan2(vx) * DEG
                };
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
            if f.abs() <= ::std::f32::EPSILON as f64 {
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
            weight: 1.0,
        }
    }
}
