use crate::bindings::{danton, turtle};
use crate::utils::convert::Geodesic;


// ===============================================================================================
//
// Geodetic coordinates system.
//
// ===============================================================================================

pub struct GeodeticCoordinates {
    pub latitude: f64,
    pub longitude: f64,
    pub altitude: f64,
}

impl GeodeticCoordinates {
    pub fn from_ecef(position: &[f64; 3], geodesic: Geodesic) -> Self {
        match geodesic {
            Geodesic::Prem => {
                const DEG: f64 = 180.0 / ::std::f64::consts::PI;
                let r = position;
                let rho2 = r[0] * r[0] + r[1] * r[1];
                let rho = rho2.sqrt();
                let theta = rho.atan2(r[2]);
                let phi = r[1].atan2(r[0]);
                let latitude = 90.0 - theta * DEG;
                let longitude = phi * DEG;
                let altitude = (rho2 + r[2] * r[2]).sqrt() - danton::PREM_EARTH_RADIUS;
                Self { latitude, longitude, altitude }
            },
            _ => {
                let mut latitude = 0.0;
                let mut longitude = 0.0;
                let mut altitude = 0.0;
                unsafe {
                    turtle::ecef_to_geodetic(
                        position as *const f64,
                        &mut latitude,
                        &mut longitude,
                        &mut altitude
                    );
                }
                Self { latitude, longitude, altitude }
            },
        }
    }

    pub fn to_ecef(&self, geodesic: Geodesic) -> [f64; 3] {
        match geodesic {
            Geodesic::Prem => {
                const RAD: f64 = ::std::f64::consts::PI / 180.0;
                let r = danton::PREM_EARTH_RADIUS + self.altitude;
                let theta = (90.0 - self.latitude) * RAD;
                let phi = self.longitude * RAD;
                let cos_theta = theta.cos();
                let sin_theta = theta.sin();
                [
                    r * sin_theta * phi.cos(),
                    r * sin_theta * phi.sin(),
                    r * cos_theta,
                ]
            },
            _ => {
                let mut position = [0_f64; 3];
                unsafe {
                    turtle::ecef_from_geodetic(
                        self.latitude,
                        self.longitude,
                        self.altitude,
                        (&mut position) as *mut f64,
                    );
                }
                position
            },
        }
    }
}


// ===============================================================================================
//
// Horizontal angular coordinates.
//
// ===============================================================================================

pub struct HorizontalCoordinates {
    pub azimuth: f64,
    pub elevation: f64,
}

impl HorizontalCoordinates {
    const DEG: f64 = 180.0 / ::std::f64::consts::PI;
    const RAD: f64 = ::std::f64::consts::PI / 180.0;

    pub fn from_ecef(
        direction: &[f64; 3],
        geodesic: Geodesic,
        origin: &GeodeticCoordinates
    ) -> Self {
        match geodesic {
            Geodesic::Prem => {
                #[allow(non_snake_case)]
                let R = Self::rotation(origin);
                let u = direction;
                let vx = R[0][0] * u[0] + R[0][1] * u[1] + R[0][2] * u[2];
                let vy = R[1][0] * u[0] + R[1][1] * u[1] + R[1][2] * u[2];
                let vz = R[2][0] * u[0] + R[2][1] * u[1] + R[2][2] * u[2];
                let rho = (vx * vx + vy * vy).sqrt();
                let elevation = 90.0 - rho.atan2(vz) * Self::DEG;
                let azimuth = if elevation.abs() == 90.0 {
                    0.0
                } else {
                    90.0 - vy.atan2(vx) * Self::DEG
                };
                Self { azimuth, elevation }
            },
            _ => {
                let mut azimuth: f64 = 0.0;
                let mut elevation: f64 = 0.0;
                unsafe {
                    turtle::ecef_to_horizontal(
                        origin.latitude,
                        origin.longitude,
                        direction as *const f64,
                        &mut azimuth,
                        &mut elevation,
                    );
                }
                Self { azimuth, elevation }
            },
        }
    }

    pub fn to_ecef(
        &self,
        geodesic: Geodesic,
        origin: &GeodeticCoordinates
    ) -> [f64; 3] {
        match geodesic {
            Geodesic::Prem => {
                #[allow(non_snake_case)]
                let R = Self::rotation(origin);
                let theta = (90.0 - self.elevation) * Self::RAD;
                let phi = (90.0 - self.azimuth) * Self::RAD;
                let cos_theta = theta.cos();
                let sin_theta = theta.sin();
                let u = [ sin_theta * phi.cos(), sin_theta * phi.sin(), cos_theta ];
                let vx = R[0][0] * u[0] + R[1][0] * u[1] + R[2][0] * u[2];
                let vy = R[0][1] * u[0] + R[1][1] * u[1] + R[2][1] * u[2];
                let vz = R[0][2] * u[0] + R[1][2] * u[1] + R[2][2] * u[2];
                [ vx, vy, vz ]
            },
            _ => {
                let mut direction = [0.0; 3];
                unsafe {
                    turtle::ecef_from_horizontal(
                        origin.latitude,
                        origin.longitude,
                        self.azimuth,
                        self.elevation,
                        (&mut direction) as *mut f64,
                    );
                }
                direction
            },
        }
    }

    fn rotation(origin: &GeodeticCoordinates) -> [[f64;3];3] {
        let theta = (90.0 - origin.latitude) * Self::RAD;
        let phi = origin.longitude * Self::RAD;
        let st = theta.sin();
        let ct = theta.cos();
        let sp = phi.sin();
        let cp = phi.cos();
        [[-sp, cp, 0. ], [-ct * cp, -ct * sp, st], [st * cp, st * sp, ct]]
    }
}
