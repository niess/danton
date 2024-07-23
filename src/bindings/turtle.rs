#![allow(unused)]

#[link(name = "danton-c")]
extern "C" {
    #[link_name="turtle_ecef_to_geodetic"]
    pub fn ecef_to_geodetic(
        ecef: *const f64,
        latitude: *mut f64,
        longitude: *mut f64,
        altitude: *mut f64,
    );

    #[link_name="turtle_ecef_to_horizontal"]
    pub fn ecef_to_horizontal(
        latitude: f64,
        longitude: f64,
        direction: *const f64,
        azimuth: *mut f64,
        elevation: *mut f64,
    );
}
