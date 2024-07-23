// XXX Use relative paths?

fn main() {
    let c_includes = [
        "deps/alouette/include",
        "deps/ent/include",
        "deps/pumas/include",
        "deps/turtle/include",
        "deps/turtle/src",
    ];

    let c_headers = [
        "src/danton.h",
        "deps/alouette/include/alouette.h",
        "deps/ent/include/ent.h",
        "deps/pumas/include/pumas.h",
        "deps/turtle/include/turtle.h",
        "deps/turtle/src/turtle/client.h",
        "deps/turtle/src/turtle/error.h",
        "deps/turtle/src/turtle/io.h",
        "deps/turtle/src/turtle/list.h",
        "deps/turtle/src/turtle/map.h",
        "deps/turtle/src/turtle/projection.h",
        "deps/turtle/src/turtle/stack.h",
        "deps/turtle/src/turtle/stepper.h",
        "deps/turtle/src/deps/jsmn.h",
        "deps/turtle/src/deps/png.h",
        "deps/turtle/src/deps/tiffio.h",
        "deps/turtle/src/deps/tinydir.h",
    ];

    let c_sources = [
        "src/danton.c",
        "deps/alouette/src/alouette.c",
        "deps/ent/src/ent.c",
        "deps/pumas/src/pumas.c",
        "deps/turtle/src/turtle/client.c",
        "deps/turtle/src/turtle/ecef.c",
        "deps/turtle/src/turtle/error.c",
        "deps/turtle/src/turtle/io.c",
        "deps/turtle/src/turtle/list.c",
        "deps/turtle/src/turtle/map.c",
        "deps/turtle/src/turtle/projection.c",
        "deps/turtle/src/turtle/stack.c",
        "deps/turtle/src/turtle/stepper.c",
        "deps/turtle/src/turtle/io/asc.c",
        "deps/turtle/src/turtle/io/geotiff16.c",
        "deps/turtle/src/turtle/io/grd.c",
        "deps/turtle/src/turtle/io/hgt.c",
        "deps/turtle/src/turtle/io/png16.c",
        "deps/turtle/src/deps/jsmn.c",
        "deps/turtle/src/deps/tinydir.c",
    ];

    let f_sources = [
        "deps/alouette/src/tauola.f",
    ];

    cc::Build::new()
        .includes(c_includes)
        .files(c_sources)
        .define("DANTON_PREFIX",  "\"danton\"")
        .define("DANTON_DEDX",    "\"data/materials/dedx\"")
        .define("DANTON_DUMP",    "\"data/materials/materials.pumas\"")
        .define("DANTON_GEOID",   "\"data/geoid/egm96.png\"")
        .define("DANTON_MDF",     "\"data/materials/materials.xml\"")
        .define("DANTON_PDF_DIR", "\"data/pdf\"")
        .define("DANTON_CS_DIR",  "\"data/cs\"")
        .warnings(false)
        .compile("danton-c");

    for path in c_headers.iter()
        .chain(c_sources.iter()) {
        println!("cargo:rerun-if-changed={}", path);
    }

    cc::Build::new()
        .compiler("gfortran")
        .flag("-fno-second-underscore")
        .flag("-fno-backslash")
        .flag("-fno-automatic")
        .flag("-ffixed-line-length-132")
        .flag("-std=legacy")
        .warnings(false)
        .files(f_sources)
        .compile("danton-f");

    for path in f_sources.iter() {
        println!("cargo:rerun-if-changed={}", path);
    }
}
