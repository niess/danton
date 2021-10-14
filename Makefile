# Generic options
DANTON_DEFAULT_PDF := $(abspath deps/ent/share/pdf/CT14nlo_0000.dat)
DANTON_DEFAULT_MDF := $(abspath share/materials/materials.xml)
DANTON_DEFAULT_DEDX := $(abspath share/materials/dedx)

# Compiler flags
CFLAGS := -O3 -std=c99 -Wall
FFLAGS := -O2 -fno-second-underscore -fno-backslash -fno-automatic             \
	-ffixed-line-length-132

# OSX additional flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
        CFLAGS += -L"$(shell dirname `gfortran --print-file-name libgfortran.dylib`)"
        CFLAGS += -Wno-tautological-compare
endif

# Main build targets
.PHONY: all bin clean lib

all: bin lib

bin: bin/danton

clean:
	@rm -rf bin build lib/*.so lib/*.a

lib: lib/libdanton.so

bin/danton: src/danton-x.c lib/libdanton.so
	@mkdir -p bin
	@$(CC) -o $@ $(CFLAGS) $(INCLUDE) $<                                   \
		-Llib -ldanton -Wl,-rpath $(PWD)/lib

OBJS := $(addprefix build/,danton.lo text.lo discrete.lo powerlaw.lo)
# ALOUETTE
OBJS += $(addprefix build/,                                                    \
	formf.lo tauola.lo curr_cleo.lo pkorb.lo f3pi.lo tauola_extras.lo      \
	f3pi_rcht.lo funct_3pi.lo funct_rpt.lo value_parameter.lo FA1RCHL.lo   \
	ffwid3pi.lo initA1TabKKpi.lo wid_a1_fit.lo initA1Tab.lo                \
	wid_a1_fitKKpi.lo gaus_integr.lo gfact.lo frho_pi_belle.lo             \
	alouette.lo)
# ENT
OBJS += build/ent.lo
# JSMN-TEA
OBJS += build/jsmn.lo build/jsmn-tea.lo
# PUMAS
OBJS += build/pumas.lo
# TURTLE
OBJS += $(addprefix build/,                                                    \
	client.lo ecef.lo error.lo io.lo list.lo map.lo projection.lo stack.lo \
	stepper.lo tinydir.lo asc.lo geotiff16.lo grd.lo hgt.lo png16.lo)

lib/libdanton.so: $(OBJS)
	@$(CC) -o $@ $(CFLAGS) -shared $(OBJS) -lgfortran -lm -ldl

# Build DANTON
INCLUDE := -Iinclude -Ideps/ent/include -Ideps/pumas/include                   \
	-Ideps/alouette/include -Ideps/jsmn-tea/include                        \
	-Ideps/roar/include -Ideps/turtle/include -Ideps/turtle/src            \
	-Ideps/turtle/src/deps

define build_c
	@mkdir -p build
	@$(CC) -o $@ $(CFLAGS) $1 -fPIC -c $<
endef

build/danton.lo: src/danton.c
	@$(call build_c,-DDANTON_DEFAULT_PDF="\"$(DANTON_DEFAULT_PDF)\""       \
		-DDANTON_DEFAULT_MDF="\"$(DANTON_DEFAULT_MDF)\""               \
		-DDANTON_DEFAULT_DEDX="\"$(DANTON_DEFAULT_DEDX)\"" $(INCLUDE))

build/%.lo: src/danton/primary/%.c
	@$(call build_c,$(INCLUDE))

build/%.lo: src/danton/recorder/%.c
	@$(call build_c,$(INCLUDE))

# Build ALOUETTE
build/%.lo: deps/alouette/src/%.c
	@$(call build_c,-Ideps/alouette/include)

define build_fortran
	@mkdir -p build
	@gfortran -o $@ $(FFLAGS) -fPIC -c $<
endef

TAUSRC = deps/alouette/src/tauola

build/%.lo: $(TAUSRC)/%.f
	@$(call build_fortran)

build/%.lo: $(TAUSRC)/new-currents/RChL-currents/rcht_3pi/%.f
	@$(call build_fortran)

build/%.lo: $(TAUSRC)/new-currents/RChL-currents/rcht_common/%.f
	@$(call build_fortran)
	
build/%.lo: $(TAUSRC)/new-currents/other-currents/%.f
	@$(call build_fortran)
	
# Build ENT
build/%.lo: deps/ent/src/%.c
	@$(call build_c,-Ideps/ent/include)

# Build JSMN-TEA
build/%.lo: deps/jsmn-tea/src/%.c deps/jsmn-tea/include/%.h
	@$(call build_c,-Ideps/jsmn-tea/include -Ideps/turtle/src/deps         \
		-Ideps/roar/include -DROAR_IMPLEMENTATION)
	
# Build PUMAS
build/%.lo: deps/pumas/src/%.c
	@$(call build_c,-Ideps/pumas/include)

# Build TURTLE
build/%.lo: deps/turtle/src/turtle/%.c deps/turtle/src/turtle/%.h
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src)

build/%.lo: deps/turtle/src/turtle/%.c
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src)

build/%.lo: deps/turtle/src/turtle/io/%.c
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src)

build/%.lo: deps/turtle/src/%.c deps/turtle/include/%.h
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src)

build/%.lo: deps/turtle/src/deps/%.c deps/turtle/src/deps/%.h
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src)
