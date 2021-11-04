# Generic options
DANTON_DEFAULT_PDF=   $(abspath share/danton/pdf/CT14nlo_0000.dat)
DANTON_DEFAULT_MDF=   $(abspath share/danton/materials/materials.xml)
DANTON_DEFAULT_DEDX=  $(abspath share/danton/materials/dedx)
DANTON_DEFAULT_GEOID= $(abspath share/danton/geoid/egm96.png)

# Compiler flags
CFLAGS= -O3 -std=c99 -Wall
FFLAGS= -O3 -fno-second-underscore -fno-backslash -fno-automatic               \
	-ffixed-line-length-132 -std=legacy

# OSX additional flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	CFLAGS += -Wno-unused-command-line-argument
        CFLAGS += -Wno-tautological-compare
	SOEXT= dylib
else
	CFLAGS += -Wno-restrict
	SOEXT= so
endif

# Main build targets
.PHONY: all bin clean lib

all: bin lib

bin: bin/danton

clean:
	@rm -rf bin build lib/*.$(SOEXT)

lib: lib/libdanton.$(SOEXT)

bin/danton: src/danton-x.c lib/libdanton.$(SOEXT)
	@mkdir -p bin
	@$(CC) -o $@ $(CFLAGS) $(INCLUDE) $<                                   \
		-Llib -ldanton -Wl,-rpath $(PWD)/lib

OBJS := $(addprefix build/,danton.lo text.lo discrete.lo powerlaw.lo)
# ALOUETTE
OBJS += $(addprefix build/,tauola.lo alouette.lo)
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

lib/libdanton.$(SOEXT): $(OBJS)
	@$(CC) -o $@ $(CFLAGS) -shared $(OBJS) -lm -ldl

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
		-DDANTON_DEFAULT_DEDX="\"$(DANTON_DEFAULT_DEDX)\""             \
		-DDANTON_DEFAULT_GEOID="\"$(DANTON_DEFAULT_GEOID)\""           \
		$(INCLUDE))

build/%.lo: src/danton/primary/%.c
	@$(call build_c,$(INCLUDE))

build/%.lo: src/danton/recorder/%.c
	@$(call build_c,$(INCLUDE))

# Build ALOUETTE
build/%.lo: deps/alouette/src/%.c
	@$(call build_c,-Ideps/alouette/include)

define build_fortran
	@mkdir -p build
	@$(FC) -o $@ $(FFLAGS) -fPIC -c $<
endef

build/%.lo: deps/alouette/src/%.f
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
