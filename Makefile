# Installation prefix and project name
PREFIX= $(abspath .)
NAME=   danton

# Defaul data paths
DANTON_DEDX=  share/$(NAME)/materials/dedx
DANTON_DUMP=  share/$(NAME)/materials/materials.pumas
DANTON_GEOID= share/$(NAME)/geoid/egm96.png
DANTON_MDF=   share/$(NAME)/materials/materials.xml
DANTON_PDF=   share/$(NAME)/pdf/CT14nlo_0000.dat

# Compiler flags
CFLAGS= -O3 -Wall
FFLAGS= -O3 -fno-second-underscore -fno-backslash -fno-automatic               \
	-ffixed-line-length-132 -std=legacy

# OSX additional flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SOEXT=  dylib
	SHARED= -dynamiclib -Wl,-install_name,@rpath/lib$(NAME).$(SOEXT)
	RPATH=  -Wl,-rpath,@loader_path/../lib
else
	SOEXT=  so
	SHARED= -shared
	RPATH=  '-Wl,-rpath,$$ORIGIN/../lib'
endif

BINDIR= $(PREFIX)/bin
EXEC=   $(BINDIR)/$(NAME)
LIBDIR= $(PREFIX)/lib
LIB=    $(LIBDIR)/lib$(NAME).$(SOEXT)

# Main build targets
.PHONY: all bin clean lib

all: bin lib

bin: $(EXEC)

clean:
	@rm -f $(EXEC) $(LIB) $(PREFIX)/$(DANTON_DUMP)
	@if [ -d $(BINDIR) ]; \
	    then rmdir --ignore-fail-on-non-empty $(BINDIR); fi
	@rm -rf build

lib: $(LIB)

$(EXEC): src/danton-x.c $(LIB)
	@mkdir -p $(BINDIR)
	@$(CC) -o $@ $(CFLAGS) $(INCLUDE) $< -L$(LIBDIR) -l$(NAME) $(RPATH)

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

$(LIB): $(OBJS)
	@$(CC) -o $@ $(CFLAGS) $(SHARED) $(OBJS) -lm -ldl

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
	@$(call build_c,-DDANTON_PREFIX="\"$(PREFIX)\""                        \
		-DDANTON_PDF="\"$(DANTON_PDF)\""                               \
		-DDANTON_MDF="\"$(DANTON_MDF)\""                               \
		-DDANTON_DEDX="\"$(DANTON_DEDX)\""                             \
		-DDANTON_GEOID="\"$(DANTON_GEOID)\""                           \
		-DDANTON_DUMP="\"$(DANTON_DUMP)\""                             \
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
