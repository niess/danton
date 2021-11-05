# Installation prefix and project name
PREFIX= $(abspath .)
NAME=   danton

# Defaul data paths
DANTON_DEDX=  share/$(NAME)/materials/dedx
DANTON_DUMP=  share/$(NAME)/materials/materials.pumas
DANTON_GEOID= share/$(NAME)/geoid/egm96.png
DANTON_MDF=   share/$(NAME)/materials/materials.xml
DANTON_PDF=   share/$(NAME)/pdf/CT14nlo_0000.dat

# Compiler options
CC=gcc
FC=gcc
CFLAGS= -O3 -Wall -std=c99
FFLAGS= -O3 -fno-second-underscore -fno-backslash -fno-automatic               \
	-ffixed-line-length-132 -std=legacy

# OS dependent flags
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
	SOEXT=  dylib
	SHARED= -dynamiclib -Wl,-install_name,@rpath/lib$(NAME).$(SOEXT)
	RPATH=  -Wl,-rpath,@loader_path/../lib
else
	SOEXT=  so
	SHARED= -shared
	RPATH=  '-Wl,-rpath,$$ORIGIN/../lib'

	TINYDIR_CFLAGS= -Wno-restrict
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

# Build the danton executable
EXE_INCLUDES= -Iinclude -Ideps/jsmn-tea/include -Ideps/roar/include            \
              -Ideps/turtle/src/deps
EXE_SRCS=     src/danton-x.c

$(EXEC): $(EXE_SRCS) $(LIB)
	@mkdir -p $(BINDIR)
	@$(CC) -o $@ $(CFLAGS) $(EXE_INCLUDES) $(EXE_SRCS)                     \
               -L$(LIBDIR) -l$(NAME) $(RPATH)

# Build the DANTON library
LIB_OBJS=  $(addprefix build/,danton.lo text.lo discrete.lo powerlaw.lo)
# ALOUETTE
LIB_OBJS+= $(addprefix build/,tauola.lo alouette.lo)
# ENT
LIB_OBJS+= build/ent.lo
# JSMN-TEA
LIB_OBJS+= build/jsmn.lo build/jsmn-tea.lo
# PUMAS
LIB_OBJS += build/pumas.lo
# TURTLE
LIB_OBJS += $(addprefix build/,                                                \
            client.lo ecef.lo error.lo io.lo list.lo map.lo projection.lo      \
            stack.lo stepper.lo tinydir.lo asc.lo geotiff16.lo grd.lo hgt.lo   \
            png16.lo)

$(LIB): $(LIB_OBJS)
	@mkdir -p $(LIBDIR)
	@$(CC) -o $@ $(CFLAGS) $(SHARED) $(LIB_OBJS) -lm -ldl

LIB_INCLUDES= -Iinclude -Ideps/ent/include -Ideps/pumas/include                \
              -Ideps/alouette/include -Ideps/jsmn-tea/include                  \
              -Ideps/roar/include -Ideps/turtle/include -Ideps/turtle/src      \
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
                $(LIB_INCLUDES))

build/%.lo: src/danton/primary/%.c
	@$(call build_c,$(LIB_INCLUDES))

build/%.lo: src/danton/recorder/%.c
	@$(call build_c,$(LIB_INCLUDES))

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
	@$(call build_c,-Ideps/turtle/include -Ideps/turtle/src                \
	                $(TINYDIR_CFLAGS))
