export DEPS_DIR := $(PWD)/deps
export LIB_DIR := lib
export PDF_DIR := $(abspath deps/ent/share/pdf)
export TURTLE_USE_PNG := 0
unexport CFLAGS

CFLAGS := -O2 -std=c99 -pedantic -Wall
INCLUDE := -Iinclude -I$(DEPS_DIR)/ent/include -I$(DEPS_DIR)/pumas/include     \
	-I$(DEPS_DIR)/alouette/include -I$(DEPS_DIR)/jsmn                      \
	-I$(DEPS_DIR)/jsmn-tea/include -I$(DEPS_DIR)/roar/include              \
	-I$(DEPS_DIR)/turtle/include
DANTON_SRC := src/danton.c src/danton/recorder/text.c                          \
	src/danton/primary/discrete.c src/danton/primary/powerlaw.c
DANTON_INC := include/danton.h include/danton/recorder/text.h                  \
	include/danton/primary/discrete.h include/danton/primary/powerlaw.h

.PHONY: bin clean lib

bin: lib bin/danton

clean:
	@rm -rf bin lib/*.so lib/*.a

lib: lib/libalouette.so lib/libdanton.so lib/libent.so lib/libjsmn-tea.a       \
	lib/libpumas.so lib/libturtle.so

bin/danton: src/danton-x.c
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< -Llib -ldanton -L$(LIB_DIR)         \
		-lalouette -ldanton -lent -ljsmn-tea -lpumas -lturtle -ltiff   \
                -lm

lib/libdanton.so: $(DANTON_SRC) $(DANTON_INC)
	@gcc -o $@ $(CFLAGS) -DPDF_DIR="\"$(PDF_DIR)\"" $(INCLUDE) -fPIC       \
		-Iinclude -shared $(DANTON_SRC)

define build_library
	echo "o Building $(1) ..."
	@$(MAKE) --directory="$(DEPS_DIR)/$(1)" -e $(2)
	@mkdir -p lib && mv $(DEPS_DIR)/$(1)/$(3)/*.$(4) lib
	@echo "--> Done"
endef

lib/lib%.so: deps/%/src deps/%/include
	@$(call build_library,$*,,lib,so)

lib/libjsmn-tea.a: deps/jsmn-tea/src deps/jsmn-tea/include
	@$(call build_library,jsmn-tea,static,lib,a)
