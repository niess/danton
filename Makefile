DANTON_DIR := $(abspath .)
CFLAGS := -O2 -std=c99 -pedantic -Wall
INCLUDE := -Iinclude -Imodules/ent/include -Imodules/pumas/include \
	-Imodules/alouette/include -Imodules/jsmn
LIBS := -Lmodules/ent/lib -lent -Lmodules/pumas/lib -lpumas \
	-Lmodules/alouette/lib -lalouette -lm -Lmodules/jsmn -ljsmn
DANTON_SRC := src/danton.c src/danton/recorder/text.c \
	src/danton/primary/discrete.c src/danton/primary/powerlaw.c
DANTON_INC := include/danton.h include/danton/recorder/text.h \
	include/danton/primary/discrete.h include/danton/primary/powerlaw.h

.PHONY: bin clean lib libclean

bin: bin/danton

clean:
	@rm -rf bin lib/libdanton.so

bin/danton: src/danton-x.c lib/libdanton.so
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< -Llib -ldanton $(LIBS)

lib/libdanton.so: $(DANTON_SRC) $(DANTON_INC) lib
	@gcc -o $@ $(CFLAGS) -DDANTON_DIR="\"$(DANTON_DIR)\"" $(INCLUDE) -fPIC \
		-Iinclude -shared $(DANTON_SRC) $(LIBS)

export CFLAGS

lib:
	@make -C "modules/pumas"
	@make -C "modules/ent"
	@make -C "modules/alouette"
	@make -C "modules/jsmn"

libclean:
	@make -C "modules/pumas" clean
	@make -C "modules/ent" clean
	@make -C "modules/alouette" clean
	@make -C "modules/jsmn" clean
