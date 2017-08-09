DANTON_DIR := $(abspath .)
CFLAGS := -O2 -std=c99 -pedantic -Wall -DDANTON_DIR="\"$(DANTON_DIR)\""
INCLUDE := -Iinclude -Ient/include -Ipumas/include -Ialouette/include
LIBS := -Lent/lib -lent -Lpumas/lib -lpumas -Lalouette/lib -lalouette -lm
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
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< lib/libdanton.so $(LIBS)

lib/libdanton.so: $(DANTON_SRC) $(DANTON_INC) lib
	@gcc -o $@ $(CFLAGS) $(INCLUDE) -fPIC -Iinclude -shared $(DANTON_SRC) $(LIBS)

lib:
	@make -C "pumas"
	@make -C "ent"
	@make -C "alouette"

libclean:
	@make -C "pumas" clean
	@make -C "ent" clean
	@make -C "alouette" clean
