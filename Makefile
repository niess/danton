DANTON_DIR := $(abspath .)
CFLAGS := -O2 -std=c99 -pedantic -Wall -DDANTON_DIR="\"$(DANTON_DIR)\""
INCLUDE := -Iinclude -Ient/include -Ipumas/include -Ialouette/include
LIBS := -Lent/lib -lent -Lpumas/lib -lpumas -Lalouette/lib -lalouette -lm

.PHONY: bin clean lib libclean

bin: bin/danton

clean:
	@rm -rf bin lib/libdanton.so

bin/danton: src/danton-x.c lib/libdanton.so
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< lib/libdanton.so $(LIBS)

lib/libdanton.so: src/danton.c include/danton.h lib
	@gcc -o $@ $(CFLAGS) $(INCLUDE) -fPIC -Iinclude -shared $< $(LIBS)

lib:
	@make -C "pumas"
	@make -C "ent"
	@make -C "alouette"

libclean:
	@make -C "pumas" clean
	@make -C "ent" clean
	@make -C "alouette" clean
