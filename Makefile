DANTON_DIR := $(abspath .)
CFLAGS := -O2 -std=c99 -pedantic -Wall -DDANTON_DIR="\"$(DANTON_DIR)\""
INCLUDE := -Ient/include -Ipumas/include -Ialouette/include
LIBS := -Lent/lib -lent -Lpumas/lib -lpumas -Lalouette/lib -lalouette

.PHONY: bin clean lib libclean

bin: bin/danton

clean:
	@rm -rf bin

bin/%: src/%.c
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< $(LIBS)

lib:
	@make -C "pumas"
	@make -C "ent"
	@make -C "alouette"

libclean:
	@make -C "pumas" clean
	@make -C "ent" clean
	@make -C "alouette" clean
