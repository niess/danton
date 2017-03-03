CFLAGS := -O2 -std=c99 -pedantic -Wall
INCLUDE := -Ient/include -Ipumas/include -Itauola-c/include
LIBS := -Lent/lib -lent -Lpumas/lib -lpumas -Ltauola-c/lib -ltauola-c -lm

.PHONY: bin clean

bin: bin/danton

clean:
	@rm -rf bin

bin/%: src/%.c
	@mkdir -p bin
	@gcc -o $@ $(CFLAGS) $(INCLUDE) $< $(LIBS)
