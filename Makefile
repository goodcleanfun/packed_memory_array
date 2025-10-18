
install:
	clib install --dev

test:
	@$(CC) test.c $(CFLAGS) -g -I src -I deps -I deps/greatest $(LDFLAGS) -o $@
	@./$@

debug:
	@$(CC) test.c $(CFLAGS) -g -I src -I deps -I deps/greatest $(LDFLAGS) -o $@
	gdb @./$@

.PHONY: install test
