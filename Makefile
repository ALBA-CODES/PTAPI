# PTAPI - shortcuts (each folder has its own Makefile)
#   make test       unit tests
#   make check      unit tests in the three modes (normal, ThreadSanitizer, AddressSanitizer)
#   make examples   builds every example and runs the quick start and one example per group of functions
#   make bench      speed of version 2 against version 1 (run it on an idle machine)
#   make clean
test:
	$(MAKE) -C tests test
check:
	$(MAKE) -C tests all
examples:
	$(MAKE) -C example all
	./example/quickstart
	$(MAKE) -C example/functions run
bench:
	$(MAKE) -C bench bench
clean:
	$(MAKE) -C example clean
	$(MAKE) -C example/functions clean
	$(MAKE) -C irace clean
	rm -rf bench/build tests/test_ptapi tests/test_ptapi_tsan tests/test_ptapi_asan

.PHONY: test check examples bench clean
