# Contributing

Bug reports and pull requests are welcome.

Three rules are checked for every change in `include/`:

1. **The API stays generic.** The core knows nothing about the problem; only the
   class derived from `Problem<S>` does.
2. **The tests pass in the three modes.** In `tests/`: `make test`, `make tsan`
   (ThreadSanitizer) and `make asan` (AddressSanitizer + UBSan). A new function
   comes with its test in `tests/test_ptapi.cpp` and, if it is part of the
   public interface, with an example in `example/functions/`.
3. **Version 2 is not slower than version 1.** Run `make bench` in `bench/`
   before and after the change, on an idle machine, and put the numbers in the
   pull request.

`make run` in `example/functions/` builds and runs every example; each one
returns 0 only if it behaved as its comments say.
