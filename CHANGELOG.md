# Changelog

## 2.0.0

Version 2 keeps the dataflow model and the interface of version 1: programs
written for version 1 compile without changes. The details, with the
measurements, are in [docs/migrating-from-v1.md](docs/migrating-from-v1.md).

**Faster core.** Threads that sleep instead of spinning, one atomic counter per
node instead of a mutex, one queue, no copies of the solution in a Markov chain
step (`neighborEval`), exchanges by `std::swap`.

**Corrected.** Random engines (every node used the same seed), the exchange
test with real-valued objectives (the energy was truncated to `int`), the
exchange rate of update type 1, ladders that drifted or went out of order with
update types 1 and 2, memory that was never freed, exceptions in worker threads.

**New.**
- Three forms of `evaluate`, `neighborEval` with delta evaluation, several moves
  with probabilities, a constructive heuristic for each replica, construction in
  parallel.
- Reproducible runs (`setSeed`), time limit with sub-second resolution, target
  value, lower bound with a gap, stagnation by rounds or by time, a stop
  condition of the user, the reason why the run ended.
- User-defined ladder, initial solutions, checkpoint and resume.
- Three more dynamic temperature updates (feedback by quantiles, stochastic
  approximation, equal rejection rates), an automatic update interval, and two
  automatic calibrations of the temperatures.
- Even-odd (non-reversible) exchange scheme, round trips, communication barrier.
- Results of the run: evaluations, times, convergence curve, acceptance rates,
  time of each component.
- Optional headers: `Experiment.h` (R executions per instance, summaries,
  convergence and time-to-target files, figures), `Comparison.h` (statistical
  comparison of configurations) and `Options.h` (command line, irace).
- Unit tests (normal, ThreadSanitizer, AddressSanitizer), a benchmark against
  version 1, and one example for each group of functions.

## 1.0.0

The version published with the thesis: `PT`, `Problem`, four temperature
distributions, dynamic updates 1 to 3, and the SSP, IPMTC and RCPMS examples.
