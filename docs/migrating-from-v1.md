# Migrating from version 1

## What stays the same

Programs written for version 1 compile with the headers of version 2 without
changes: the constructor of `PT`, `start()`, `getBestSol()`,
`getBestPTLSolIndex()`, `getPTLEnd()`, `setTimeEnd()` and `setPTLStopEnd()` are
unchanged, and so are the three steps to use the API. Replace the folder
`include/` and compile with `-std=c++17` (or newer).

To also reproduce the *numbers* of version 1, call
`pt.setV1Compatibility(true)`: it restores the exchange test with truncated
energies, the exchange rate of update type 1 and the rules of update types 1
and 2 described below. Example:
[`exampleV1Compatibility.cpp`](../example/functions/exampleV1Compatibility.cpp).

## Worth changing in old code

| In version 1 | In version 2 | Why |
|---|---|---|
| `std::random_device` and an engine created inside `neighbor()` | `ptapi::rng()`, `ptapi::randInt()`, `ptapi::rand01()` | faster, and the run becomes reproducible with `setSeed` |
| `double evaluate(S sol)` | `double evaluate(const S& sol)` | no copy of the solution |
| copying `Nup` and `Ndown` in the solution | nothing | the API stores the round-trip labels itself |
| `neighbor` + `evaluate` | `neighborEval` with delta evaluation | no copies in a Markov chain step |
| `setTimeEnd(int)` | `setTimeLimit(double)` | fractions of a second |

## What changed inside

| Topic | v1 | v2 |
|---|---|---|
| Idle threads | busy wait (100% CPU), including the thread that called `start()` | spin for a few microseconds, then sleep |
| Number of threads | always `thN` | `min(thN, replicas)`; `thN <= 0` means all cores |
| Node dependencies | mutex + flag vector scanned on every notification | one atomic counter per node |
| Queues | one mutex-protected queue per thread, scanned in a loop | one queue; a worker runs the next ready node itself |
| Copies per MC step | 3–4 deep copies of the solution | 1 (0 with `neighborEval`) |
| Exchange of replicas | 3–4 deep copies | `std::swap` |
| Best solution | collected from the nodes when they finish | global, updated when it is found |
| Random numbers | every node seeded with the same default seed | one engine per node, derived from one seed |
| Memory | nodes and pool never freed | freed at the end of `start()` |

## Behaviour that is different from v1

- **Exchange criterion with real-valued objectives.** v1 converted the energy
  of a replica to `int` before the exchange test (and tracked the best value as
  `int`). v2 uses `double`. Nothing changes for integer objectives.
- **Exchange rate used by update type 1 (23%).** v1 divided the exchanges
  accepted since the last update by the total number of rounds; v2 divides by
  the size of the window, as described in the thesis.
- **Dynamic temperature updates.** Types 1 and 2 now keep the first and the
  last temperature and only redistribute the others between them. In version 1
  the rule of type 1 could only bring temperatures closer, so the ladder
  drifted down and the maximum temperature was lost; type 2 accepted any
  positive value, which could leave the ladder out of order. Types 4, 5 and 6
  are new ([Temperatures](temperatures.md)).
- `setV1Compatibility(true)` restores the v1 behaviours above (exchange test,
  exchange rate, and the rules of types 1 and 2), to reproduce published results.
- **Exceptions** thrown by the problem code inside a worker thread stop the
  run and are rethrown by `start()` (v1 terminated the program).
- **Round-trip labels** (update types 3 and 4, and `getRoundTrips()`) are stored by the API. `Nup`/`Ndown`
  remain in `solution` only so that old code compiles; they are ignored.
- **Time limit** is enforced by the thread that called `start()`; MC nodes stop
  at the next step instead of after one or two more rounds.
- **Invalid parameters** (temperature distribution, update type, zero replicas)
  throw `std::invalid_argument`. One replica is accepted (plain Metropolis).
- `start()` can be called again on the same object.
- The API headers can be included from more than one `.cpp` file.
- `using namespace std;` is still exported by `PT.h` for compatibility; define
  `PTAPI_NO_USING_NAMESPACE_STD` before the include to turn it off.

The examples published with v1 (`double evaluate(S)`) compile and run without
changes. The copies in `example/` are written in the v2 style: `ptapi::rng()`,
`double evaluate(const S&)`, no `Nup`/`Ndown`, the move applied to the argument
of `neighbor()` without an extra copy, and the options `--SEED`, `--TIME_LIMIT`,
`--CALIBRATE` and `--SWAP_SCHEME` in the main programs, which also print the
statistics of the run.

Temperature distribution 2 (inverse linear) goes from the maximum to the
minimum temperature, as in v1; the other three are ascending.


## Measured speed

`make bench` in `bench/` builds the same program with the two versions and runs
them interleaved (TSP of 100 cities, cheap evaluation, so the cost of the API
itself shows). Wall time in seconds on a 24-thread machine:

| Replicas, MCL, PTL, threads | Version 1 | Version 2 |
|---|---|---|
| 8, 500, 2000, 8 | 0.70 | 0.61 |
| 8, 500, 2000, 24 | 0.81 | 0.60 |
| 30, 500, 500, 24 | 0.32 | 0.27 |
| 16, 20, 20000, 16 | 0.52 | 0.48 |

Version 2 also uses less CPU time for the same wall time, because idle threads
sleep instead of spinning. These numbers are for the classic
`neighbor` + `evaluate`; a problem that overrides `neighborEval` gains more.
Run the benchmark on an idle machine: on a loaded one the comparison is not
meaningful.

## Tests

Run `make test` in `tests/` after any change in the core, and `make tsan` /
`make asan` before a release. The tests use three small problems with
different encodings (permutation, binary, continuous) and cover: the four
temperature distributions, every stop criterion, reproducibility with a seed
for any number of threads, the equivalence of all the forms of `evaluate`, the
six dynamic temperature updates (with a fixed and with an automatic interval),
the two automatic calibrations of the temperatures, moves with probabilities,
checkpoint and resume, the construction of each replica (also in parallel), the
time of each component, the experiment manager, the command-line options, the
statistics, error handling, and a statistical check that each temperature
samples the Boltzmann distribution.

Run `make bench` in `bench/` to compare the speed with the previous version
(same problem, same parameters, both versions built from the same file).

