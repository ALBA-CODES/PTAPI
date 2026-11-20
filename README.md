# PTAPI — Parallel Tempering API

[![tests](https://github.com/ALBA-CODES/PTAPI/actions/workflows/ci.yml/badge.svg)](https://github.com/ALBA-CODES/PTAPI/actions/workflows/ci.yml)
[![license: CC BY-NC 4.0](https://img.shields.io/badge/license-CC%20BY--NC%204.0-blue.svg)](LICENSE)

A C++ library, headers only, that runs **Parallel Tempering** on multi-core
machines for any optimization problem. You write one struct (the solution) and
one class (how to build, change and evaluate a solution); the API runs the
replicas in parallel, exchanges them, adjusts the temperatures and reports the
results.

This is **version 2**. Programs written for version 1 compile without changes;
see [what is new](CHANGELOG.md) and
[Migrating from version 1](docs/migrating-from-v1.md).

## Quick start

```cpp
#include "include/PT.h"

// Step 1: the solution
struct solTSP: public solution{ std::vector<int> sol; };

// Step 2: the problem
class TSP: public Problem<solTSP>{
public:
    solTSP construction();                  // an initial solution, already evaluated
    solTSP neighbor(solTSP sol);            // one move
    double evaluate(const solTSP& sol);     // the cost (minimized)
};

// Step 3: the PT
int main(){
    TSP prob(/* instance */);
    //              Tmin  Tmax  replicas MCL  PTL  distribution update interval
    PT<solTSP> pt(1.0f, 60.0f, 16,     500, 1000, 4,          0,     50);
    pt.setSeed(42);
    pt.setTimeLimit(30.0);
    solTSP best = pt.start(0, &prob);       // 0 = all the cores
    std::cout << best.evalSol << "\n";
}
```

```
g++ main.cpp -std=c++17 -O3 -pthread -o main
```

A complete program is [`example/quickstart.cpp`](example/quickstart.cpp),
explained line by line in [Getting started](docs/getting-started.md).

## What it does

- **Parallel by construction.** The run is a dataflow graph executed by a
  thread pool; with a seed, the result does not depend on the number of threads.
- **Generic.** The core knows nothing about the problem. Three forms of
  `evaluate`, optional delta evaluation without copies (`neighborEval`),
  several moves with probabilities, a constructive heuristic for each replica.
- **Stop criteria.** Rounds, time, target value, lower bound with a gap,
  stagnation, a condition of your own.
- **Temperatures.** Four initial distributions, seven dynamic updates, two
  automatic calibrations, the even-odd (non-reversible) exchange scheme, and
  diagnostics of the ladder (exchange rates, round trips, communication
  barrier).
- **Long runs.** Checkpoint and resume.
- **Experiments.** R executions per instance with summaries in CSV, text and
  LaTeX, convergence and time-to-target files with ready figures, and the
  statistical comparison of configurations (Mann-Whitney, Wilcoxon, Friedman).
- **Tuning.** Every parameter on the command line and a ready scenario for
  [irace](https://mlopez-ibanez.github.io/irace/).

## Documentation

The [manual](docs/README.md):
[getting started](docs/getting-started.md) ·
[the problem class](docs/problem.md) ·
[configuring the run](docs/configuration.md) ·
[stop criteria](docs/stopping.md) ·
[temperatures](docs/temperatures.md) ·
[results](docs/results.md) ·
[experiments](docs/experiments.md) ·
[command line and irace](docs/tuning.md) ·
[migrating from version 1](docs/migrating-from-v1.md) ·
[reference](docs/reference.md)

## Repository

```
include/    the API (headers only); Experiment.h, Comparison.h and Options.h are optional
example/    quickstart.cpp, exampleTSP.cpp and the SSP, IPMTC and RCPMS examples
  functions/  one short example for each group of functions
instances/  instance files of the examples
irace/      tuning the parameters with irace
tests/      unit tests
bench/      speed of version 2 against version 1
v1/         the core of version 1, kept for the benchmark
docs/       the manual and the technical report of version 2 (in Portuguese)
```

```
make test       # unit tests
make check      # the same, also under ThreadSanitizer and AddressSanitizer
make examples   # builds and runs the examples
make bench      # version 2 against version 1 (on an idle machine)
```

Requirements: a C++17 compiler (GCC or Clang) and POSIX threads.

## Citing

If you use PTAPI in your research, please cite the thesis in which it was
proposed; the data are in [`CITATION.cff`](CITATION.cff) (the "Cite this
repository" button of GitHub).

## License

[Creative Commons Attribution-NonCommercial 4.0 International](LICENSE)
(CC BY-NC 4.0), the same license of version 1: PTAPI can be used and adapted
for non-commercial purposes, with proper credit given.
