# PTAPI manual

PTAPI is a C++ library (headers only) that runs Parallel Tempering on a
multi-core machine for any optimization problem the user describes with one
struct and one class.

Read [Getting started](getting-started.md) first; the other chapters are
independent of each other.

| Chapter | What it covers |
|---|---|
| [Getting started](getting-started.md) | requirements, how the method and the API work, a complete program explained line by line, how to choose the parameters |
| [The problem class](problem.md) | solution struct, `construction`, `neighbor`, the forms of `evaluate`, `neighborEval` and delta evaluation, random numbers, several moves, a heuristic for each replica |
| [Configuring the run](configuration.md) | the constructor, threads, seed, temperature ladder, initial solutions, exchange scheme, checkpoint |
| [Stop criteria](stopping.md) | rounds, time, target, lower bound, stagnation, a condition of the user, stopping from another thread |
| [Temperatures](temperatures.md) | the seven dynamic updates, the update interval, the two automatic calibrations, how to read the diagnostics of a ladder |
| [Results of a run](results.md) | best solution, times, evaluations, convergence, acceptance rates, time of each component |
| [Experiments](experiments.md) | R executions per instance, files and summaries, figures, several machines, statistical comparison of configurations |
| [Command line and irace](tuning.md) | `PTOptions`, tuning the parameters with irace |
| [Migrating from version 1](migrating-from-v1.md) | what compiles unchanged, what behaves differently, what changed inside, measured speed |
| [Reference](reference.md) | every public function on one page |
| [References](references.md) | the papers behind each method |

Other material:

- [`example/quickstart.cpp`](../example/quickstart.cpp) — the program of the first chapter.
- [`example/functions/`](../example/functions) — one short program for each group of functions; `make run` builds and runs them all.
- [`report/PTAPI-v2.pdf`](report/PTAPI-v2.pdf) — technical report of version 2, in Portuguese: every change, the source of each method and the measurements.
