# Configuring the run

## The constructor

```cpp
PT<S> pt(tempMin, tempMax, replicas, MCL, PTL, distribution, updateType, updateInterval);
```

| # | Argument | Meaning |
|---|---|---|
| 1 | `float tempMin` | first (lowest) temperature |
| 2 | `float tempMax` | last (highest) temperature |
| 3 | `int replicas` | number of replicas = number of temperatures (1 is accepted: plain Metropolis) |
| 4 | `int MCL` | moves of each replica between two exchange proposals |
| 5 | `int PTL` | number of rounds |
| 6 | `int distribution` | initial distribution of the temperatures, 1 to 4 (below) |
| 7 | `int updateType` | dynamic update of the temperatures, 0 to 6 ([Temperatures](temperatures.md)); 0 = fixed |
| 8 | `int updateInterval` | exchange proposals between two updates; `PT_UPDATE_AUTO` (0) lets the API choose |

Invalid values throw `std::invalid_argument`.

| `distribution` | Temperatures `T[i]`, `i = 0 .. n-1`, `x = i/(n-1)` |
|---|---|
| 1 linear | `Tmin + (Tmax - Tmin) x` |
| 2 inverse linear | `1/T[i]` linear between `1/Tmax` and `1/Tmin`: the ladder goes **from `Tmax` down to `Tmin`** |
| 3 exponential | `log T[i]` linear between `log Tmin` and `log Tmax` |
| 4 geometric | `Tmin (Tmax/Tmin)^x` |

Distributions 3 and 4 are two ways of writing the same ladder (constant ratio
between adjacent temperatures); they are both kept because version 1 has them.
Update types 1 and 5 need an ascending ladder (distributions 1, 3 and 4).

## Running

```cpp
S best = pt.start(threads, &prob);
```

- `threads <= 0` uses all the cores. The API never creates more threads than
  replicas; `getThreadsUsed()` tells how many were created.
- `&prob` must be a pointer of the user's class, not of `Problem<S>`.
- `start()` can be called again on the same object; the functions below are
  called before it.
- The best solution is also available afterwards with `getBestSol()`.

## Seed

```cpp
pt.setSeed(42);
```

Every replica and every exchange node gets its own random engine, derived from
the seed. Without `setSeed` a seed is drawn; `getSeed()` returns the one that
was used, so that any run can be repeated.

A run is reproducible, for any number of threads, when the problem class only
uses the random functions of the API, the run stops by PTL, and the update
type is 0, 1, 2 or 5 (type 1 with an update interval of at least 2). Stops that
depend on the clock, and update types 3, 4 and 6 (which read statistics of
other temperatures while they run), are not reproducible.
Example: [`exampleSeed.cpp`](../example/functions/exampleSeed.cpp).

## A ladder of your own

```cpp
pt.setTemperatures({1.0, 1.5, 2.5, 4.0, 7.0, 12.0, 20.0, 35.0});
```

Replaces the ladder of the constructor; the number of replicas becomes the size
of the list. `getInitialTemperatures()` returns the ladder a run starts from
(also after a calibration). Example:
[`exampleTemperatures.cpp`](../example/functions/exampleTemperatures.cpp).

## Initial solutions

```cpp
pt.setInitialSolutions(solutions);     // std::vector<S>
```

Warm start: replica `i` starts from `solutions[i % solutions.size()]` instead
of calling `construction()`. The API evaluates them, so the class needs
`evaluate`. To *build* a different solution for each replica, override
`constructionReplica` ([The problem class](problem.md)); to build them in
parallel, `pt.setParallelConstruction(true)`.

## Exchange scheme

```cpp
pt.setSwapScheme(PT_SWAP_EVEN_ODD);    // default: PT_SWAP_SEQUENTIAL
```

| Scheme | Pairs that propose an exchange |
|---|---|
| `PT_SWAP_SEQUENTIAL` | every pair on every round, from the first temperature up (version 1) |
| `PT_SWAP_EVEN_ODD` | pairs (0,1), (2,3), ... on odd rounds and (1,2), (3,4), ... on even rounds: the non-reversible scheme of Okabe et al. (2001) and Syed et al. (2019) |

With the even-odd scheme each pair proposes on half of the rounds. Example:
[`exampleSwapScheme.cpp`](../example/functions/exampleSwapScheme.cpp).

## Checkpoint and resume

The state of a run (temperatures, solution of every replica, best solution,
rounds and time already used) can be saved and continued later. The API writes
the objective values; the user gives the two functions that write and read the
rest of one solution.

```cpp
auto write = [](std::ostream& out, const solTSP& s){ for(int c : s.sol) out << c << ' '; };
auto read  = [n](std::istream& in, solTSP& s){ s.sol.resize(n); for(int& c : s.sol) in >> c; };

pt.loadCheckpoint("run.ckpt", read);        // continues, if the file exists (returns false otherwise)
pt.setCheckpoint("run.ckpt", write, 60.0);  // saves every 60 s and at the end of start()
solTSP best = pt.start(0, &prob);
```

`saveCheckpoint(file, writer)` writes the state of the last `start()` on demand.

A resumed execution runs the remaining PTL and the remaining time. It is a
continuation of the search, not a bit-identical replay: the random engines and
the statistics of the dynamic updates start again, and `getConvergence()` only
lists the improvements found after the resume. The periodic copy is taken by
each replica at the end of its round, so it can be up to one interval old.
Example: [`exampleCheckpoint.cpp`](../example/functions/exampleCheckpoint.cpp).

## Compile-time options

| Macro | Effect |
|---|---|
| `-DPTAPI_SPIN=n` | iterations an idle thread spins before it sleeps (default 5000) |
| `-DPTAPI_NO_USING_NAMESPACE_STD` | `PT.h` does not export `using namespace std;` (it does by default, because programs written for version 1 rely on it) |
