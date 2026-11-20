# Getting started

## Requirements

- A C++17 compiler (GCC or Clang) and POSIX threads.
- Nothing to install: the API is the folder `include/`.

```
g++ main.cpp MyProblem.cpp -std=c++17 -O3 -pthread -o main
```

To check the copy you downloaded: `make test` (unit tests) and `make examples`
at the root of the repository.

## Parallel Tempering in one paragraph

Several copies of the search (the *replicas*) run at the same time, each one at
its own *temperature*. A replica at temperature `T` proposes a neighbor of its
solution and accepts it with probability `min(1, exp(-delta/T))`, where `delta`
is the change of the cost: cold replicas accept almost only improvements
(intensification) and hot ones accept almost everything (diversification).
After a fixed number of moves, replicas at adjacent temperatures propose to
*exchange* their solutions, with probability
`min(1, exp((1/T1 - 1/T2)(E1 - E2)))`. Good solutions found by the hot replicas
travel down to the cold ones, and solutions trapped in a local optimum travel
up and are released. The API minimizes: a smaller `evalSol` is better.

## How the API runs it

| Term | Meaning | Constructor argument |
|---|---|---|
| replica | one Markov chain at one temperature | 3rd (number of replicas) |
| MCL | moves of a replica between two exchange proposals | 4th |
| round | MCL moves of every replica followed by the exchange proposals | — |
| PTL | number of rounds | 5th |

A run of `R` replicas evaluates `R x MCL x PTL` neighbors.

The run is a *dataflow graph*: one node for each chain of MCL moves, one for
each exchange proposal, and the nodes of dynamic temperature updates. A node
runs as soon as the nodes it depends on have finished, on any thread of a pool,
so the replicas advance in parallel without a global barrier at each round.
Each replica has its own random engine: with a seed, the result does not
depend on the number of threads.

## The three steps

The complete program is [`example/quickstart.cpp`](../example/quickstart.cpp):
divide 60 numbers in two groups whose sums are as close as possible.

### Step 1: the solution

A struct that extends `solution`. The base struct brings `evalSol`, the cost
of the solution; add whatever encodes a solution of your problem.

```cpp
#include "include/PT.h"

struct solPart: public solution{
    std::vector<char> side;   // side[i] = 0 or 1: the group of number i
    double diff = 0;          // sum of group 0 minus sum of group 1
};
```

### Step 2: the problem

A class that extends `Problem<S>` and says how to build a solution, how to
change it, and how much it costs. Nothing else in the API depends on the
problem.

```cpp
class Partition: public Problem<solPart>{
    std::vector<double> num;
public:
    Partition(int n, unsigned instanceSeed);

    solPart construction(){                 // an initial solution, already evaluated
        solPart s;
        for(size_t i = 0; i < num.size(); i++) s.side.push_back((char)ptapi::randInt(0, 1));
        s.evalSol = evaluate(s);
        s.diff = difference(s);
        return s;
    }
    solPart neighbor(solPart s){            // one move: a number changes group
        int i = ptapi::randInt(0, (int)num.size()-1);
        s.side[i] = 1 - s.side[i];
        s.diff = difference(s);
        return s;
    }
    double evaluate(const solPart& s){      // the cost (minimized)
        return std::fabs(difference(s));
    }
    double difference(const solPart& s) const;
};
```

`ptapi::randInt` draws from the random engine of the replica that is running;
use it (or `ptapi::rng()`, `ptapi::rand01()`) instead of creating your own
engine. See [The problem class](problem.md) for the optional functions, among
them `neighborEval`, which makes a move without copying the solution and with
delta evaluation.

### Step 3: the PT

```cpp
int main(){
    Partition prob(60, 2026);

    // PT(Tmin, Tmax, replicas, MCL, PTL, distribution, updateType, updateInterval)
    PT<solPart> pt(10.0f, 1000000.0f, 16, 400, 2000, 4, 0, 20);
    pt.setSeed(42);                        // reproducible run
    pt.setTimeLimit(5.0);                  // at most 5 seconds

    solPart best = pt.start(0, &prob);     // 0 = all the cores

    printf("best difference: %.0f\n", best.evalSol);
}
```

Pass to `start()` a pointer of your own class (`Partition*`), not of
`Problem<solPart>`: this is how the API finds your `evaluate`.

Output of the program, which also prints some of the
[results of the run](results.md):

```
best difference: 1 (round 1411 of 2000, 104.0 ms of 149.9 ms)
12800000 evaluations on 16 threads, seed 42, stop: PTL
round trips: 12   barrier: 5.96
final ladder : 10 22 46 100 215 464 1000 2154 4642 10000 21544 46416 100000 215443 464159 1000000
exchange rate: 0.40 0.47 0.51 0.59 0.65 0.71 0.62 0.66 0.62 0.63 0.65 0.61 0.62 0.64 0.65
```

## Choosing the parameters

There is no setting that is right for every problem; this is a way to start.

1. **Temperatures.** `Tmin` should accept almost no worsening move and `Tmax`
   a good part of them. If you do not know the size of a typical worsening,
   let the API measure it: `pt.calibrateTemperatures(&prob)` before `start()`
   ([Temperatures](temperatures.md)). Use distribution 4 (geometric).
2. **Replicas.** As many as the cores you have is a good first value. After a
   run, look at `getSwapAcceptanceRates()`: rates near zero between two
   temperatures mean that solutions do not cross that gap (more replicas, or a
   narrower interval); rates above 50% everywhere mean replicas to spare.
   `getRoundTrips()` equal to zero means that no solution travelled from one
   end of the ladder to the other.
3. **MCL.** Enough moves for a replica to change its solution between two
   exchanges; a few times the size of the solution is a common start. Very
   small values make the exchange nodes and the threads waiting for them a
   large share of the time (`setProfiling(true)` shows it).
4. **PTL and the stop.** Give a large PTL and stop by time, by a target or by
   stagnation ([Stop criteria](stopping.md)).
5. **Update type.** Start with 0 (fixed ladder). A dynamic update helps when
   the ladder is badly placed; it did not beat a well-placed fixed ladder in
   the problems measured so far.
6. To search the parameters automatically, see
   [Command line and irace](tuning.md).

## Folders

```
include/    the API (headers only); Experiment.h, Comparison.h and Options.h are optional
example/    quickstart.cpp, exampleTSP.cpp and the SSP, IPMTC and RCPMS examples
  functions/  one short example for each group of functions (make run)
instances/  instance files of the examples
irace/      tuning the parameters with irace
tests/      unit tests (make test | make tsan | make asan)
bench/      speed of version 2 against version 1 (make bench)
v1/         the core of version 1, kept for the benchmark
docs/       this manual
```
