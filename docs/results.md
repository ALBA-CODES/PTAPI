# Results of a run

`start()` returns the best solution. Everything else about the last run is
read from the `PT` object afterwards.

| Function | Returns |
|---|---|
| `getBestSol()` | the best solution |
| `getBestPTLSolIndex()` | round in which it was found |
| `getPTLEnd()` | round in which the run ended |
| `getStopReason()`, `getStopReasonName()` | why it ended ([Stop criteria](stopping.md)) |
| `getElapsedMs()`, `getTimeToBestMs()` | total time and time to the best solution |
| `getEvaluations()` | number of neighbors evaluated |
| `getThreadsUsed()` | threads actually created |
| `getSeed()` | seed that was used |
| `getConvergence()` | every new best solution: `PTImprovement{timeMs, ptl, value}` |
| `getReplicas()` | final solution of each temperature |
| `getInitialTemperatures()`, `getTemperatures()` | the ladder at the beginning and at the end of the run |
| `getMoveAcceptanceRates()` | fraction of accepted moves, one per temperature |
| `getSwapAcceptanceRates()` | fraction of accepted exchanges, one per adjacent pair |
| `getRoundTrips()`, `getCommunicationBarrier()` | diagnostics of the ladder ([Temperatures](temperatures.md)) |
| `getMoveStatistics()` | `PTMoveStat{proposed, accepted}` for each move registered with `setMoveProbabilities` |
| `getTiming()` | where the time went (below) |

Example: [`exampleStatistics.cpp`](../example/functions/exampleStatistics.cpp).

## Convergence curve

```cpp
for(const PTImprovement& p : pt.getConvergence())
    printf("%.2f ms  round %d  %.0f\n", p.timeMs, p.ptl, p.value);
```

The list has one entry for each new best solution, in the order they were
found. The replicas are not all in the same round at the same instant, so the
rounds of the list are not always ascending; the times are. For the curve of
an instance over several runs, see [Experiments](experiments.md).

## Time of each component

```cpp
pt.setProfiling(true);
pt.start(0, &prob);
const PTTiming& t = pt.getTiming();
```

| Field | Meaning | Measured |
|---|---|---|
| `wallMs` | the whole `start()` | always |
| `constructionMs` | initial solutions (wall time) | always |
| `chainsMs` | Markov chains: the problem code (`neighbor`/`evaluate`) and the acceptance test | with profiling |
| `exchangesMs` | exchange proposals | with profiling |
| `updatesMs` | dynamic temperature updates | with profiling |
| `waitingMs` | threads waiting for a node to become ready | with profiling |
| `threads` | threads of the run | always |

The last four are summed over the threads, so they add up to
`threads x (wallMs - constructionMs)`. Measuring costs two clock readings per
node execution (not per move) and is off by default.

A large share of waiting means that the chains are too short for the number of
threads (raise MCL) or that there are more threads than the problem can use. A
share of the chains close to 100% means the time is in your `neighbor` and
`evaluate`: that is where to optimize (see `neighborEval` in
[The problem class](problem.md)). Example:
[`exampleTiming.cpp`](../example/functions/exampleTiming.cpp).
