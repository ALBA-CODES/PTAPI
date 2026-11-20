# Stop criteria

A run always ends when it completes the PTL rounds of the constructor. The
criteria below can be added and combined; the first one reached stops the run.
To stop only by one of them, give a very large PTL.

| Function | The run stops when |
|---|---|
| `setTimeLimit(double seconds)` | the time is over (fractions of a second are accepted) |
| `setTimeEnd(int seconds)` | the same, in whole seconds (version 1) |
| `setTargetValue(double value)` | the best value is `<= value` |
| `setLowerBound(double bound, double gap = 0)` | the best value is within `gap` (relative; 0.01 = 1%) of a known lower bound: `best <= bound + gap*|bound|` |
| `setStagnationRounds(int rounds)` | `rounds` rounds passed without a new best solution |
| `setStagnationTime(double seconds)` | `seconds` passed without a new best solution |
| `setPTLStopEnd(float fraction)` | `fraction x PTL` rounds passed without a new best solution (version 1) |
| `setStopCondition(f)` | `f(const PTStatus&)` returns `true` |
| `requestStop()` | it is called, from any thread or from the callback of `setOnImprove` |

`setTargetValue` and `setLowerBound` set the same target: the last one called
is the one in effect.

## A condition of your own

```cpp
pt.setStopCondition([](const PTStatus& st){
    return st.round >= 200 && st.bestValue < 7500.0;
});
```

The function is called at the end of every round, from a worker thread, one
call at a time. `PTStatus` has:

| Field | Meaning |
|---|---|
| `round` | round that has just ended |
| `elapsedMs` | time since `start()` |
| `bestValue` | best value so far |
| `bestRound`, `bestMs` | round and time in which it was found |

## Following the run

```cpp
pt.setOnImprove([&](const mySolution& s, int round, double ms){
    printf("%.1f ms, round %d: %.0f\n", ms, round, s.evalSol);
    if(goodEnough(s)) pt.requestStop();
});
```

The callback is called on every new best solution, from a worker thread, one
call at a time. Keep it short: the replica that found the solution waits for it.

## Why it stopped

```cpp
int why = pt.getStopReason();           // PTStopReason
const char* name = pt.getStopReasonName();
```

| Value | Name | Reason |
|---|---|---|
| `PT_STOP_PTL` | `PTL` | all the rounds were executed |
| `PT_STOP_TIME` | `time` | time limit |
| `PT_STOP_TARGET` | `target` | target value or lower bound reached |
| `PT_STOP_STAGNATION` | `stagnation` | rounds or time without a new best solution |
| `PT_STOP_REQUEST` | `request` | `requestStop()` |
| `PT_STOP_CONDITION` | `condition` | the function of `setStopCondition` |
| `PT_STOP_ERROR` | `error` | an exception in the problem code (thrown again by `start()`) |

`getPTLEnd()` gives the round in which the run ended and
`getBestPTLSolIndex()` the round of the best solution.

The two time criteria are checked by the thread that called `start()` and the
chains stop at their next move, so the run ends a few milliseconds after the
limit at most. The target is checked whenever a new best solution is found, and
the other criteria at the end of a round.

Example: [`exampleStop.cpp`](../example/functions/exampleStop.cpp).
