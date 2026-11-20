# Temperatures

The temperature ladder decides how well Parallel Tempering works: solutions
must be able to travel between the two ends of it. This chapter covers the
tools the API has to place the temperatures: the dynamic updates during the
run, the automatic calibrations before it, and the diagnostics to judge a
ladder.

## Diagnostics of a ladder

| Function | What it tells |
|---|---|
| `getSwapAcceptanceRates()` | exchange rate of each adjacent pair. Near zero: solutions do not cross that gap. Near one: the two temperatures are doing the same work |
| `getMoveAcceptanceRates()` | fraction of the moves each temperature accepted. Near zero at the first temperature and well above zero at the last is the usual aim |
| `getRoundTrips()` | replicas that travelled from the first temperature to the last and back. Zero means the two ends never communicated |
| `getCommunicationBarrier()` | sum of the exchange rejection rates (Syed et al., 2019): about half of the number of replicas that the temperature interval needs |
| `getTemperatures()` | the ladder at the end of the run, after the dynamic updates |

Example: [`exampleSwapScheme.cpp`](../example/functions/exampleSwapScheme.cpp).

## Dynamic updates

The 7th argument of the constructor chooses a rule that moves the temperatures
during the run; the 8th is the number of exchange proposals between two
updates.

| Type | Rule | First and last temperature | Reproducible with a seed |
|---|---|---|---|
| 0 | fixed ladder | — | yes |
| 1 | exchange rate of 23% in every pair | fixed | yes (interval >= 2) |
| 2 | equal exchange rates (Hukushima, 1999) | fixed | yes |
| 3 | feedback-optimized (Katzgraber et al., 2006), as in the thesis | fixed | no |
| 4 | feedback-optimized, by quantiles of the density | fixed | no |
| 5 | stochastic approximation towards 23.4% (Miasojedow, Moulines and Vihola, 2013) | the first is fixed; the last is a result | yes |
| 6 | equal rejection rates (Syed et al., 2019) | fixed | no |

- **Type 1.** Each pair changes its gap towards the one that gives 23%
  (narrower when its rate is below 20%, wider above 26%); the two temperatures
  of the pair share the change and the ends never move. Where 23% cannot be
  reached by every pair (too few or too many replicas for the interval), the
  deviation is spread. Needs an ascending ladder.
- **Type 2.** Each temperature moves to the point that equalizes the exchange
  probability with its two neighbours, computed from the mean energies. A new
  temperature is only accepted between its two neighbours.
- **Type 3.** The temperatures follow the density that maximizes the flow of
  replicas between the two ends, measured from the fraction of replicas that
  last visited the first temperature. Implemented as described in the thesis.
- **Type 4.** The same idea computed for the whole ladder at once: the new
  temperatures divide the integral of the optimal density in equal parts. The
  ends never move and the ladder stays ordered; half of the step is applied at
  each update.
- **Type 5.** The ladder is written as `T[i+1] = T[i]*exp(exp(rho[i]))` and, at
  every update, each pair does
  `rho[i] += (n+1)^-0.6 * (average exchange probability - 0.234)`, where `n`
  counts the updates. The last temperature is limited to the highest
  temperature of the initial ladder. Needs an ascending ladder.
- **Type 6.** The cumulative sum of the average rejection probabilities of the
  pairs is interpolated in `1/T` and the temperatures are placed at equal steps
  of it, so that every pair rejects the same fraction of the exchanges. Half of
  the step is applied at each update.

In the problems measured so far (TSP, SSP, IPMTC, RCPMS) no dynamic rule beat
a well-placed fixed ladder consistently; their value is in not needing that
ladder to begin with. Example:
[`exampleTempUpdate.cpp`](../example/functions/exampleTempUpdate.cpp).

### Update interval

The 8th argument of the constructor is the number of exchange proposals
between two updates. With `PT_UPDATE_AUTO` (0) the API chooses it: the first
update happens after 10 proposals and each interval is twice the previous one
(10, 20, 40, ...), as in the tuning rounds of Syed et al. (2019). The ladder is
corrected often at the beginning and with statistics of longer windows later,
and there is one parameter less to tune.

```cpp
PT<S> pt(1.0f, 60.0f, 16, 500, 1000, 4, 5, PT_UPDATE_AUTO);
pt.setAutoUpdate(20, 1.5);     // optional: first interval and growth (defaults 10 and 2)
```

Measured on the TSP (20 seeds, 16 replicas, 1000 rounds), the automatic
interval was equivalent to a fixed interval of 50 within the noise, with a
tendency to be slightly worse. It removes a parameter; it is not a way of
getting better solutions. Example:
[`exampleAutoUpdate.cpp`](../example/functions/exampleAutoUpdate.cpp).

## Calibration before the run

### First and last temperature: `calibrateTemperatures`

```cpp
std::pair<double,double> t = pt.calibrateTemperatures(&prob);   // {Tmin, Tmax}
// all the arguments: (prob, samples = 5000, acceptHot = 0.02, acceptCold = 0.001)
```

Runs a descent of `samples` moves from a constructed solution and keeps the
worsening moves of its second half, when the solution is close to a local
optimum. With `A(T)` = average of `exp(-delta/T)` over these moves (the
fraction of them that a replica at temperature `T` accepts), it solves
`A(Tmin) = acceptCold` and `A(Tmax) = acceptHot` and rebuilds the ladder with
the distribution and the number of replicas of the constructor.

It costs `samples` evaluations, which are not counted by `getEvaluations()`.
The defaults are a starting point that worked on the problems of `example/`;
the best temperatures still depend on the problem and on the budget. For a
target on the acceptance of the moves at the last temperature, raise
`acceptHot` (for example 0.23). Example:
[`exampleCalibrateTemperatures.cpp`](../example/functions/exampleCalibrateTemperatures.cpp).

### The whole ladder for an exchange rate: `calibrateLadder`

```cpp
std::vector<double> ladder = pt.calibrateLadder(&prob);
// all the arguments: (prob, threads = 0, swapRate = 0.23, rounds = 200, iterations = 10, tMin = 0)
```

Keeps the first temperature (`tMin`, or the one already defined) and searches
all the others with pilot runs of the PT itself, so that every adjacent pair
exchanges at `swapRate`. The last temperature is a result: where the given
number of replicas reaches with that rate.

The pilot runs use the stochastic approximation of update type 5 aimed at
`swapRate`: every 10 exchange proposals each pair corrects its gap with a step
that decreases along the whole calibration. `iterations` pilot runs of `rounds`
rounds are made, each one continuing from the replicas and the ladder of the
previous one, and the result is the mean over the second half of them. If
there are more replicas than the rate needs, the ladder stops at the first
temperature that accepts 90% of the moves (a random walk) and the remaining
replicas stay close to it: a sign that fewer replicas are enough.

The pilot runs cost `iterations x rounds x replicas x MCL` evaluations, which
are not part of the run that follows, and ignore the stop criteria, the
callback and the checkpoint of the object.

**How close it gets** (TSP of 100 cities, 10 replicas, 10 seeds, rates measured
in a run of 5000 rounds after the calibration, target 23%):

| First temperature | Mean deviation of a pair from the target | Variation of the rate between seeds with one fixed ladder |
|---|---|---|
| 20 (the replicas keep moving) | 0.025 (0.021 with 30 pilot runs of 1000 rounds) | 0.02 |
| 1 (the cold replicas are trapped in local optima) | 0.09 | 0.10 |

Where the replicas keep moving the result is as close as the rate can be
measured. Where they are trapped, the exchange rate of a pair depends on the
two solutions it holds, not on the ladder, and rises as the search freezes: no
ladder meets the target in every run. In that case start from a higher first
temperature or use a dynamic update during the run.

The two calibrations can be combined: `calibrateTemperatures` first (it defines
the first temperature), then `calibrateLadder`. Example:
[`exampleCalibrateLadder.cpp`](../example/functions/exampleCalibrateLadder.cpp).
