# Command line and irace

The optional header `Options.h` reads the parameters of the PT from the
command line, which is what an automatic tuner needs.

## `PTOptions`

`PTOptions` reads every parameter of the PT in the form `--NAME value`. The
first argument that is not an option is the instance (`opt.instance`). Options
it does not know are kept for the problem (`opt.get("--ALPHA", 0.5)`,
`opt.has("--ALPHA")`), so the parameters of the problem are tuned in the same
way.

```cpp
#include "include/Options.h"

int main(int argc, char* argv[]){
    PTOptions opt(argc, argv);                    // opt.instance = first argument
    MyProblem prob(opt.instance);
    mySolution best = opt.run<mySolution>(&prob); // builds the PT, calibrates if asked, runs
    if(!opt.printForTuner()) std::cout << best.evalSol << "\n";
}
```

| Option | Meaning | Default |
|---|---|---|
| `--TEMP_INIT`, `--TEMP_FIM` | first and last temperature | 1, 10 |
| `--N_REPLICAS` | number of replicas | 8 |
| `--MCL`, `--PTL` | moves between exchanges, rounds | 100, 1000 |
| `--TEMP_DIST` | temperature distribution, 1 to 4 | 4 |
| `--TYPE_UPDATE` | dynamic update, 0 to 6 | 0 |
| `--TEMP_UPDATE` | proposals between updates; 0 = automatic | 50 |
| `--THREAD_USED` | threads; 0 = all the cores | 0 |
| `--SEED` | seed | drawn |
| `--TIME_LIMIT` | seconds | none |
| `--TARGET` | stop when the best value is `<=` it | none |
| `--STAGNATION_ROUNDS`, `--STAGNATION_TIME`, `--STAGNATION_PTL` | stop without a new best: rounds, seconds, fraction of PTL | off |
| `--SWAP_SCHEME` | 0 sequential, 1 even-odd | 0 |
| `--CALIBRATE` | 1 `calibrateTemperatures`, 2 `calibrateLadder`, 3 both | 0 |
| `--PARALLEL_CONSTRUCTION`, `--PROFILE` | 0 or 1 | 0 |
| `--IRACE`, `--IRACE_COST`, `--IRACE_TIME_WEIGHT`, `--IRACE_PENALTY` | output for the tuner (below) | off |

The names are those of the examples published with version 1, plus the
additions of version 2. `opt.make<S>()` returns the configured `PT` object
without running it, `opt.apply(pt)` configures an object you built, and after
`opt.run()` the results are in `opt.result` (a `PTRunResult`, see
[Experiments](experiments.md)). Example:
[`exampleOptions.cpp`](../example/functions/exampleOptions.cpp).

## Tuning with irace

With `--IRACE 1` the program prints only the cost of the run, which is what
[irace](https://mlopez-ibanez.github.io/irace/) reads (`--IRACE 2` prints
`cost time`). `--IRACE_COST` chooses what is minimized:

| `--IRACE_COST` | Cost |
|---|---|
| `value` (default) | the objective value at the end of the run; give every configuration the same `--TIME_LIMIT`, so that quality is compared under the same time |
| `time` | seconds to reach `--TARGET`; a run that does not reach it within `--TIME_LIMIT` costs `--IRACE_PENALTY` (default 10) times the limit |
| `value+time` | `value + w x seconds`, with `w` given by `--IRACE_TIME_WEIGHT`: a price for each second, for runs that stop by themselves (stagnation, PTL) |

The folder `irace/` has a complete scenario for the SSP example:
`parameters.txt` (the parameters and their ranges), `scenario.txt`,
`target-runner` (the script irace calls; the budget of each run is set there)
and `ptSSP.cpp` (20 lines, of which two depend on the problem). `make check`
calls the target runner once as irace does and `make tune` runs irace.

To tune your own problem, copy the folder and change: the two lines of
`ptSSP.cpp` that name the problem class, the instances folder in
`scenario.txt`, the ranges in `parameters.txt`, and the budget in
`target-runner`. irace is an R package: `install.packages("irace")`.

Give every configuration the same budget (the line `FIXED` of `target-runner`);
otherwise irace simply chooses the configuration that runs longer. The scenario
was run end to end with irace 4.5.1 (300 experiments on 10 SSP instances) to
check the integration; its result is not a recommendation of parameters.
