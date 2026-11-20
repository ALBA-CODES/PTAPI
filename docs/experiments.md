# Experiments

Two optional headers: `Experiment.h` runs the experiment of a paper and
`Comparison.h` compares several of them with statistical tests. They are not
included by `PT.h`.

## Running an experiment (`Experiment.h`)

The experiment of a paper — R executions for each instance, one file per
execution, and the best, mean, median and standard deviation of each instance
— is written once in the API. The user gives the function that runs one
execution:

```cpp
#include "include/Experiment.h"

PTExperiment exp("results");           // output folder
exp.setRuns(30);                       // executions per instance
exp.setFirstSeed(1);                   // execution k uses the seed 1 + k
exp.addInstance("instances/a.txt");    // what the runner receives
exp.setReference(1234.0);              // optional: best known value / lower bound of that instance
exp.addInstancesFromFile("list.txt");  // or one "path [reference]" per line

exp.run([](const std::string& instance, unsigned long long seed){
    MyProblem prob(instance);
    PT<mySolution> algo(1.0, 60.0, 16, 500, 1000, 4, 0, 50);
    algo.setSeed(seed);
    mySolution best = algo.start(0, &prob);
    return ptResult(algo, best);       // PTRunResult: everything the API measured
});
```

Files written under the output folder:

| File | Content |
|---|---|
| `<instance>/run_<k>.txt` | one execution: seed, value, times, rounds, evaluations, stop reason, final ladder, acceptance rates, convergence curve and, if the user fills `PTRunResult::solution`, the solution |
| `runs.csv` | one line per execution |
| `summary.csv` | one line per instance: runs, best, mean, median, standard deviation (sample), worst, mean time, mean time to best, mean evaluations and, with a reference, the gaps and the number of executions that reached it |
| `summary.txt` | the same table, aligned |
| `summary.tex` | the same table for LaTeX (`booktabs`) |
| `<instance>/convergence.csv` | convergence of the instance over its executions, by round: for each point of a grid (default 100 points, closer at the beginning; `setCurvePoints`), the best, the median, the mean and the worst of the values the executions had by then |
| `<instance>/convergence_time.csv` | the same by time (seconds) |
| `<instance>/time_to_target.csv` | instances with a reference: the time each execution took to reach it, in ascending order, with the empirical probability `(i - 0.5)/R` of having reached it by then (the data of a time-to-target plot, Aiex, Resende and Ribeiro, 2007) |

The three last files are plain CSV made to be plotted, and the folder gets two
scripts that draw them: `figures.gp` (`gnuplot figures.gp` writes
`<instance>/convergence.pdf` and `<instance>/time_to_target.pdf`) and
`figures.tex` (`lualatex figures.tex` writes `figures.pdf`, one figure per
page, with pgfplots; each `tikzpicture` can be copied into a paper). The API
does not run them. The convergence figure has the median, the mean and the
band from the best to the worst execution, with the rounds in logarithmic
scale. The curves come from the list of new best solutions of
each execution (`getConvergence`), so the first points can change a little
between two executions with the same seed: the order in which the replicas
report their first solutions depends on the threads. `setReference` also
accepts a target value instead of the best known one; `summary.csv` then has
the number of executions that reached it and the mean time they took.

The columns of `runs.csv` can be chosen. The default is `instance`, `run`,
`seed`, `value`, `seconds`, `seconds_to_best`, `rounds` (the PTL at which the
execution ended), `best_round`, `evaluations`, `threads`, `stop`,
`round_trips`. `addCsvColumns({...})` adds columns and `setCsvColumns({...})`
replaces the list; the names the API can fill by itself are those of
`PTExperiment::csvColumnNames()`: also `path`, `evaluations_per_second`,
`improvements`, `replicas`, `temp_first`, `temp_last`, `temperatures`, `temperatures_initial`,
`moves_proposed`, `moves_accepted`, `moves_rate` (one value per registered
move), `seconds_to_reference`, `swap_rate_mean`, `swap_rate_min`, `swap_rate_max`, `barrier`, `swap_rates`,
`move_rates`, `construction_seconds`, `chains_ms`, `exchanges_ms`,
`updates_ms`, `waiting_ms`, `reference`, `gap_percent` and `solution`. Values
of the user go in with `result.add("name", value)` inside the runner: each
name becomes one more column, and the value is kept in the file of the
execution.

An interrupted experiment continues where it stopped: executions whose file
exists are read back instead of repeated (`setResume(false)` runs everything
again). The summaries are rewritten after each instance. `getSummary()` and
`getRuns(i)` give the same numbers to the program. The executions run one after
the other, because each one already uses all the threads it is given.

### Several machines

`setShard(index, count)` divides the experiment in
`count` parts: the process only runs the executions of its part
(`index` from 0), so the same program can be started on `count` machines, or
`count` times on one, writing in the same folder or in folders that are joined
afterwards. `report()` runs nothing: it reads the executions that are in the
folder and writes all the reports (also useful to rebuild `runs.csv` with other
columns).

```cpp
exp.setShard(atoi(argv[1]), 4);   // this process: part argv[1] of 4
exp.run(runner);
// afterwards, once, in the folder with every file:
exp.report();
```

Example: [`exampleExperiment.cpp`](../example/functions/exampleExperiment.cpp).

## Comparing configurations (`Comparison.h`)

Each configuration is one experiment over the same instances.
`PTComparison` puts them side by side with the tests usually asked for in a
paper. The first configuration added is the control; the others are compared
with it.

```cpp
#include "include/Comparison.h"

PTExperiment fixed("results/fixed"), rate23("results/rate23");   // run() or report() on each
PTComparison cmp;
cmp.add("fixed", fixed);          // control
cmp.add("rate23", rate23);
cmp.setAlpha(0.05);               // level of the marks +, -, = (default)
cmp.write("results");             // comparison.txt, .csv, _total.csv, .tex
std::cout << cmp.table();
```

| Level | Test | Result |
|---|---|---|
| each instance | Mann-Whitney U (Wilcoxon rank-sum) on the values of the executions of the two configurations | p-value and a mark: `+` better than the control, `-` worse, `=` no significant difference |
| all the instances | Wilcoxon signed-rank on the means of the instances (paired by instance) | p-value of each configuration against the control; with more than one comparison, also corrected by the method of Holm |
| three or more configurations | Friedman on the ranks of the means inside each instance | mean rank of each configuration and one p-value |

The table also has, for each configuration, the mean gap (to the reference of
the instances if all have one, otherwise to the best mean of the instance),
the number of instances in which it has the best mean, and the count of `+`,
`=` and `-`. All the tests are two-sided and use ranks, so they do not suppose
normal values. The p-values are exact for small samples (Mann-Whitney up to 60
values in total, signed-rank up to 50 instances, ties included) and use the
normal approximation above that; Friedman uses the chi-square approximation.
They were checked against SciPy and against the enumeration of every case
(differences below 1e-14). With fewer than 6 instances the signed-rank test
cannot give p < 0.05, whatever the results. The tests are also functions of
their own: `ptstat::mannWhitney(a, b)`, `ptstat::wilcoxonSignedRank(a, b)`,
`ptstat::friedman(matrix, &meanRanks)`, `ptstat::holm(pValues)`.

Output of [`exampleComparison.cpp`](../example/functions/exampleComparison.cpp) (six TSP instances, ten
executions; 12 replicas with a fixed ladder as the control, 12 replicas with
update type 1, and 4 replicas with a fixed ladder):

```
instance          fixed         rate23         p    four_replicas         p
tsp50           5189.00        5189.00    1.0000 =        5298.30    0.0108 -
tsp70           6302.10        6311.00    0.6955 =        6453.20   <0.0001 -
tsp90           7130.10        7170.60    0.0785 =        7297.50   <0.0001 -
tsp110          7973.90        8025.30    0.3246 =        8136.40    0.0068 -
tsp130          9008.50        9005.00    0.8113 =        9172.50    0.0052 -
tsp150          9837.80        9778.60    0.1596 =        9910.70    0.3423 =

configuration    mean rank     gap %  best mean      +    =    -         p      Holm
fixed                 1.42      0.11          4                    control
rate23                1.58      0.23          3      0    6    0    0.8125    0.8125
four_replicas         3.00      2.02          0      0    1    5    0.0312    0.0625

Friedman test (the configurations have the same mean rank): p = 0.0087
```
