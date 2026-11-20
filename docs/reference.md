# Reference

Every public function, by header. `S` is the solution struct of the user.

## `Problem.h`

```cpp
struct solution{ double evalSol; };                     // base of the solution struct: the cost

template<typename S> class Problem{
    virtual S construction() = 0;                       // required: an evaluated initial solution
    virtual S neighbor(S sol) = 0;                      // required: one move
    virtual S constructionReplica(int index, int count);          // default: construction()
    virtual void neighborEval(const S& current, S& candidate);    // default: neighbor() + evaluate()

    void setMoveProbabilities(const std::vector<double>& weights);
    int  chooseMove();
    int  getMoveCount() const;
    std::vector<double> getMoveProbabilities() const;
};

// declared by the user class, in one of the forms:
S      evaluate(S sol);
double evaluate(S sol);
double evaluate(const S& sol);

ptapi::Rng& ptapi::rng();            // std::mt19937_64 of the replica running on this thread
int    ptapi::randInt(int lo, int hi);   // [lo, hi]
double ptapi::rand01();                  // [0, 1)
```

Chapter: [The problem class](problem.md).

## `PT.h`

### Construction and execution

| Function | |
|---|---|
| `PT(float tempMin, float tempMax, int replicas, int MCL, int PTL, int distribution, int updateType, int updateInterval)` | [Configuring the run](configuration.md) |
| `S start(int threads, P* problem)` | runs and returns the best solution; `threads <= 0`: all the cores |
| `S getBestSol()` | best solution of the last run |

### Configuration

| Function | Chapter |
|---|---|
| `void setSeed(unsigned long long)` | [configuration](configuration.md#seed) |
| `void setTemperatures(const std::vector<double>&)` | [configuration](configuration.md#a-ladder-of-your-own) |
| `void setInitialSolutions(const std::vector<S>&)` | [configuration](configuration.md#initial-solutions) |
| `void setParallelConstruction(bool)` | [configuration](configuration.md#initial-solutions) |
| `void setSwapScheme(int)` — `PT_SWAP_SEQUENTIAL`, `PT_SWAP_EVEN_ODD` | [configuration](configuration.md#exchange-scheme) |
| `void setCheckpoint(const std::string& file, SolutionWriter, double everySeconds = 0)` | [configuration](configuration.md#checkpoint-and-resume) |
| `bool loadCheckpoint(const std::string& file, SolutionReader)` | [configuration](configuration.md#checkpoint-and-resume) |
| `void saveCheckpoint(const std::string& file, SolutionWriter) const` | [configuration](configuration.md#checkpoint-and-resume) |
| `void setV1Compatibility(bool)` | [migrating](migrating-from-v1.md) |
| `void setProfiling(bool)` | [results](results.md#time-of-each-component) |

`SolutionWriter` is `std::function<void(std::ostream&, const S&)>` and
`SolutionReader` is `std::function<void(std::istream&, S&)>`.

### Stop criteria

| Function | |
|---|---|
| `void setTimeLimit(double seconds)` | |
| `void setTimeEnd(int seconds)` | version 1 |
| `void setTargetValue(double value)` | |
| `void setLowerBound(double bound, double gap = 0.0)` | |
| `void setStagnationRounds(int rounds)` | |
| `void setStagnationTime(double seconds)` | |
| `void setPTLStopEnd(float fraction)` | version 1 |
| `void setStopCondition(std::function<bool(const PTStatus&)>)` | `PTStatus{round, elapsedMs, bestValue, bestRound, bestMs}` |
| `void setOnImprove(std::function<void(const S&, int round, double ms)>)` | callback on every new best solution |
| `void requestStop()` | thread-safe |
| `int getStopReason() const`, `const char* getStopReasonName() const`, `static const char* stopReasonName(int)` | `PT_STOP_PTL`, `_TIME`, `_TARGET`, `_STAGNATION`, `_REQUEST`, `_CONDITION`, `_ERROR` |

Chapter: [Stop criteria](stopping.md).

### Temperatures

| Function | |
|---|---|
| `void setAutoUpdate(int firstInterval, double growth = 2.0)` | used with `PT_UPDATE_AUTO` as the update interval |
| `std::pair<double,double> calibrateTemperatures(P* problem, int samples = 5000, double acceptHot = 0.02, double acceptCold = 0.001)` | returns `{Tmin, Tmax}` |
| `std::vector<double> calibrateLadder(P* problem, int threads = 0, double swapRate = 0.23, int rounds = 200, int iterations = 10, double tMin = 0.0)` | returns the ladder |
| `std::deque<double> tempLinear / tempILinear / tempExp / tempPG(float tempMin, float tempMax, int replicas)` | the four distributions (version 1) |

Chapter: [Temperatures](temperatures.md).

### Results

| Function | Type |
|---|---|
| `getBestPTLSolIndex()`, `getPTLEnd()` | `int` |
| `getSeed()` | `unsigned long long` |
| `getEvaluations()`, `getRoundTrips()` | `long long` |
| `getElapsedMs()`, `getTimeToBestMs()`, `getCommunicationBarrier()` | `double` |
| `getThreadsUsed()` | `int` |
| `getInitialTemperatures()`, `getTemperatures()`, `getMoveAcceptanceRates()`, `getSwapAcceptanceRates()` | `std::vector<double>` |
| `getReplicas()` | `std::vector<S>` |
| `getConvergence()` | `std::vector<PTImprovement>` — `{timeMs, ptl, value}` |
| `getMoveStatistics()` | `std::vector<PTMoveStat>` — `{proposed, accepted}` |
| `getTiming()` | `PTTiming` — `{wallMs, constructionMs, chainsMs, exchangesMs, updatesMs, waitingMs, threads}` |

Chapter: [Results of a run](results.md). `PTAPI_VERSION_MAJOR` and
`PTAPI_VERSION_MINOR` give the version of the headers.

## `Experiment.h` (optional)

```cpp
PTRunResult ptResult(PT<S>& pt, const S& best);      // everything the API measured in the last start()
```

`PTRunResult` fields: `value`, `seed`, `elapsedMs`, `timeToBestMs`, `rounds`,
`bestRound`, `evaluations`, `threads`, `stopReason`, `roundTrips`, `timing`,
`convergence`, `temperatures`, `initialTemperatures`, `swapRates`, `moveRates`,
`moves`, `improvements`, `solution` (text written by the user), and
`add(name, value)` for values of the user.

| `PTExperiment` | |
|---|---|
| `PTExperiment(const std::string& outputDir)` | |
| `void addInstance(const std::string& path, const std::string& name = "")` | |
| `void setReference(double value)` | best known value, lower bound or target of the last instance added |
| `void addInstancesFromFile(const std::string& listFile)` | one `path [reference]` per line |
| `void setRuns(int)`, `void setFirstSeed(unsigned long long)` | defaults 10 and 1 |
| `void setResume(bool)`, `void setVerbose(bool)` | defaults `true` |
| `void setCsvSeparator(char)` | default `,` |
| `void setCsvColumns(const std::vector<std::string>&)`, `void addCsvColumns(...)`, `static std::vector<std::string> csvColumnNames()` | columns of `runs.csv` |
| `void setSaveConvergence(bool)`, `void setCurvePoints(int)` | convergence files; defaults `true` and 100 |
| `void setShard(int index, int count)` | this process runs part `index` of `count` |
| `void run(Runner)` | `Runner` is `std::function<PTRunResult(const std::string& instance, unsigned long long seed)>` |
| `void report()` | writes the reports from the files of the folder, without running |
| `getSummary()` | `std::vector<PTInstanceSummary>` — `name`, `runs`, `best`, `mean`, `median`, `sd`, `worst`, `meanSeconds`, `meanSecondsToBest`, `meanEvaluations`, `hasReference`, `reference`, `gapBest`, `gapMean`, `hits`, `meanSecondsToReference` |
| `getRuns(size_t i)` | `std::vector<PTRunResult>` of instance `i` |
| `summaryTable()`, `getOutputDir()` | `std::string` |

## `Comparison.h` (optional)

| `PTComparison` | |
|---|---|
| `void add(const std::string& name, const PTExperiment&)` | the first one added is the control |
| `void setAlpha(double)` | default 0.05 |
| `void setCsvSeparator(char)` | |
| `void write(const std::string& dir)` | `comparison.txt`, `comparison.csv`, `comparison_total.csv`, `comparison.tex` |
| `std::string table()` | |
| `getInstances()` | `std::vector<std::string>` |
| `getCells()` | `[instance][configuration]` of `PTComparisonCell{runs, best, mean, median, sd, meanSeconds, pValue, sign}` |
| `getTotals()` | `PTComparisonTotal{name, meanRank, meanGap, bestMean, wins, ties, losses, pValue, pHolm}` per configuration |
| `getFriedmanP()` | `double` |

```cpp
double ptstat::mannWhitney(const std::vector<double>& a, const std::vector<double>& b);         // two-sided p
double ptstat::wilcoxonSignedRank(const std::vector<double>& a, const std::vector<double>& b);  // paired, two-sided p
double ptstat::friedman(const std::vector<std::vector<double>>& values, std::vector<double>* meanRanks = nullptr);
std::vector<double> ptstat::holm(const std::vector<double>& p);
std::vector<double> ptstat::ranks(const std::vector<double>& v, double* tieTerm = nullptr);
double ptstat::chiSquareUpper(double x, double df);
```

Chapter: [Experiments](experiments.md).

## `Options.h` (optional)

| `PTOptions` | |
|---|---|
| `PTOptions(int argc, char* argv[])` | throws `std::invalid_argument` on a bad value |
| `instance` | first argument that is not an option |
| `bool has(name)`, `get(name, default)` | options of the problem (`std::string`, `double` or `int` default) |
| `std::unique_ptr<PT<S>> make<S>() const` | the configured object, not started |
| `void apply(PT<S>&) const` | configures an object built by the user |
| `S run<S>(P* problem)` | builds, calibrates if asked, runs; fills `result` |
| `double tunerCostValue() const`, `bool printForTuner() const` | cost for irace |

Chapter: [Command line and irace](tuning.md).
