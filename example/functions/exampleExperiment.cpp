/*
 * PTAPI v2 - management of experiments: PTExperiment, PTRunResult and ptResult()
 * (header ../../include/Experiment.h).
 *
 * The usual experiment of a paper: R executions for each instance, one file
 * per execution and a summary with the best, the mean, the median and the
 * standard deviation of each instance. The user writes only the function that
 * runs ONE execution; the API does the loop, the files and the statistics.
 *
 * Files written under the output folder:
 *   <instance>/run_<k>.txt   one execution: seed, value, times, rounds, stop
 *                            reason, ladder, rates, convergence curve, solution
 *   <instance>/convergence.csv        convergence of the instance over its
 *                            executions, by round: best, median, mean and worst
 *                            of the values the executions had at each point
 *   <instance>/convergence_time.csv   the same by time (seconds)
 *   <instance>/time_to_target.csv     instances with a reference: time each
 *                            execution took to reach it and the empirical
 *                            probability of having reached it by then
 *   figures.gp, figures.tex  scripts that draw the three files above, with
 *                            gnuplot (gnuplot figures.gp) or with LaTeX and
 *                            pgfplots (lualatex figures.tex), run in the folder
 *   runs.csv                 one line per execution
 *   summary.csv              one line per instance
 *   summary.txt              the same table, aligned
 *   summary.tex              the same table, for LaTeX (booktabs)
 *
 * An interrupted experiment continues where it stopped: executions whose file
 * exists are not repeated. The experiment can also be divided among several
 * processes or machines (setShard) and the reports rebuilt from the files of
 * the folder (report), as shown at the end.
 *
 *   g++ exampleExperiment.cpp -std=c++17 -O3 -pthread -o exampleExperiment
 *   ./exampleExperiment [outputFolder] [runs]
 */
#include <sstream>
#include "tsp.h"
#include "../../include/Experiment.h"

int main(int argc, char* argv[]){
	const std::string out = argc > 1 ? argv[1] : "exampleExperiment.out";
	const int runs = argc > 2 ? atoi(argv[2]) : 10;
	const bool keep = argc > 1;   // the default folder is removed at the end ('make run')

	PTExperiment exp(out);
	exp.setRuns(runs);            // executions per instance
	exp.setFirstSeed(1);          // execution k of every instance uses the seed 1 + k

	// Instances: what the runner receives (usually the path of the instance file;
	// here the number of cities) and, optionally, the name used in the reports.
	exp.addInstance("60", "tsp60");
	exp.setReference(6020.0);     // a target value also works: time_to_target.csv tells how long the executions took to reach it
	exp.addInstance("100", "tsp100");
	exp.setReference(7350.0);     // best known value or lower bound of the last instance: adds the gap
	exp.addInstance("150", "tsp150");
	// exp.addInstancesFromFile("instances.txt");   // one "path [reference]" per line

	// Columns of runs.csv. The default is instance, run, seed, value, seconds,
	// seconds_to_best, rounds, best_round, evaluations, threads, stop, round_trips;
	// addCsvColumns() adds others and setCsvColumns() replaces the whole list.
	// PTExperiment::csvColumnNames() lists every column the API can fill.
	exp.addCsvColumns({"improvements", "barrier", "swap_rate_mean", "temp_last", "seconds_to_reference"});

	// One execution: build the problem, configure the PT, run, return the result.
	exp.run([](const std::string& instance, unsigned long long seed){
		TSP prob(std::stoi(instance));
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, 300, 4, 0, 50);
		algo.setSeed(seed);
		solTSP best = algo.start(0, &prob);

		PTRunResult r = ptResult(algo, best);      // everything the API measured
		std::ostringstream tour;                    // optional: the solution, as text
		for(int c : best.sol) tour << c << ' ';
		r.solution = tour.str();
		r.add("cities", (double)prob.size());      // optional: values of the user, one more column each
		return r;
	});

	// The same numbers are available in the program
	bool ok = true;
	for(const PTInstanceSummary& s : exp.getSummary()){
		ok = ok && s.runs == runs && s.best <= s.median && s.median <= s.worst && s.sd >= 0.0;
		if(s.hasReference) printf("%s: gap of the best value to the reference %.0f: %.2f%%\n", s.name.c_str(), s.reference, s.gapBest);
	}

	// first lines of runs.csv
	{
		std::ifstream csv((out + "/runs.csv").c_str());
		std::string line;
		printf("\nruns.csv:\n");
		for(int i = 0; i < 3 && std::getline(csv, line); i++) printf("   %s\n", line.c_str());
	}

	// convergence of an instance and times to the target, ready to be plotted
	for(const char* name : {"/tsp100/convergence.csv", "/tsp60/time_to_target.csv"}){
		std::ifstream f((out + name).c_str());
		std::string line;
		printf("\n%s:\n", name + 1);
		for(int i = 0; i < 4 && std::getline(f, line); i++) printf("   %s\n", line.c_str());
		ok = ok && f.good();
	}
	for(const PTInstanceSummary& s : exp.getSummary())
		if(s.hasReference && s.hits) printf("%s: %d of %d executions reached %.0f, in %.4f s on average\n", s.name.c_str(), s.hits, s.runs, s.reference, s.meanSecondsToReference);

	// The same experiment divided in two parts. Each part runs only its executions
	// (here one after the other; in practice two processes or two machines writing
	// in the same folder, or in folders that are joined afterwards). report() runs
	// nothing: it reads the files of the folder and writes the complete reports.
	auto describe = [&](PTExperiment& e){
		e.setRuns(runs);
		e.setVerbose(false);
		e.addInstance("60", "tsp60");
		e.addInstance("100", "tsp100");
	};
	auto one = [](const std::string& instance, unsigned long long seed){
		TSP prob(std::stoi(instance));
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, 300, 4, 0, 50);
		algo.setSeed(seed);
		solTSP best = algo.start(0, &prob);
		return ptResult(algo, best);
	};
	const std::string parts = out + "/parts";
	PTExperiment part0(parts), part1(parts), whole(parts);
	describe(part0); part0.setShard(0, 2); part0.run(one);
	describe(part1); part1.setShard(1, 2); part1.run(one);
	describe(whole); whole.report();
	printf("\ntwo parts: %d + %d executions; report() found %d\n", part0.getSummary()[0].runs + part0.getSummary()[1].runs,
	       part1.getSummary()[0].runs + part1.getSummary()[1].runs - part0.getSummary()[0].runs - part0.getSummary()[1].runs,
	       whole.getSummary()[0].runs + whole.getSummary()[1].runs);
	ok = ok && whole.getSummary()[0].runs == runs && whole.getSummary()[1].runs == runs;

	if(!keep) std::filesystem::remove_all(out);
	return ok ? 0 : 1;
}
