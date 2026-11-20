/*
 * PTAPI v2 - comparison of configurations with statistical tests: PTComparison
 * (header ../../include/Comparison.h).
 *
 * Each configuration is one experiment (PTExperiment) over the same instances.
 * The first one added to the comparison is the control; the others are
 * compared with it:
 *
 *   in each instance      Mann-Whitney U test on the values of the executions
 *                         (+ better than the control, - worse, = no difference)
 *   over the instances    Wilcoxon signed-rank test on the means of the
 *                         instances, also corrected by the method of Holm when
 *                         there is more than one comparison
 *   three or more         mean rank of each configuration and Friedman test
 *
 * Files written: comparison.txt, comparison.csv, comparison_total.csv and
 * comparison.tex. The tests are also functions of their own (namespace ptstat).
 *
 *   g++ exampleComparison.cpp -std=c++17 -O3 -pthread -o exampleComparison
 *   ./exampleComparison [outputFolder] [runs]
 */
#include "tsp.h"
#include "../../include/Comparison.h"

// One experiment: the same instances and seeds, a different configuration of the PT.
static void runExperiment(PTExperiment& exp, int runs, int replicas, int updateType){
	exp.setRuns(runs);
	exp.setVerbose(false);
	for(int cities : {50, 70, 90, 110, 130, 150}) exp.addInstance(std::to_string(cities), "tsp" + std::to_string(cities));
	exp.run([=](const std::string& instance, unsigned long long seed){
		TSP prob(std::stoi(instance));
		PT<solTSP> algo(1.0f, 60.0f, replicas, 200, 300, 4, updateType, 50);
		algo.setSeed(seed);
		solTSP best = algo.start(0, &prob);
		return ptResult(algo, best);
	});
}

int main(int argc, char* argv[]){
	const std::string out = argc > 1 ? argv[1] : "exampleComparison.out";
	const int runs = argc > 2 ? atoi(argv[2]) : 10;
	const bool keep = argc > 1;   // the default folder is removed at the end ('make run')

	// One folder for each configuration (the experiments can also come from
	// different programs or days: PTExperiment::report() reads a folder back).
	PTExperiment fixed(out + "/fixed"), rate23(out + "/rate23"), few(out + "/four_replicas");
	runExperiment(fixed, runs, 12, 0);    // 12 replicas, fixed ladder
	runExperiment(rate23, runs, 12, 1);   // 12 replicas, update type 1
	runExperiment(few, runs, 4, 0);       // 4 replicas, fixed ladder

	PTComparison cmp;
	cmp.add("fixed", fixed);              // the first is the control
	cmp.add("rate23", rate23);
	cmp.add("four_replicas", few);
	cmp.setAlpha(0.05);                   // level of the marks +, -, = (default 0.05)
	cmp.write(out);                       // comparison.txt, .csv, _total.csv and .tex
	std::cout << cmp.table();

	// The same numbers in the program
	bool ok = cmp.getInstances().size() == 6 && cmp.getTotals().size() == 3;
	for(const PTComparisonTotal& t : cmp.getTotals()) ok = ok && t.meanRank >= 1.0 && t.meanRank <= 3.0;
	const PTComparisonTotal& f = cmp.getTotals()[2];
	ok = ok && f.wins + f.ties + f.losses == 6 && f.pValue >= 0.0 && f.pValue <= 1.0;

	// The tests on their own: two samples of executions, and paired values
	std::vector<double> a, b;
	for(const PTRunResult& r : fixed.getRuns(5)) a.push_back(r.value);
	for(const PTRunResult& r : few.getRuns(5)) b.push_back(r.value);
	printf("\ntsp150, fixed against four_replicas: Mann-Whitney p = %.4f\n", ptstat::mannWhitney(a, b));
	printf("paired example: Wilcoxon signed-rank p = %.4f\n", ptstat::wilcoxonSignedRank({10, 12, 9, 14, 11, 13}, {11, 14, 12, 18, 16, 19}));
	ok = ok && std::fabs(ptstat::wilcoxonSignedRank({10, 12, 9, 14, 11, 13}, {11, 14, 12, 18, 16, 19}) - 2.0/64.0) < 1e-12;

	if(!keep) std::filesystem::remove_all(out);
	return ok ? 0 : 1;
}
