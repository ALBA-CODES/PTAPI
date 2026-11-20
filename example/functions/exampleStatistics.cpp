/*
 * PTAPI v2 - results of a run: getConvergence(), getEvaluations(),
 * getElapsedMs(), getTimeToBestMs(), getThreadsUsed(), getTemperatures(),
 * getMoveAcceptanceRates(), getSwapAcceptanceRates() and getReplicas().
 *
 * They describe the last start() and cost nothing during the run.
 *
 *   g++ exampleStatistics.cpp -std=c++17 -O3 -pthread -o exampleStatistics
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	PT<solTSP> algo(1.0f, 60.0f, 10, 200, 1000, 4, 0, 50);
	algo.setSeed(8);
	solTSP best = algo.start(0, &prob);

	// summary
	printf("best %.0f, found at round %d of %d, after %.1f ms of %.1f ms\n", best.evalSol, algo.getBestPTLSolIndex(), algo.getPTLEnd(),
		algo.getTimeToBestMs(), algo.getElapsedMs());
	printf("%lld evaluations (%.1f million/s) with %d threads, seed %llu\n\n", algo.getEvaluations(),
		algo.getEvaluations()/algo.getElapsedMs()/1000.0, algo.getThreadsUsed(), algo.getSeed());

	// one column per temperature
	std::vector<double> cost;
	for(const solTSP& r : algo.getReplicas()) cost.push_back(r.evalSol);
	printRow("temperature      ", algo.getTemperatures(), " %8.2f");
	printRow("accepted moves   ", algo.getMoveAcceptanceRates(), " %8.3f");
	printRow("final replica    ", cost, " %8.0f");
	printf(  "accepted exchanges    ");   // one per pair of adjacent temperatures
	for(double r : algo.getSwapAcceptanceRates()) printf(" %8.3f", r);
	printf("\n\n");

	// convergence curve: every new best solution {timeMs, ptl, value}
	const std::vector<PTImprovement>& curve = algo.getConvergence();
	printf("%zu improvements; some points of the curve (time, round, value):\n", curve.size());
	for(size_t i = 0; i < curve.size(); i += std::max<size_t>(1, curve.size()/8))
		printf("   %8.2f ms   round %4d   %7.0f\n", curve[i].timeMs, curve[i].ptl, curve[i].value);
	printf("   %8.2f ms   round %4d   %7.0f\n", curve.back().timeMs, curve.back().ptl, curve.back().value);

	bool ok = prob.valid(best) && curve.back().value == best.evalSol && algo.getEvaluations() == 10LL*200*1000;
	return ok ? 0 : 1;
}
