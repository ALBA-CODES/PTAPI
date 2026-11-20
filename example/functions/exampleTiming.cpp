/*
 * PTAPI v2 - time of each component: setProfiling() and getTiming().
 *
 * Where does the time of a run go: the problem code (the chains), the exchange
 * proposals, the temperature updates, or threads waiting for work? The answer
 * tells what is worth optimizing: evaluate()/neighbor(), or the parameters
 * (a small MCL makes the threads wait more for each other).
 *
 * The times of the components are summed over the threads ("thread-ms").
 * Measuring costs two clock readings per node execution (not per move).
 *
 *   g++ exampleTiming.cpp -std=c++17 -O3 -pthread -o exampleTiming
 */
#include "tsp.h"

static bool run(const char* title, int replicas, int MCL, int PTL, int threads){
	TSP prob(100);
	PT<solTSP> algo(1.0f, 60.0f, replicas, MCL, PTL, 4, 6, 20);
	algo.setSeed(5);
	algo.setProfiling(true);
	solTSP s = algo.start(threads, &prob);

	const PTTiming& t = algo.getTiming();
	const double total = t.chainsMs + t.exchangesMs + t.updatesMs + t.waitingMs;
	printf("%s (%d replicas, MCL %d, %d threads): best %.0f in %.1f ms\n", title, replicas, MCL, t.threads, s.evalSol, t.wallMs);
	printf("   construction          %9.2f ms\n", t.constructionMs);
	printf("   chains (problem code) %9.2f thread-ms  %5.1f%%\n", t.chainsMs, 100.0*t.chainsMs/total);
	printf("   exchanges             %9.2f thread-ms  %5.1f%%\n", t.exchangesMs, 100.0*t.exchangesMs/total);
	printf("   temperature updates   %9.2f thread-ms  %5.1f%%\n", t.updatesMs, 100.0*t.updatesMs/total);
	printf("   waiting + scheduling  %9.2f thread-ms  %5.1f%%\n", t.waitingMs, 100.0*t.waitingMs/total);
	return prob.valid(s) && t.chainsMs > 0.0 && t.wallMs > 0.0;
}

int main(){
	bool ok = true;
	// the same number of evaluations in the three runs
	ok = run("long chains", 8, 2000, 250, 8) && ok;
	ok = run("\nshort chains", 8, 20, 25000, 8) && ok;
	ok = run("\nfewer threads than replicas", 8, 2000, 250, 2) && ok;

	printf("\nWith short chains the threads meet at the exchanges much more often, and the\n"
	       "share of the time spent waiting grows: MCL is also a parallel efficiency parameter.\n");
	return ok ? 0 : 1;
}
