/*
 * PTAPI v2 - reproducible runs: setSeed(), getSeed() and the random functions
 * of the API (ptapi::rng(), ptapi::randInt(), ptapi::rand01()).
 *
 * The problem class takes its random numbers from the API (see tsp.h) instead
 * of creating its own engine. Each replica has its own engine, derived from
 * one seed, so a run that stops by PTL gives the same result for any number
 * of threads.
 *
 *   g++ exampleSeed.cpp -std=c++17 -O3 -pthread -o exampleSeed
 */
#include "tsp.h"

int main(){
	TSP prob(100);

	// PT(tempMin, tempMax, replicas, MCL, PTL, distribution, updateType, updateInterval)
	PT<solTSP> algo(1.0f, 60.0f, 12, 200, 500, 4, 0, 50);

	// same seed, different numbers of threads: same solution
	algo.setSeed(42);
	solTSP a = algo.start(1, &prob);
	solTSP b = algo.start(4, &prob);
	solTSP c = algo.start(0, &prob);   // 0 = all cores
	printf("seed %llu   1 thread: %.0f   4 threads: %.0f   all cores (%d): %.0f\n",
		algo.getSeed(), a.evalSol, b.evalSol, algo.getThreadsUsed(), c.evalSol);

	// another seed: another run
	algo.setSeed(43);
	solTSP d = algo.start(0, &prob);
	printf("seed %llu   %.0f\n", algo.getSeed(), d.evalSol);

	// without setSeed() the API draws a seed; getSeed() tells which one, so
	// that an interesting run can be repeated later
	PT<solTSP> free(1.0f, 60.0f, 12, 200, 500, 4, 0, 50);
	solTSP e = free.start(0, &prob);
	unsigned long long drawn = free.getSeed();
	free.setSeed(drawn);
	solTSP f = free.start(0, &prob);
	printf("drawn seed %llu   %.0f   repeated: %.0f\n", drawn, e.evalSol, f.evalSol);

	bool ok = (a.sol == b.sol) && (b.sol == c.sol) && (e.sol == f.sol);
	printf("%s\n", ok ? "the runs with the same seed are identical" : "ERROR: runs with the same seed differ");
	return ok ? 0 : 1;
}
