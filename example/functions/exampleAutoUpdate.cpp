/*
 * PTAPI v2 - automatic interval between temperature updates: PT_UPDATE_AUTO
 * and setAutoUpdate().
 *
 * The 8th argument of the constructor is the number of exchange proposals
 * between two updates of the temperatures. With PT_UPDATE_AUTO (0) the API
 * chooses: the first update happens after 10 proposals and each interval is
 * twice the previous one (10, 20, 40, ...). The ladder is corrected often at
 * the beginning, when it is far from a good one, and later with statistics of
 * longer windows; there is one parameter less to tune.
 *
 *   g++ exampleAutoUpdate.cpp -std=c++17 -O3 -pthread -o exampleAutoUpdate
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	bool ok = true;
	const int upType = 6;   // equal rejection rates; any type from 1 to 6 works the same way

	printf("interval          best   ladder after the run\n");
	for(int interval : {5, 20, 100, 500, PT_UPDATE_AUTO}){
		PT<solTSP> algo(1.0f, 60.0f, 10, 200, 2000, 1, upType, interval);   // linear ladder: a bad start
		algo.setSeed(11);
		solTSP s = algo.start(0, &prob);
		if(interval == PT_UPDATE_AUTO) printf("automatic      "); else printf("fixed, %-7d ", interval);
		printf("%7.0f  ", s.evalSol);
		for(double t : algo.getTemperatures()) printf(" %6.2f", t);
		printf("\n");
		ok = ok && prob.valid(s);
	}

	// the first interval and the growth factor can be changed
	PT<solTSP> algo(1.0f, 60.0f, 10, 200, 2000, 1, upType, PT_UPDATE_AUTO);
	algo.setAutoUpdate(50, 1.5);   // 50, 75, 113, ...
	algo.setSeed(11);
	solTSP s = algo.start(0, &prob);
	printf("automatic 50x1.5 %5.0f  ", s.evalSol);
	for(double t : algo.getTemperatures()) printf(" %6.2f", t);
	printf("\n\nOne run of each: to compare the options, average several seeds (see exampleExperiment).\n");
	return (ok && prob.valid(s)) ? 0 : 1;
}
