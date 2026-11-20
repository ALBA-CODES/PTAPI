/*
 * PTAPI v2 - dynamic temperature updates: the seven options of the 7th
 * argument of the constructor (the 8th is the number of exchange rounds
 * between two updates).
 *
 *   0  fixed temperatures
 *   1  exchange rate of 23%
 *   2  equal exchange rates (Hukushima)
 *   3  feedback-optimized, as in the thesis
 *   4  feedback-optimized by quantiles
 *   5  stochastic approximation towards 23.4% of exchanges
 *   6  equal rejection rates
 *
 * The example starts from a ladder that is deliberately bad (linear from 1 to
 * 60: almost no exchanges at the cold end) and prints what each rule does to
 * it. Types 1, 2, 4 and 6 keep the first and the last temperature; type 5
 * keeps the first one and never exceeds the initial maximum.
 *
 *   g++ exampleTempUpdate.cpp -std=c++17 -O3 -pthread -o exampleTempUpdate
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	const char* names[] = {"fixed", "23%", "equal rates", "feedback (thesis)", "feedback (quantiles)", "stochastic approximation", "equal rejection"};
	const int replicas = 10, updateInterval = 20;
	bool ok = true;

	for(int upType = 0; upType <= 6; upType++){
		PT<solTSP> algo(1.0f, 60.0f, replicas, 200, 2000, 1, upType, updateInterval);
		algo.setSeed(11);
		solTSP s = algo.start(0, &prob);

		double lo = 1.0, hi = 0.0;
		for(double r : algo.getSwapAcceptanceRates()){ lo = std::min(lo, r); hi = std::max(hi, r); }

		printf("type %d - %s: best %.0f, exchange rates from %.2f to %.2f\n", upType, names[upType], s.evalSol, lo, hi);
		if(upType == 0) printRow("   initial ladder ", algo.getInitialTemperatures());
		printRow("   final ladder   ", algo.getTemperatures());
		printRow("   exchange rate  ", algo.getSwapAcceptanceRates(), " %7.3f");
		ok = ok && prob.valid(s);
	}

	printf("\nThe exchange rates are those of the whole run, including the rounds before the ladder settled.\n");
	return ok ? 0 : 1;
}
