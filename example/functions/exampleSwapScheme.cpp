/*
 * PTAPI v2 - exchange scheme and the diagnostics of the ladder: setSwapScheme(),
 * getRoundTrips() and getCommunicationBarrier().
 *
 *   PT_SWAP_SEQUENTIAL  every pair proposes an exchange on every round (v1, default)
 *   PT_SWAP_EVEN_ODD    pairs (0,1),(2,3),... and (1,2),(3,4),... on alternate
 *                       rounds (non-reversible PT)
 *
 * A round trip is a replica that goes from the first temperature to the last
 * and comes back: it is what carries a solution found at high temperature to
 * the low ones. The communication barrier is the sum of the rejection rates
 * of the exchanges: it measures how hard the temperature interval is to cross.
 *
 *   g++ exampleSwapScheme.cpp -std=c++17 -O3 -pthread -o exampleSwapScheme
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	bool ok = true;

	printf("replicas  scheme       best   round trips   barrier   exchange rates (min..max)\n");
	for(int replicas : {6, 12, 24, 48}){
		for(int scheme : {PT_SWAP_SEQUENTIAL, PT_SWAP_EVEN_ODD}){
			PT<solTSP> algo(1.0f, 60.0f, replicas, 200, 2000, 4, 0, 50);
			algo.setSeed(17);
			algo.setSwapScheme(scheme);
			solTSP s = algo.start(0, &prob);

			double lo = 1.0, hi = 0.0;
			for(double r : algo.getSwapAcceptanceRates()){ lo = std::min(lo, r); hi = std::max(hi, r); }
			printf("%8d  %-10s %6.0f   %11lld   %7.2f   %.2f..%.2f\n", replicas, scheme == PT_SWAP_SEQUENTIAL ? "sequential" : "even-odd",
				s.evalSol, algo.getRoundTrips(), algo.getCommunicationBarrier(), lo, hi);
			ok = ok && prob.valid(s);
		}
	}

	printf("\nWith few replicas almost every exchange is rejected: the barrier is close to\n"
	       "its ceiling (replicas - 1) and there are no round trips. When the barrier\n"
	       "stops growing with the number of replicas, more replicas no longer help.\n");
	return ok ? 0 : 1;
}
