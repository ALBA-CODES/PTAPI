/*
 * PTAPI - quick start: number partitioning.
 *
 * Divide N numbers in two groups whose sums are as close as possible. The whole
 * program is here: the three steps of the API and a few results of the run.
 * It is the program explained line by line in docs/getting-started.md.
 *
 *   g++ quickstart.cpp -std=c++17 -O3 -pthread -o quickstart
 */
#include <cstdio>
#include <cmath>
#include "../include/PT.h"

// Step 1: the solution. It extends 'solution', which has the field evalSol (the cost).
struct solPart: public solution{
	std::vector<char> side;   // side[i] = 0 or 1: the group of number i
	double diff = 0;          // sum of group 0 minus sum of group 1
};

// Step 2: the problem. Everything that depends on the problem is in this class.
class Partition: public Problem<solPart>{
	private:
		std::vector<double> num;
	public:
		Partition(int n, unsigned instanceSeed){
			std::mt19937 g(instanceSeed);
			std::uniform_int_distribution<int> u(1, 1000000);
			for(int i = 0; i < n; i++) num.push_back(u(g));
		}

		// initial solution, already evaluated; ptapi::randInt uses the engine of the replica
		solPart construction(){
			solPart s;
			for(size_t i = 0; i < num.size(); i++) s.side.push_back((char)ptapi::randInt(0, 1));
			s.evalSol = evaluate(s);
			s.diff = difference(s);
			return s;
		}

		// neighbor: one number changes group
		solPart neighbor(solPart s){
			int i = ptapi::randInt(0, (int)num.size()-1);
			s.side[i] = 1 - s.side[i];
			s.diff = difference(s);
			return s;
		}

		// cost: |sum of group 0 - sum of group 1| (the API minimizes)
		double evaluate(const solPart& s){ return std::fabs(difference(s)); }

		double difference(const solPart& s) const {
			double d = 0;
			for(size_t i = 0; i < num.size(); i++) d += s.side[i] ? -num[i] : num[i];
			return d;
		}

		// Optional: one step without copies and with delta evaluation, O(1) instead of O(N)
		void neighborEval(const solPart& current, solPart& cand){
			int i = ptapi::randInt(0, (int)num.size()-1);
			cand.side = current.side;                     // reuses the memory of 'cand'
			cand.side[i] = 1 - cand.side[i];
			cand.diff = current.diff + (current.side[i] ? 2.0 : -2.0)*num[i];
			cand.evalSol = std::fabs(cand.diff);
		}
};

int main(){
	Partition prob(60, 2026);

	// Step 3: the PT.
	// PT(Tmin, Tmax, replicas, MCL, PTL, distribution, updateType, updateInterval)
	PT<solPart> pt(10.0f, 1000000.0f, 16, 400, 2000, 4, 0, 20);
	pt.setSeed(42);                                       // reproducible run
	pt.setTimeLimit(5.0);                                 // at most 5 seconds

	solPart best = pt.start(0, &prob);                    // 0 = all the cores

	printf("best difference: %.0f (round %d of %d, %.1f ms of %.1f ms)\n", best.evalSol,
		pt.getBestPTLSolIndex(), pt.getPTLEnd(), pt.getTimeToBestMs(), pt.getElapsedMs());
	printf("%lld evaluations on %d threads, seed %llu, stop: %s\n", pt.getEvaluations(), pt.getThreadsUsed(), pt.getSeed(), pt.getStopReasonName());
	printf("round trips: %lld   barrier: %.2f\n", pt.getRoundTrips(), pt.getCommunicationBarrier());
	printf("final ladder :"); for(double x : pt.getTemperatures()) printf(" %.0f", x); printf("\n");
	printf("exchange rate:"); for(double x : pt.getSwapAcceptanceRates()) printf(" %.2f", x); printf("\n");
	return std::fabs(prob.evaluate(best) - best.evalSol) < 1e-6 ? 0 : 1;
}
