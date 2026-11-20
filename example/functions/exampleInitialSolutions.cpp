/*
 * PTAPI v2 - initial solutions: setInitialSolutions(), constructionReplica()
 * and setParallelConstruction().
 *
 * Two ways of starting from something better than random solutions:
 *
 *  - setInitialSolutions(): the replicas start from solutions given by the
 *    user (from a heuristic, from a previous run, from a file). The API
 *    evaluates them; if there are fewer solutions than replicas they are
 *    assigned cyclically.
 *
 *  - constructionReplica(index, count): the problem class builds the solution
 *    of each replica knowing its position in the ladder, so the heuristic is
 *    part of the problem and nothing has to be prepared in main(). With
 *    setParallelConstruction(true) the replicas are built at the same time,
 *    which pays off when the heuristic is expensive.
 *
 *   g++ exampleInitialSolutions.cpp -std=c++17 -O3 -pthread -o exampleInitialSolutions
 */
#include "tsp.h"

// The problem of tsp.h plus a constructive heuristic
class TSPGreedy: public TSP{
	public:
		using TSP::TSP;

		// nearest neighbour tour starting at city 'first'
		solTSP nearestNeighbour(int first){
			solTSP s;
			std::vector<char> used(n, 0);
			int cur = first;
			for(int k = 0; k < n; k++){
				s.sol.push_back(cur);
				used[cur] = 1;
				int next = -1;
				for(int j = 0; j < n; j++) if(!used[j] && (next < 0 || dist(cur, j) < dist(cur, next))) next = j;
				cur = next;
			}
			return s;   // not evaluated: the API does it
		}
};

// The heuristic inside the problem: the colder half of the ladder starts from
// nearest neighbour tours (a different first city for each replica) and the
// hotter half from random tours, which keeps diversity.
class TSPHeuristicStart: public TSPGreedy{
	public:
		using TSPGreedy::TSPGreedy;

		solTSP constructionReplica(int index, int count) override {
			if(index >= count/2) return construction();      // random (TSP::construction)
			solTSP s = nearestNeighbour(index % n);
			s.evalSol = evaluate(s);                          // it must return an evaluated solution
			return s;
		}
};

int main(){
	TSPGreedy prob(200);
	const int replicas = 12, MCL = 200, PTL = 150;   // a short run, where the start matters

	// random start (construction() of the problem)
	PT<solTSP> cold(1.0f, 30.0f, replicas, MCL, PTL, 4, 0, 50);
	cold.setSeed(5);
	solTSP a = cold.start(0, &prob);

	// one nearest neighbour tour for each replica
	std::vector<solTSP> starts;
	for(int i = 0; i < replicas; i++) starts.push_back(prob.nearestNeighbour(i));
	PT<solTSP> warm(1.0f, 30.0f, replicas, MCL, PTL, 4, 0, 50);
	warm.setSeed(5);
	warm.setInitialSolutions(starts);
	solTSP b = warm.start(0, &prob);

	// a second stage that starts from the final replicas of the first one
	// (getReplicas() returns the solution of each temperature)
	PT<solTSP> stage2(1.0f, 30.0f, replicas, MCL, PTL, 4, 0, 50);
	stage2.setSeed(6);
	stage2.setInitialSolutions(warm.getReplicas());
	solTSP c = stage2.start(0, &prob);

	printf("nearest neighbour tour (city 0)  %6.0f\n", prob.evaluate(starts[0]));
	printf("random start                     %6.0f   first best at %.0f\n", a.evalSol, cold.getConvergence().front().value);
	printf("start from nearest neighbour     %6.0f   first best at %.0f\n", b.evalSol, warm.getConvergence().front().value);
	printf("second stage from its replicas   %6.0f\n", c.evalSol);

	// the same idea with the heuristic inside the problem class, built in parallel
	TSPHeuristicStart hprob(200);
	PT<solTSP> heur(1.0f, 30.0f, replicas, MCL, PTL, 4, 0, 50);
	heur.setSeed(5);
	heur.setParallelConstruction(true);
	solTSP d = heur.start(0, &hprob);
	printf("constructionReplica (parallel)   %6.0f   first best at %.0f   construction took %.2f ms\n", d.evalSol,
		heur.getConvergence().front().value, heur.getTiming().constructionMs);

	bool ok = prob.valid(a) && prob.valid(b) && prob.valid(c) && hprob.valid(d);
	return ok ? 0 : 1;
}
