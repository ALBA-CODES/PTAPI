/*
 * PTAPI v2 - several moves with probabilities: setMoveProbabilities(),
 * chooseMove() and getMoveStatistics().
 *
 * The problem registers one weight for each kind of move and asks the API
 * which one to apply. After the run the API reports how many times each move
 * was proposed and accepted.
 *
 *   g++ exampleMoves.cpp -std=c++17 -O3 -pthread -o exampleMoves
 */
#include "tsp.h"

// TSP with three moves: 2-opt (reverse a segment), swap of two cities and
// insertion of a city in another position
class TSPMoves: public Problem<solTSP>{
	private:
		TSP base;
		int n;
	public:
		TSPMoves(int n_, const std::vector<double>& weights): base(n_), n(n_){
			setMoveProbabilities(weights);   // weights: they are normalized by the API
		}

		solTSP construction(){ return base.construction(); }
		double evaluate(const solTSP& s){ return base.evaluate(s); }
		bool valid(const solTSP& s){ return base.valid(s); }

		solTSP neighbor(solTSP s){
			int a, b;
			do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);

			switch(chooseMove()){
				case 0:   // 2-opt
					if(a > b) std::swap(a, b);
					std::reverse(s.sol.begin()+a, s.sol.begin()+b+1);
					break;
				case 1:   // swap
					std::swap(s.sol[a], s.sol[b]);
					break;
				default:{ // insertion: the city of position a goes to position b
					int city = s.sol[a];
					s.sol.erase(s.sol.begin()+a);
					s.sol.insert(s.sol.begin()+b, city);
				}
			}
			return s;
		}
};

static solTSP run(const char* title, const std::vector<double>& weights, bool& ok){
	const char* names[] = {"2-opt", "swap", "insertion"};
	TSPMoves prob(100, weights);
	PT<solTSP> algo(1.0f, 60.0f, 12, 200, 1000, 4, 0, 50);
	algo.setSeed(21);
	solTSP s = algo.start(0, &prob);

	printf("%s: best %.0f\n", title, s.evalSol);
	std::vector<double> p = prob.getMoveProbabilities();
	const std::vector<PTMoveStat>& st = algo.getMoveStatistics();
	long long total = 0;
	for(size_t m = 0; m < st.size(); m++){
		printf("   %-10s probability %.2f   proposed %9lld   accepted %8lld (%4.1f%%)\n", names[m], p[m],
			st[m].proposed, st[m].accepted, st[m].proposed ? 100.0*st[m].accepted/st[m].proposed : 0.0);
		total += st[m].proposed;
	}
	ok = ok && prob.valid(s) && total == algo.getEvaluations();
	return s;
}

int main(){
	bool ok = true;
	run("only 2-opt", {1, 0, 0}, ok);
	run("\nequal probabilities", {1, 1, 1}, ok);
	run("\n60% 2-opt, 10% swap, 30% insertion", {0.6, 0.1, 0.3}, ok);

	printf("\nThe acceptance counts all the temperatures. A move that is rarely accepted\n"
	       "may not be worth its share of the evaluations.\n");
	return ok ? 0 : 1;
}
