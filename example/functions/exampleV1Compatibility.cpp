/*
 * PTAPI v2 - code written for v1 and setV1Compatibility().
 *
 * 1) A problem class written exactly as in v1 (double evaluate(S), its own
 *    random engine, Nup/Ndown copied by hand) compiles and runs unchanged.
 * 2) setV1Compatibility(true) restores the v1 behaviours that v2 changed: the
 *    energies truncated to int in the exchange test, and the rules of the
 *    temperature update types 1 and 2 that did not keep the ends of the
 *    ladder. Use it to reproduce results published with v1.
 *
 *   g++ exampleV1Compatibility.cpp -std=c++17 -O3 -pthread -o exampleV1Compatibility
 */
#include "tsp.h"

// A class in the v1 style
class TSPv1: public Problem<solTSP>{
	private:
		TSP base;
		int n;
	public:
		TSPv1(int n_): base(n_), n(n_){}

		solTSP construction(){
			solTSP s;
			std::random_device rnd_device;
			std::mt19937 mersenne_engine{rnd_device()};
			for(int i = 0; i < n; i++) s.sol.push_back(i);
			std::shuffle(s.sol.begin(), s.sol.end(), mersenne_engine);
			s.evalSol = evaluate(s);
			return s;
		}

		solTSP neighbor(solTSP sol){
			solTSP s;
			std::random_device rnd_device;
			std::mt19937 mersenne_engine{rnd_device()};
			std::uniform_int_distribution<int> dist(0, n-1);
			int a = dist(mersenne_engine), b = dist(mersenne_engine);
			s.sol = sol.sol;
			s.Nup = sol.Nup;       // needed in v1; ignored in v2
			s.Ndown = sol.Ndown;
			if(a > b) std::swap(a, b);
			std::reverse(s.sol.begin()+a, s.sol.begin()+b+1);
			return s;
		}

		double evaluate(solTSP s){ return base.evaluate(s); }
};

int main(){
	bool ok = true;

	// 1) v1 code, v1 calls
	TSPv1* prob = new TSPv1(100);
	PT<solTSP>* algo = new PT<solTSP>(1.0f, 60.0f, 12, 200, 500, 4, 1, 50);
	algo->setPTLStopEnd(0.5f);
	solTSP s = algo->start(4, prob);
	printf("class written for v1:        best %.0f, rounds %d, best at round %d\n", s.evalSol, algo->getPTLEnd(), algo->getBestPTLSolIndex());
	delete algo;
	delete prob;

	// 2) update type 1 (23%) with the v2 rule and with the v1 rule
	TSP p2(100);
	for(int v1 = 0; v1 <= 1; v1++){
		PT<solTSP> a(1.0f, 60.0f, 12, 200, 2000, 4, 1, 20);
		a.setSeed(4);
		a.setV1Compatibility(v1 == 1);
		s = a.start(0, &p2);
		printf("\nupdate type 1, %s: best %.0f\n", v1 ? "v1 rule (setV1Compatibility)" : "v2 rule", s.evalSol);
		printRow("   initial ladder", a.getInitialTemperatures());
		printRow("   final ladder  ", a.getTemperatures());
		ok = ok && p2.valid(s);
		if(!v1) ok = ok && std::fabs(a.getTemperatures().back() - 60.0) < 1e-3;
	}
	printf("\nWith the v1 rule the ladder drifts down and the highest temperature is lost.\n");
	return ok ? 0 : 1;
}
