/*
 * PTAPI v2 - self-contained example (random Euclidean TSP, 2-opt move).
 *
 * It shows the three steps of the API (solution, problem, PT) plus the v2
 * additions: ptapi::rng(), neighborEval() with delta evaluation, seed, stop
 * criteria, automatic calibration of the temperatures, checkpoint and the
 * statistics of the run.
 *
 *   g++ exampleTSP.cpp -std=c++17 -O3 -pthread -o exampleTSP
 *   ./exampleTSP [cities] [replicas] [MCL] [PTL] [threads] [seed] [fast 0/1] [tempUpdateType 0-6] [updateInterval]
 *
 * Environment variables: PT_TMIN, PT_TMAX (temperatures; default 1 and 60),
 * PT_CALIBRATE=1 (temperatures chosen by the API), PT_LADDER=0.23 (ladder
 * built for that exchange rate), PT_CHECKPOINT=file (resume
 * and save), PT_TIME (s), PT_TARGET, PT_STAG, PT_VERBOSE=1.
 */
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "../include/PT.h"
#include "../include/ExecTime.h"

// Step 1: solution encoding
struct solTSP: public solution{
	std::vector<int> sol;
};

// Step 2: problem dependent components
class TSP: public Problem<solTSP>{
	protected:
		int n;
		std::vector<double> d;
		double dist(int a, int b) const { return d[a*n + b]; }
	public:
		TSP(int n_, unsigned seed): n(n_), d(n_*n_){
			std::mt19937 g(seed);
			std::uniform_real_distribution<double> u(0, 1000);
			std::vector<double> x(n), y(n);
			for(int i = 0; i < n; i++){ x[i] = u(g); y[i] = u(g); }
			for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) d[i*n+j] = std::round(std::hypot(x[i]-x[j], y[i]-y[j]));
		}

		solTSP construction(){
			solTSP s;
			for(int i = 0; i < n; i++) s.sol.push_back(i);
			std::shuffle(s.sol.begin(), s.sol.end(), ptapi::rng()); // engine of the replica
			return evaluate(s);
		}

		solTSP neighbor(solTSP s){
			int a, b;
			do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);
			if(a > b) std::swap(a, b);
			std::reverse(s.sol.begin()+a, s.sol.begin()+b+1);
			return s;
		}

		solTSP evaluate(solTSP s){
			double c = 0;
			for(int i = 0; i < n; i++) c += dist(s.sol[i], s.sol[(i+1)%n]);
			s.evalSol = c;
			return s;
		}
};

// Optional (v2): the same problem with a copy-free step and delta evaluation
class TSPFast: public TSP{
	public:
		using TSP::TSP;

		void neighborEval(const solTSP& cur, solTSP& cand){
			int a, b;
			do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);
			if(a > b) std::swap(a, b);

			cand.sol = cur.sol; // reuses the memory of 'cand'
			std::reverse(cand.sol.begin()+a, cand.sol.begin()+b+1);

			if(a == 0 && b == n-1){ cand.evalSol = cur.evalSol; return; }
			int pa = cur.sol[(a+n-1)%n], ca = cur.sol[a], cb = cur.sol[b], nb = cur.sol[(b+1)%n];
			cand.evalSol = cur.evalSol - dist(pa,ca) - dist(cb,nb) + dist(pa,cb) + dist(ca,nb);
		}
};

int main(int argc, char* argv[]){
	int cities  = argc > 1 ? atoi(argv[1]) : 100;
	int rep     = argc > 2 ? atoi(argv[2]) : 16;
	int MCL     = argc > 3 ? atoi(argv[3]) : 500;
	int PTL     = argc > 4 ? atoi(argv[4]) : 1000;
	int threads = argc > 5 ? atoi(argv[5]) : 0;   // 0 = all cores
	long seed   = argc > 6 ? atol(argv[6]) : -1;
	int fast    = argc > 7 ? atoi(argv[7]) : 1;
	int upType  = argc > 8 ? atoi(argv[8]) : 0;   // 0 fixed, 1 23%, 2 equal rates, 3 feedback (thesis), 4 feedback (quantiles), 5 stochastic approximation, 6 equal rejection
	int upEvery = argc > 9 ? atoi(argv[9]) : 50;

	TSP* prob = fast ? new TSPFast(cities, 12345) : new TSP(cities, 12345);

	// Step 3: PT (same constructor of v1)
	float tMin = getenv("PT_TMIN") ? atof(getenv("PT_TMIN")) : 1.0f;
	float tMax = getenv("PT_TMAX") ? atof(getenv("PT_TMAX")) : 60.0f;
	PT<solTSP> algo(tMin, tMax, rep, MCL, PTL, 4, upType, upEvery);
	if(seed >= 0) algo.setSeed((unsigned long long)seed);
	if(getenv("PT_TIME"))   algo.setTimeLimit(atof(getenv("PT_TIME")));
	if(getenv("PT_TARGET")) algo.setTargetValue(atof(getenv("PT_TARGET")));
	if(getenv("PT_STAG"))   algo.setPTLStopEnd(atof(getenv("PT_STAG")));
	if(getenv("PT_SWAP"))   algo.setSwapScheme(atoi(getenv("PT_SWAP")));   // 0 sequential, 1 even-odd

	// temperatures chosen by the API from a sample of moves (instead of 1.0 and 60.0)
	if(getenv("PT_CALIBRATE")){
		std::pair<double,double> t = algo.calibrateTemperatures(prob);
		printf("calibrated temperatures: %.3f .. %.3f\n", t.first, t.second);
	}
	// every temperature searched for an exchange rate of 23% with the previous one
	// (after PT_CALIBRATE, it starts from the calibrated first temperature)
	if(getenv("PT_LADDER")){
		std::vector<double> l = algo.calibrateLadder(prob, threads, atof(getenv("PT_LADDER")) > 0 ? atof(getenv("PT_LADDER")) : 0.23);
		printf("calibrated ladder: %.3f .. %.3f\n", l.front(), l.back());
	}

	// checkpoint: continues the execution saved in the file, if it exists, and
	// saves the state every second and at the end
	if(getenv("PT_CHECKPOINT")){
		auto write = [](std::ostream& out, const solTSP& sol){ for(int c : sol.sol) out << c << ' '; };
		auto read  = [cities](std::istream& in, solTSP& sol){ sol.sol.resize(cities); for(int& c : sol.sol) in >> c; };
		if(algo.loadCheckpoint(getenv("PT_CHECKPOINT"), read)) printf("continuing from %s\n", getenv("PT_CHECKPOINT"));
		algo.setCheckpoint(getenv("PT_CHECKPOINT"), write, 1.0);
	}

	solTSP s = algo.start(threads, prob);

	// the returned solution must be a permutation with the reported cost
	std::vector<int> chk = s.sol; std::sort(chk.begin(), chk.end());
	bool perm = true; for(int i = 0; i < cities; i++) perm = perm && (chk[i] == i);
	bool cost = std::fabs(prob->evaluate(s).evalSol - s.evalSol) < 1e-6;

	printf("best=%.0f valid=%d bestPTL=%d endPTL=%d timeToBest=%.1fms elapsed=%.1fms evals=%lld threads=%d seed=%llu improvements=%zu\n",
		s.evalSol, (int)(perm && cost), algo.getBestPTLSolIndex(), algo.getPTLEnd(), algo.getTimeToBestMs(), algo.getElapsedMs(),
		algo.getEvaluations(), algo.getThreadsUsed(), algo.getSeed(), algo.getConvergence().size());
	if(getenv("PT_VERBOSE")){
		printf("temp      :"); for(double t : algo.getTemperatures()) printf(" %7.2f", t); printf("\n");
		printf("move rate :"); for(double r : algo.getMoveAcceptanceRates()) printf(" %7.3f", r); printf("\n");
		printf("swap rate :"); for(double r : algo.getSwapAcceptanceRates()) printf(" %7.3f", r); printf("\n");
		printf("round trips: %lld   communication barrier: %.2f\n", algo.getRoundTrips(), algo.getCommunicationBarrier());
	}
	delete prob;
	return (perm && cost) ? 0 : 1;
}
