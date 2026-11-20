/*
 * PTAPI v2 - checkpoint and resume: setCheckpoint(), loadCheckpoint() and
 * saveCheckpoint().
 *
 * The API saves the temperatures, the solution of every replica, the best
 * solution and the rounds and time already used. The user gives two functions:
 * one that writes a solution and one that reads it back (the objective value
 * is handled by the API).
 *
 * Run it twice to see the usual pattern of a long job:
 *   ./exampleCheckpoint          runs for 0.15 s and is "interrupted"
 *   ./exampleCheckpoint          continues from the file until the end
 *   ./exampleCheckpoint demo     both parts in one execution (default of 'make run')
 *
 *   g++ exampleCheckpoint.cpp -std=c++17 -O3 -pthread -o exampleCheckpoint
 */
#include <cstring>
#include "tsp.h"

static const char* FILE_NAME = "exampleCheckpoint.ckpt";
static const int CITIES = 200, PTL = 3000;

// how one solution goes to the file and comes back
static void writeSol(std::ostream& out, const solTSP& s){ for(int c : s.sol) out << c << ' '; }
static void readSol(std::istream& in, solTSP& s){ s.sol.resize(CITIES); for(int& c : s.sol) in >> c; }

// One part of the job. 'interruptAfter' > 0 simulates the machine going down.
static solTSP part(TSP& prob, double interruptAfter){
	PT<solTSP> algo(1.0f, 60.0f, 12, 300, PTL, 4, 0, 50);

	// continues the saved execution, if there is one
	bool resumed = algo.loadCheckpoint(FILE_NAME, readSol);

	// saves every 0.1 s while it runs, and at the end
	algo.setCheckpoint(FILE_NAME, writeSol, 0.1);

	// a time limit is used here only to interrupt the run in the middle; a real
	// interruption (kill, power failure) leaves the last periodic copy in the file
	if(interruptAfter > 0.0) algo.setTimeLimit(interruptAfter);

	solTSP s = algo.start(0, &prob);
	printf("%-9s rounds %4d of %d   best %6.0f   time of the job so far %7.1f ms   evaluations of this part %lld\n",
		resumed ? "resumed:" : "new run:", algo.getPTLEnd(), PTL, s.evalSol, algo.getElapsedMs(), algo.getEvaluations());
	return s;
}

int main(int argc, char* argv[]){
	TSP prob(CITIES);
	bool ok = true;

	if(argc > 1 && strcmp(argv[1], "demo") == 0){
		std::remove(FILE_NAME);
		solTSP a = part(prob, 0.15);   // interrupted
		solTSP b = part(prob, 0.0);    // continues until the PTL of the constructor
		ok = prob.valid(b) && b.evalSol <= a.evalSol;

		// saveCheckpoint(): the state of the last start(), written on demand
		PT<solTSP> algo(1.0f, 60.0f, 12, 300, 200, 4, 0, 50);
		algo.setSeed(1);
		algo.start(0, &prob);
		algo.saveCheckpoint(FILE_NAME, writeSol);
		PT<solTSP> longer(1.0f, 60.0f, 12, 300, 400, 4, 0, 50);   // the same job with a larger PTL
		ok = ok && longer.loadCheckpoint(FILE_NAME, readSol);
		solTSP c = longer.start(0, &prob);
		printf("saved after 200 rounds and extended to 400: best %.0f, rounds %d\n", c.evalSol, longer.getPTLEnd());
		ok = ok && longer.getPTLEnd() == 400;
		std::remove(FILE_NAME);
		return ok ? 0 : 1;
	}

	// first execution: interrupted; next executions: continue from the file
	std::ifstream exists(FILE_NAME);
	if(!exists.good()){
		part(prob, 0.15);
		printf("state saved in %s; run the program again to continue\n", FILE_NAME);
	}else{
		exists.close();
		solTSP s = part(prob, 0.0);
		ok = prob.valid(s);
		std::remove(FILE_NAME);
		printf("finished; %s removed\n", FILE_NAME);
	}
	return ok ? 0 : 1;
}
