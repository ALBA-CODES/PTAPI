/*
 * PTAPI v2 - stop criteria: setTimeLimit(), setTargetValue(), setLowerBound(),
 * setPTLStopEnd(), setStagnationRounds(), setStagnationTime(), setStopCondition(),
 * setOnImprove() + requestStop(), and getStopReason().
 *
 * The PTL of the constructor is always a limit; the others are optional and
 * can be combined (the first one that is reached stops the run).
 * getStopReasonName() tells which one it was.
 *
 *   g++ exampleStop.cpp -std=c++17 -O3 -pthread -o exampleStop
 */
#include <thread>
#include <chrono>
#include "tsp.h"

static void report(const char* what, PT<solTSP>& algo, const solTSP& s){
	printf("%-30s best %6.0f   round %7d   %7.1f ms   stop: %s\n", what, s.evalSol, algo.getPTLEnd(), algo.getElapsedMs(), algo.getStopReasonName());
}

int main(){
	TSP prob(100);
	const int manyRounds = 2000000000;   // "no limit": another criterion stops the run
	bool ok = true;

	{	// 1) time limit in seconds, with fractions (v1: setTimeEnd(int), still accepted)
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setTimeLimit(0.25);
		solTSP s = algo.start(0, &prob);
		report("time limit of 0.25 s", algo, s);
		ok = ok && algo.getElapsedMs() < 1000.0;
	}

	{	// 2) target value: stops as soon as the best solution is <= the target
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setTargetValue(9000.0);
		algo.setTimeLimit(10.0);   // safety net
		solTSP s = algo.start(0, &prob);
		report("target value of 9000", algo, s);
		ok = ok && s.evalSol <= 9000.0;
	}

	{	// 3) stagnation (as in v1): stops after 10% of the PTL without a new best
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, 20000, 4, 0, 50);
		algo.setPTLStopEnd(0.10f);
		solTSP s = algo.start(0, &prob);
		report("10% of PTL without improving", algo, s);
		printf("%-30s last improvement at round %d\n", "", algo.getBestPTLSolIndex());
		ok = ok && algo.getPTLEnd() <= 20000;
	}

	{	// 4) callback on every new best solution; it may stop the run
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		int improvements = 0;
		algo.setOnImprove([&](const solTSP& sol, int ptl, double ms){
			// called from a worker thread, one call at a time
			++improvements;
			if(improvements % 50 == 0) printf("   improvement %4d: %6.0f at round %d (%.1f ms)\n", improvements, sol.evalSol, ptl, ms);
			if(improvements == 250) algo.requestStop();
		});
		algo.setTimeLimit(10.0);
		solTSP s = algo.start(0, &prob);
		report("stopped by the callback", algo, s);
		ok = ok && improvements >= 250 && algo.getElapsedMs() < 5000.0;
	}

	{	// 5) requestStop() from another thread (a signal handler, a GUI, ...)
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		std::thread other([&]{
			std::this_thread::sleep_for(std::chrono::milliseconds(150));
			algo.requestStop();
		});
		solTSP s = algo.start(0, &prob);
		other.join();
		report("stopped by another thread", algo, s);
		ok = ok && prob.valid(s);
	}

	{	// 6) rounds without a new best solution (an absolute number, not a fraction of PTL)
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setStagnationRounds(800);
		algo.setTimeLimit(10.0);
		solTSP s = algo.start(0, &prob);
		report("800 rounds without improving", algo, s);
		ok = ok && algo.getStopReason() == PT_STOP_STAGNATION && algo.getPTLEnd() > algo.getBestPTLSolIndex() + 800;
	}

	{	// 7) seconds without a new best solution
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setStagnationTime(0.15);
		algo.setTimeLimit(10.0);
		solTSP s = algo.start(0, &prob);
		report("0.15 s without improving", algo, s);
		ok = ok && algo.getStopReason() == PT_STOP_STAGNATION;
	}

	{	// 8) lower bound with a tolerance: stops when best <= bound + 5% of the bound
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setLowerBound(7300.0, 0.05);
		algo.setTimeLimit(10.0);
		solTSP s = algo.start(0, &prob);
		report("within 5% of the bound 7300", algo, s);
		ok = ok && algo.getStopReason() == PT_STOP_TARGET && s.evalSol <= 7300.0*1.05;
	}

	{	// 9) any other rule: a function called at the end of every round.
		//    Here: stop when the search spent 3 times as long without improving as
		//    it took to find the current best (and at least 50 ms have passed).
		PT<solTSP> algo(1.0f, 60.0f, 12, 200, manyRounds, 4, 0, 50);
		algo.setStopCondition([](const PTStatus& st){
			return st.elapsedMs > 50.0 && (st.elapsedMs - st.bestMs) > 3.0*st.bestMs;
		});
		algo.setTimeLimit(10.0);
		solTSP s = algo.start(0, &prob);
		report("user-defined condition", algo, s);
		ok = ok && algo.getStopReason() == PT_STOP_CONDITION;
	}

	return ok ? 0 : 1;
}
