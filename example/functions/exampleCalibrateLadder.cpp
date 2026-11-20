/*
 * PTAPI v2 - ladder built for an exchange rate: calibrateLadder().
 *
 * The first temperature is kept and all the others are searched with pilot
 * runs of the PT itself, in which the distance between adjacent temperatures
 * is corrected, with a decreasing step, until every pair exchanges at the
 * given rate (default 23%). The last temperature is a result: where the number
 * of replicas reaches with that rate.
 *
 * The rate is met where it is a property of the ladder: at temperatures in
 * which the replicas keep moving (case 2). At temperatures so low that the
 * replicas are trapped in local optima, the rate of a pair depends on the two
 * solutions it holds and changes from run to run with the same ladder; there
 * the calibration can only be right on average (case 1).
 *
 *   calibrateLadder(problem, threads, swapRate, rounds, iterations, tMin)
 *
 *   g++ exampleCalibrateLadder.cpp -std=c++17 -O3 -pthread -o exampleCalibrateLadder
 */
#include "tsp.h"

static void show(const char* title, PT<solTSP>& algo, const solTSP& s){
	double mean = 0.0;
	for(double r : algo.getSwapAcceptanceRates()) mean += r/algo.getSwapAcceptanceRates().size();
	printf("%s: best %.0f, round trips %lld, mean exchange rate %.2f\n", title, s.evalSol, algo.getRoundTrips(), mean);
	printRow("   ladder        ", algo.getInitialTemperatures());
	printRow("   exchange rate ", algo.getSwapAcceptanceRates(), " %7.3f");
}

int main(){
	TSP prob(100);
	const int replicas = 10, MCL = 200, PTL = 1000;
	bool ok = true;

	// reference: geometric ladder from 1 to 60
	PT<solTSP> hand(1.0f, 60.0f, replicas, MCL, PTL, 4, 0, 50);
	hand.setSeed(9);
	solTSP s = hand.start(0, &prob);
	show("geometric ladder from 1 to 60", hand, s);

	// 1) ladder for 23% of exchanges, starting at the first temperature of the constructor
	PT<solTSP> algo(1.0f, 60.0f, replicas, MCL, PTL, 4, 0, 50);
	algo.setSeed(9);
	std::vector<double> ladder = algo.calibrateLadder(&prob);
	s = algo.start(0, &prob);
	show("\ncalibrateLadder, 23%", algo, s);
	ok = ok && (int)ladder.size() == replicas && std::fabs(ladder.front() - 1.0) < 1e-6 && prob.valid(s);

	// 2) all the parameters: threads (0 = all cores), target rate, rounds of each
	//    pilot run, number of pilot runs and first temperature. From T = 20 the
	//    replicas are not trapped and every pair stays close to the target.
	PT<solTSP> longer(1.0f, 60.0f, replicas, MCL, 5*PTL, 4, 0, 50);
	longer.setSeed(9);
	longer.calibrateLadder(&prob, 0, 0.23, 200, 10, 20.0);
	s = longer.start(0, &prob);
	show("\ncalibrateLadder, 23%, 10 pilot runs of 200 rounds, from T = 20", longer, s);
	double worst = 0.0;
	for(double r : longer.getSwapAcceptanceRates()) worst = std::max(worst, std::fabs(r - 0.23));
	printf("   largest deviation of a pair from 23%%: %.3f\n", worst);
	ok = ok && std::fabs(longer.getInitialTemperatures().front() - 20.0) < 1e-6 && prob.valid(s) && worst < 0.10;

	// 3) the two calibrations together: the first temperature from the sample of
	//    moves, the others from the exchange rate
	PT<solTSP> both(1.0f, 60.0f, replicas, MCL, PTL, 4, 0, 50);
	both.setSeed(9);
	both.calibrateTemperatures(&prob);
	both.calibrateLadder(&prob);
	s = both.start(0, &prob);
	show("\ncalibrateTemperatures + calibrateLadder", both, s);
	ok = ok && prob.valid(s);

	printf("\nThe pilot runs cost iterations x rounds x replicas x MCL evaluations (%d with\n"
	       "the defaults here), which are not part of the run that follows.\n", 10*200*replicas*MCL);
	return ok ? 0 : 1;
}
