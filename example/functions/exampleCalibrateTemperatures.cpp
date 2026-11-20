/*
 * PTAPI v2 - first and last temperature chosen by the API: calibrateTemperatures().
 *
 * The temperatures of the constructor are replaced by values measured on the
 * problem: a short descent gives a sample of worsening moves near a local
 * optimum, and the API solves for the temperatures at which a replica accepts
 * 'acceptCold' and 'acceptHot' of them (defaults 0.001 and 0.02). The ladder
 * is rebuilt with the distribution and the number of replicas of the
 * constructor.
 *
 *   g++ exampleCalibrateTemperatures.cpp -std=c++17 -O3 -pthread -o exampleCalibrateTemperatures
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	bool ok = true;

	// 1) defaults: the temperatures of the constructor are only placeholders
	PT<solTSP> algo(1.0f, 1.0f, 12, 200, 1000, 4, 0, 50);
	algo.setSeed(3);   // the calibration uses the seed too, so it is reproducible
	std::pair<double,double> t = algo.calibrateTemperatures(&prob);
	solTSP s = algo.start(0, &prob);
	printf("defaults                       Tmin %7.3f  Tmax %8.3f   best %.0f\n", t.first, t.second, s.evalSol);
	printRow("   ladder        ", algo.getInitialTemperatures());
	printRow("   accepted moves", algo.getMoveAcceptanceRates(), " %7.3f");
	ok = ok && t.first > 0.0 && t.first < t.second && prob.valid(s);

	// 2) other targets: (problem, samples, acceptHot, acceptCold)
	//    a hotter last temperature and a colder first one
	PT<solTSP> wide(1.0f, 1.0f, 12, 200, 1000, 4, 0, 50);
	wide.setSeed(3);
	t = wide.calibrateTemperatures(&prob, 20000, 0.10, 0.0003);
	s = wide.start(0, &prob);
	printf("acceptHot 0.10, cold 0.0003    Tmin %7.3f  Tmax %8.3f   best %.0f\n", t.first, t.second, s.evalSol);
	ok = ok && prob.valid(s);

	// 3) hottest replica accepting 23% of the worsening moves
	PT<solTSP> hot(1.0f, 1.0f, 12, 200, 1000, 4, 0, 50);
	hot.setSeed(3);
	t = hot.calibrateTemperatures(&prob, 5000, 0.23);
	s = hot.start(0, &prob);
	printf("acceptHot 0.23                 Tmin %7.3f  Tmax %8.3f   best %.0f\n", t.first, t.second, s.evalSol);
	ok = ok && prob.valid(s);

	// reference: temperatures chosen by hand for this problem
	PT<solTSP> hand(1.0f, 60.0f, 12, 200, 1000, 4, 0, 50);
	hand.setSeed(3);
	s = hand.start(0, &prob);
	printf("by hand                        Tmin %7.3f  Tmax %8.3f   best %.0f\n", 1.0, 60.0, s.evalSol);

	printf("\nOne run of each: to compare the options, average several seeds.\n");
	return ok ? 0 : 1;
}
