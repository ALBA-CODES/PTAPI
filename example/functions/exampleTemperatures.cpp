/*
 * PTAPI v2 - temperature ladder: the four distributions of the constructor,
 * setTemperatures() and getInitialTemperatures().
 *
 *   g++ exampleTemperatures.cpp -std=c++17 -O3 -pthread -o exampleTemperatures
 */
#include "tsp.h"

int main(){
	TSP prob(100);
	bool ok = true;

	// 1) the distributions of the constructor (6th argument), from 1 to 60 with 8 replicas
	const char* names[] = {"", "1 linear        ", "2 inverse linear", "3 exponential   ", "4 geometric     "};
	for(int dist = 1; dist <= 4; dist++){
		PT<solTSP> algo(1.0f, 60.0f, 8, 200, 500, dist, 0, 50);
		printRow(names[dist], algo.getInitialTemperatures());
	}

	// 2) a ladder given by the user: any values, any number of replicas
	//    (the temperatures and the number of replicas of the constructor are replaced)
	PT<solTSP> algo(1.0f, 60.0f, 8, 200, 500, 4, 0, 50);
	algo.setTemperatures({1.0, 1.5, 2.5, 4.0, 7.0, 12.0, 20.0, 35.0, 60.0, 100.0});
	algo.setSeed(1);
	solTSP s = algo.start(0, &prob);

	printf("\nuser-defined ladder, %zu replicas, best %.0f\n", algo.getInitialTemperatures().size(), s.evalSol);
	printRow("temperature     ", algo.getTemperatures());
	printRow("accepted moves  ", algo.getMoveAcceptanceRates(), " %7.3f");
	printRow("accepted exch.  ", algo.getSwapAcceptanceRates(), " %7.3f");
	ok = ok && algo.getTemperatures().size() == 10 && prob.valid(s);

	// 3) one replica is accepted: plain Metropolis (simulated annealing at a fixed temperature)
	PT<solTSP> single(1.0f, 60.0f, 8, 200, 500, 4, 0, 50);
	single.setTemperatures({3.0});
	single.setSeed(1);
	solTSP m = single.start(0, &prob);
	printf("\none replica at T = 3: best %.0f (%d thread)\n", m.evalSol, single.getThreadsUsed());
	ok = ok && prob.valid(m);

	return ok ? 0 : 1;
}
