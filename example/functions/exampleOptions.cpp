/*
 * PTAPI v2 - parameters from the command line: PTOptions (header
 * ../../include/Options.h). It is what a tuner such as irace needs; the folder
 * ../../irace has a complete scenario.
 *
 * Every parameter of the PT is an option "--NAME value". Options the API does
 * not know are kept for the problem, so its parameters can be tuned too.
 *
 *   g++ exampleOptions.cpp -std=c++17 -O3 -pthread -o exampleOptions
 *   ./exampleOptions 100 --TEMP_INIT 1 --TEMP_FIM 60 --N_REPLICAS 12 --MCL 200 --PTL 500 --SEED 3
 *   ./exampleOptions 100 --N_REPLICAS 12 --MCL 200 --TIME_LIMIT 0.2 --CALIBRATE 1 --IRACE 1
 *
 * Options: --TEMP_INIT --TEMP_FIM --N_REPLICAS --MCL --PTL --TEMP_DIST --TYPE_UPDATE
 * --TEMP_UPDATE --THREAD_USED --SEED --TIME_LIMIT --TARGET --STAGNATION_ROUNDS
 * --STAGNATION_TIME --STAGNATION_PTL --SWAP_SCHEME --CALIBRATE --PARALLEL_CONSTRUCTION
 * --PROFILE --IRACE --IRACE_COST --IRACE_TIME_WEIGHT --IRACE_PENALTY
 */
#include "tsp.h"
#include "../../include/Options.h"

int main(int argc, char* argv[]){
	// without arguments (make run): the example supplies a command line
	const char* demo[] = {"exampleOptions", "100", "--TEMP_INIT", "1", "--TEMP_FIM", "60", "--N_REPLICAS", "12", "--MCL", "200",
	                      "--PTL", "500", "--TYPE_UPDATE", "5", "--TEMP_UPDATE", "0", "--SEED", "3", "--INSTANCE_SEED", "777"};
	if(argc < 2){ argc = (int)(sizeof(demo)/sizeof(demo[0])); argv = const_cast<char**>(demo); }

	try{
		PTOptions opt(argc, argv);

		// the first argument is the instance; --INSTANCE_SEED is an option of the problem
		TSP prob(std::stoi(opt.instance), (unsigned)opt.get("--INSTANCE_SEED", 12345));

		solTSP best = opt.run<solTSP>(&prob);   // builds the PT, calibrates if asked, runs

		// with --IRACE 1|2 the program prints only what the tuner reads
		if(!opt.printForTuner()){
			const PTRunResult& r = opt.result;
			printf("%d replicas, MCL %d, temperatures %.2f to %.2f, update type %d\n", opt.replicas, opt.MCL, opt.tempMin, opt.tempMax, opt.updateType);
			printf("best %.0f   rounds %d   %.1f ms   stop: %s   seed %llu\n", best.evalSol, r.rounds, r.elapsedMs,
				PT<solTSP>::stopReasonName(r.stopReason), r.seed);
		}
		return prob.valid(best) ? 0 : 1;
	}catch(const std::exception& e){
		fprintf(stderr, "error: %s\n", e.what());
		return 1;
	}
}
