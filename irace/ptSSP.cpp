/*
 * PTAPI v2 - program used by irace to tune the PT on the SSP example.
 *
 * All the parameters of the PT come from the command line through PTOptions
 * (../include/Options.h), so this file has nothing specific to irace: the same
 * pattern works for any problem. Only the two lines marked below depend on the
 * problem.
 *
 *   g++ ptSSP.cpp ../example/SSP.cpp -std=c++17 -O3 -pthread -o ptSSP
 *   ./ptSSP ../instances/SSP/F1001.txt --TEMP_INIT 0.2 --TEMP_FIM 1 --N_REPLICAS 16 --MCL 400 --PTL 400 --SEED 1
 *   ./ptSSP ../instances/SSP/F1001.txt --N_REPLICAS 16 --MCL 400 --TIME_LIMIT 1 --IRACE 1      (what irace calls)
 */
#include <iostream>
#include "../example/SSP.h"            // <-- the problem
#include "../include/Options.h"

int main(int argc, char* argv[]){
	try{
		PTOptions opt(argc, argv);
		if(opt.instance.empty()){ std::cerr << "usage: ptSSP instance [--OPTION value]...\n"; return 2; }

		SSP prob(opt.instance);        // <-- the problem
		solSSP best = opt.run<solSSP>(&prob);

		// with --IRACE 1|2 only the cost (and the time) goes to the standard output
		if(!opt.printForTuner()){
			const PTRunResult& r = opt.result;
			std::cout << "best " << best.evalSol << "  rounds " << r.rounds << "  seconds " << r.elapsedMs/1000.0 << "  seconds to best " << r.timeToBestMs/1000.0
			          << "  evaluations " << r.evaluations << "  stop " << PT<solSSP>::stopReasonName(r.stopReason) << "  seed " << r.seed << "\n";
		}
		return 0;
	}catch(const std::exception& e){
		std::cerr << "error: " << e.what() << "\n";
		return 1;
	}
}
