#ifndef __PTAPI_OPTIONS_H__
#define __PTAPI_OPTIONS_H__

// ---------------------------------------------------------------------------
// PTAPI v2 - parameters from the command line (optional header)
//
// Every parameter of the PT can be given as "--NAME value", which is what a
// tuner such as irace needs (see the folder irace/). The names are those of the
// examples published with v1, plus the v2 additions. Options that the API does
// not know are kept for the problem (get()), so the parameters of the problem
// can be tuned in the same way.
//
//     int main(int argc, char* argv[]){
//         PTOptions opt(argc, argv);                 // opt.instance = first argument
//         MyProblem prob(opt.instance, opt.get("--ALPHA", 0.5));
//         mySolution best = opt.run<mySolution>(&prob);
//         if(!opt.printForTuner()) std::cout << best.evalSol << "\n";
//     }
// ---------------------------------------------------------------------------

#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstdio>
#include <stdexcept>
#include "PT.h"
#include "Experiment.h"

struct PTOptions{
	// constructor of PT
	float tempMin = 1.0f;          // --TEMP_INIT
	float tempMax = 10.0f;         // --TEMP_FIM
	int replicas = 8;              // --N_REPLICAS
	int MCL = 100;                 // --MCL
	int PTL = 1000;                // --PTL
	int tempDist = 4;              // --TEMP_DIST     1 linear, 2 inverse linear, 3 exponential, 4 geometric
	int updateType = 0;            // --TYPE_UPDATE   0..6
	int updateInterval = 50;       // --TEMP_UPDATE   proposals between updates; 0 = automatic
	// execution
	int threads = 0;               // --THREAD_USED   0 = all cores
	long long seed = -1;           // --SEED          < 0: drawn by the API
	double timeLimit = 0.0;        // --TIME_LIMIT    seconds
	bool hasTarget = false;        // --TARGET        stop when best <= value
	double target = 0.0;
	int stagnationRounds = 0;      // --STAGNATION_ROUNDS
	double stagnationTime = 0.0;   // --STAGNATION_TIME (seconds)
	float stagnationPTL = 0.0f;    // --STAGNATION_PTL  fraction of PTL (setPTLStopEnd)
	int swapScheme = 0;            // --SWAP_SCHEME   0 sequential, 1 even-odd
	int calibrate = 0;             // --CALIBRATE     1 first/last temperature, 2 ladder, 3 both
	bool parallelConstruction = false; // --PARALLEL_CONSTRUCTION 0|1
	bool profile = false;          // --PROFILE 0|1
	// tuner
	int tuner = 0;                 // --IRACE 1: print only the cost; 2: print "cost time"
	std::string tunerCost = "value";   // --IRACE_COST value | time | value+time
	double timeWeight = 0.0;       // --IRACE_TIME_WEIGHT  w in  value + w*seconds
	double penalty = 10.0;         // --IRACE_PENALTY      "time": a run that misses the target costs penalty*TIME_LIMIT

	std::string instance;                       // first argument that is not an option
	std::map<std::string, std::string> other;   // options the API does not know
	PTRunResult result;                         // filled by run()

	PTOptions(){}
	PTOptions(int argc, char* argv[]){ parse(argc, argv); }

	void parse(int argc, char* argv[]){
		std::vector<std::string> a(argv + 1, argv + argc);
		for(size_t i = 0; i < a.size(); ++i){
			const std::string& k = a[i];
			if(k.rfind("--", 0) != 0){
				if(instance.empty()) instance = k;
				else throw std::invalid_argument("PTAPI: unexpected argument " + k);
				continue;
			}
			if(i + 1 >= a.size()) throw std::invalid_argument("PTAPI: option " + k + " without a value");
			const std::string& v = a[++i];
			try{
				if(k == "--TEMP_INIT") tempMin = std::stof(v);
				else if(k == "--TEMP_FIM") tempMax = std::stof(v);
				else if(k == "--N_REPLICAS") replicas = std::stoi(v);
				else if(k == "--MCL") MCL = std::stoi(v);
				else if(k == "--PTL") PTL = std::stoi(v);
				else if(k == "--TEMP_DIST") tempDist = std::stoi(v);
				else if(k == "--TYPE_UPDATE") updateType = std::stoi(v);
				else if(k == "--TEMP_UPDATE") updateInterval = std::stoi(v);
				else if(k == "--THREAD_USED") threads = std::stoi(v);
				else if(k == "--SEED") seed = std::stoll(v);
				else if(k == "--TIME_LIMIT") timeLimit = std::stod(v);
				else if(k == "--TARGET"){ hasTarget = true; target = std::stod(v); }
				else if(k == "--STAGNATION_ROUNDS") stagnationRounds = std::stoi(v);
				else if(k == "--STAGNATION_TIME") stagnationTime = std::stod(v);
				else if(k == "--STAGNATION_PTL") stagnationPTL = std::stof(v);
				else if(k == "--SWAP_SCHEME") swapScheme = std::stoi(v);
				else if(k == "--CALIBRATE") calibrate = std::stoi(v);
				else if(k == "--PARALLEL_CONSTRUCTION") parallelConstruction = (std::stoi(v) != 0);
				else if(k == "--PROFILE") profile = (std::stoi(v) != 0);
				else if(k == "--IRACE") tuner = std::stoi(v);
				else if(k == "--IRACE_COST") tunerCost = v;
				else if(k == "--IRACE_TIME_WEIGHT") timeWeight = std::stod(v);
				else if(k == "--IRACE_PENALTY") penalty = std::stod(v);
				else other[k] = v;
			}catch(const std::invalid_argument&){
				throw std::invalid_argument("PTAPI: invalid value for " + k + ": " + v);
			}catch(const std::out_of_range&){
				throw std::invalid_argument("PTAPI: value out of range for " + k + ": " + v);
			}
		}
		if(tunerCost != "value" && tunerCost != "time" && tunerCost != "value+time") throw std::invalid_argument("PTAPI: --IRACE_COST must be value, time or value+time");
		if(tunerCost == "time" && (!hasTarget || !(timeLimit > 0.0))) throw std::invalid_argument("PTAPI: --IRACE_COST time needs --TARGET and --TIME_LIMIT");
	}

	// Options of the problem (those the API does not know).
	bool has(const std::string& name) const { return other.count(name) > 0; }
	std::string get(const std::string& name, const std::string& def) const { auto it = other.find(name); return it == other.end() ? def : it->second; }
	std::string get(const std::string& name, const char* def) const { return get(name, std::string(def)); }
	double get(const std::string& name, double def) const { auto it = other.find(name); return it == other.end() ? def : std::stod(it->second); }
	int get(const std::string& name, int def) const { auto it = other.find(name); return it == other.end() ? def : std::stoi(it->second); }

	// The PT object described by the options (not started).
	template<typename S>
	std::unique_ptr<PT<S>> make() const {
		std::unique_ptr<PT<S>> pt(new PT<S>(tempMin, tempMax, replicas, MCL, PTL, tempDist, updateType, updateInterval));
		apply(*pt);
		return pt;
	}

	// The same for an object built by the user.
	template<typename S>
	void apply(PT<S>& pt) const {
		if(seed >= 0) pt.setSeed((unsigned long long)seed);
		if(timeLimit > 0.0) pt.setTimeLimit(timeLimit);
		if(hasTarget) pt.setTargetValue(target);
		if(stagnationRounds > 0) pt.setStagnationRounds(stagnationRounds);
		if(stagnationTime > 0.0) pt.setStagnationTime(stagnationTime);
		if(stagnationPTL > 0.0f) pt.setPTLStopEnd(stagnationPTL);
		pt.setSwapScheme(swapScheme);
		pt.setParallelConstruction(parallelConstruction);
		pt.setProfiling(profile);
	}

	// Builds the PT, calibrates the temperatures if asked, runs it and keeps the
	// results in 'result'. Returns the best solution.
	template<typename S, typename P>
	S run(P* prob){
		std::unique_ptr<PT<S>> pt = make<S>();
		if(calibrate & 1) pt->calibrateTemperatures(prob);
		if(calibrate & 2) pt->calibrateLadder(prob, threads);
		S best = pt->start(threads, prob);
		result = ptResult(*pt, best);
		return best;
	}

	// Cost of the last run() for the tuner (smaller is better):
	//   value        the objective value;
	//   time         seconds to reach --TARGET; penalty*TIME_LIMIT if it was not reached;
	//   value+time   value + timeWeight*seconds.
	double tunerCostValue() const {
		const double seconds = result.elapsedMs/1000.0;
		if(tunerCost == "time") return (result.value <= target) ? result.timeToBestMs/1000.0 : penalty*timeLimit;
		if(tunerCost == "value+time") return result.value + timeWeight*seconds;
		return result.value;
	}

	// With --IRACE 1|2 prints the line the tuner reads ("cost" or "cost time") and
	// returns true: the program must then print nothing else on the standard output.
	bool printForTuner() const {
		if(tuner <= 0) return false;
		if(tuner == 1) std::printf("%.10g\n", tunerCostValue());
		else std::printf("%.10g %.6f\n", tunerCostValue(), result.elapsedMs/1000.0);
		return true;
	}
};

#endif
