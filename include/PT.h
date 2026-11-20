#ifndef __PT_H__
#define __PT_H__

// ---------------------------------------------------------------------------
// PTAPI v2 - Parallel Tempering (dataflow implementation for CPU)
// ---------------------------------------------------------------------------

#include <vector>
#include <deque>
#include <thread>
#include <memory>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <cmath>
#include <chrono>
#include <random>
#include <algorithm>
#include <stdexcept>
#include <functional>
#include <exception>
#include <mutex>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <utility>
#include "Problem.h"
#include "Consumer.h"
#include "NodeMCMC.h"
#include "NodeSwap.h"

#define PTAPI_VERSION_MAJOR 2
#define PTAPI_VERSION_MINOR 0

// Exchange schemes (PT::setSwapScheme)
enum PTSwapScheme{
	PT_SWAP_SEQUENTIAL = 0,  // v1: every pair proposes on every round, from the first temperature up
	PT_SWAP_EVEN_ODD = 1     // pairs (0,1),(2,3),... and (1,2),(3,4),... on alternate rounds (non-reversible PT)
};

// Value of the 8th argument of the constructor (interval between temperature
// updates) that lets the API choose the interval: the first update happens after
// 10 exchange proposals and each interval is twice the previous one (10, 20, 40,
// ...), so the ladder is corrected often at the beginning, when it is far from
// the right one, and with better statistics later. See PT::setAutoUpdate.
const int PT_UPDATE_AUTO = 0;
// exchange proposals between two corrections of the ladder in the pilot runs of calibrateLadder
const int PT_CALIBRATION_INTERVAL = 10;

// Where the time of the last start() went (PT::getTiming). 'chainsMs',
// 'exchangesMs', 'updatesMs' and 'waitingMs' are summed over the threads and are
// only measured after PT::setProfiling(true).
struct PTTiming{
	double wallMs;          // whole start()
	double constructionMs;  // initial solutions (wall time)
	double chainsMs;        // Markov chains: the problem code (neighbor/evaluate) + acceptance test
	double exchangesMs;     // exchange proposals
	double updatesMs;       // dynamic temperature updates
	double waitingMs;       // threads waiting for a node to become ready + scheduling
	int threads;
};

// Use of one of the moves registered with Problem::setMoveProbabilities.
struct PTMoveStat{
	long long proposed;
	long long accepted;
};

template<typename S>
class PT{
	private:
		// PT variables
		float tempMin_; // init temp
		float tempMax_; // end temp
		float perPTLEnd = 0.0; // PT end PTL best sol
		int tempL_; // temp size
		int MKL_; // MC length
		int PTL_; // number of swaps
		int ut; // update type: 0 none, 1 23%, 2 equal rates, 3 feedback (thesis), 4 feedback (quantiles), 5 stochastic approximation, 6 equal rejection
		int tempUpdate; // number of swaps to update the temp
		int tempD_; // initial temperature distribution
		int indexBestPTL = 0; // Best PTL solution
		int timePTEnd = 0; // Time to finish the PT execution (s)
		int endPTL = 0; // End PTL
		std::deque<double> allTemps;
		S best;

		// v2
		double timeLimit = 0.0; // seconds (0 = no limit)
		bool hasTarget = false;
		double target = 0.0;
		bool hasSeed = false;
		unsigned long long seed = 0;
		unsigned long long usedSeed = 0;
		long long evaluations = 0;
		double elapsedMs = 0.0;
		double bestMs = 0.0;
		int threadsUsed = 0;
		std::vector<S> initialSols;
		std::vector<S> finalSols;
		std::vector<double> finalTemps;
		std::vector<double> moveRates;
		std::vector<double> swapRates;
		std::vector<PTImprovement> history;
		std::function<void(const S&, int, double)> onImprove;
		std::mutex runMtx;              // protects 'running'
		Consumer<S>* running = nullptr; // execution in progress (for requestStop)
		bool v1Compat = false;
		std::vector<PTMoveStat> moveStats;
		int swapScheme = PT_SWAP_SEQUENTIAL;
		long long roundTrips = 0;
		int roundsDone = 0;             // rounds completed by every replica in the last start()
		int stagRounds = 0;             // stop after this number of rounds without a new best (0 = off)
		double stagSeconds = 0.0;       // stop after this time without a new best (0 = off)
		std::function<bool(const PTStatus&)> stopCondition;
		int stopReason = PT_STOP_PTL;
		bool profiling = false;
		PTTiming timing{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0};
		bool parallelConstruction = false;
		int autoFirst = 10;             // automatic update interval: first interval and growth
		double autoGrowth = 2.0;

	public:
		typedef std::function<void(std::ostream&, const S&)> SolutionWriter;
		typedef std::function<void(std::istream&, S&)> SolutionReader;
	private:
		// checkpoint
		std::string ckptFile;
		double ckptEvery = 0.0;
		SolutionWriter ckptWriter;
		struct Resume{
			bool active = false;
			std::vector<double> temps;
			std::vector<S> sols;
			S best;
			int bestIndex = 0;
			double bestMs = 0.0;
			int rounds = 0;
			double elapsedMs = 0.0;
		} resume_;
		// pilot runs of calibrateLadder (update type 5 with another target and no upper limit)
		struct Pilot{
			bool active = false;
			double target = 0.234;
			long long updates = 0;          // updates made by the previous pilot runs
			std::vector<double> meanRho;    // result of the last one: mean log(log(T[i+1]/T[i])) of each pair
		} pilot_;

		std::deque<double> makeLadder(float tMin, float tMax, int tempL){
			switch (tempD_){
				case 1: return tempLinear(tMin,tMax,tempL);
				case 2: return tempILinear(tMin,tMax,tempL);
				case 3: return tempExp(tMin,tMax,tempL);
				case 4: return tempPG(tMin,tMax,tempL);
				default: throw std::invalid_argument("PTAPI: invalid initial temperature distribution (1..4)");
			}
		}

		static void writeCheckpointFile(const std::string& file, const SolutionWriter& writer, const std::vector<double>& temps,
		                                const std::vector<S>& sols, const S& bestSol, int bestIndex, double bestMsV, int rounds, double elapsed);

		static ptapi::Rng makeRng(unsigned long long s, unsigned long long stream){
			std::seed_seq sq{(std::uint32_t)(s & 0xffffffffu), (std::uint32_t)(s >> 32), (std::uint32_t)(stream & 0xffffffffu), (std::uint32_t)(stream >> 32)};
			return ptapi::Rng(sq);
		}
	public:
		PT(float tempMin, float tempMax, int tempL, int MKL, int PTL, int tempD, int upType, int tempUp); // constructor
		~PT();

		// ---- v1 interface ----
		template<typename P> S start(int thN, P* prob);
		S getBestSol();
		int getBestPTLSolIndex();
		int getPTLEnd();
		void setTimeEnd(int timeEnd);
		void setPTLStopEnd(float per);
		std::deque<double> tempPG(float tempMin, float tempMax, int tempL);
		std::deque<double> tempExp(float tempMin, float tempMax, int tempL);
		std::deque<double> tempLinear(float tempMin, float tempMax, int tempL);
		std::deque<double> tempILinear(float tempMin, float tempMax, int tempL);

		// ---- v2: configuration ----
		void setSeed(unsigned long long s){ hasSeed = true; seed = s; }          // reproducible runs
		void setTimeLimit(double seconds){ timeLimit = seconds; timePTEnd = 0; } // sub-second resolution
		void setTargetValue(double value){ hasTarget = true; target = value; }   // stop when best <= value
		void setTemperatures(const std::vector<double>& temps);                  // user-defined ladder
		void setInitialSolutions(const std::vector<S>& sols){ initialSols = sols; } // warm start (evaluated by the API)
		void setOnImprove(std::function<void(const S&, int, double)> f){ onImprove = std::move(f); } // (sol, PTL index, ms)
		// true = reproduce the v1 behaviours: energies truncated to int in the exchange
		// test, exchange rate of update type 1 divided by the total number of rounds
		// and the rules of update types 1 and 2 that do not keep the ladder anchored
		void setV1Compatibility(bool on){ v1Compat = on; }
		// Which pairs of temperatures propose an exchange on each round (see PTSwapScheme)
		void setSwapScheme(int scheme){
			if(scheme != PT_SWAP_SEQUENTIAL && scheme != PT_SWAP_EVEN_ODD) throw std::invalid_argument("PTAPI: invalid exchange scheme");
			swapScheme = scheme;
		}
		void requestStop(){ std::lock_guard<std::mutex> lock(runMtx); if(running) running->requestStop(); } // thread-safe

		// ---- v2: more stop criteria (all can be combined; the first one reached stops the run) ----
		// Stop after 'rounds' exchange rounds without a new best solution (0 = off).
		void setStagnationRounds(int rounds){ stagRounds = rounds; }
		// Stop after 'seconds' without a new best solution (0 = off).
		void setStagnationTime(double seconds){ stagSeconds = seconds; }
		// Stop when the best value is within 'gap' (relative, e.g. 0.01 = 1%) of a known
		// lower bound: best <= bound + gap*|bound|.
		void setLowerBound(double bound, double gap = 0.0){ setTargetValue(bound + gap*std::fabs(bound)); }
		// Any other rule: f(status) is called at the end of every round (from a worker
		// thread, one call at a time) and stops the run when it returns true.
		void setStopCondition(std::function<bool(const PTStatus&)> f){ stopCondition = std::move(f); }
		// Why the last start() ended (PTStopReason) and its name.
		int getStopReason() const { return stopReason; }
		static const char* stopReasonName(int reason){
			switch(reason){
				case PT_STOP_PTL: return "PTL";
				case PT_STOP_TIME: return "time";
				case PT_STOP_TARGET: return "target";
				case PT_STOP_STAGNATION: return "stagnation";
				case PT_STOP_REQUEST: return "request";
				case PT_STOP_CONDITION: return "condition";
				case PT_STOP_ERROR: return "error";
				default: return "?";
			}
		}
		const char* getStopReasonName() const { return stopReasonName(stopReason); }

		// ---- v2: time of each component ----
		// true: the time spent in the chains, the exchanges and the temperature updates
		// is measured (two clock readings per node execution, not per move).
		void setProfiling(bool on){ profiling = on; }
		const PTTiming& getTiming() const { return timing; }

		// ---- v2: initial solutions ----
		// true: the initial solutions are built by several threads at the same time
		// (worth it when construction() is an expensive heuristic). The construction
		// code must then be safe to call from several threads, as neighbor() already is.
		void setParallelConstruction(bool on){ parallelConstruction = on; }

		// ---- v2: automatic interval between temperature updates ----
		// Used when the 8th argument of the constructor is PT_UPDATE_AUTO (0): the first
		// update happens after 'firstInterval' proposals and each interval is 'growth'
		// times the previous one.
		void setAutoUpdate(int firstInterval, double growth = 2.0){
			if(firstInterval < 2 || !(growth >= 1.0)) throw std::invalid_argument("PTAPI: setAutoUpdate needs firstInterval >= 2 and growth >= 1");
			autoFirst = firstInterval; autoGrowth = growth;
		}

		// Chooses the first and the last temperature from a sample of moves of the
		// problem and rebuilds the ladder (same distribution and number of replicas).
		// Returns {Tmin, Tmax}. See docs/temperatures.md for the rule.
		template<typename P> std::pair<double,double> calibrateTemperatures(P* prob, int samples = 5000, double acceptHot = 0.02, double acceptCold = 0.001);

		// Builds the ladder from the first temperature with pilot runs of the PT: the
		// temperatures are corrected until adjacent temperatures exchange at 'swapRate'.
		// The last temperature is a result. tMin <= 0 keeps the first temperature
		// already defined. Returns the ladder.
		template<typename P> std::vector<double> calibrateLadder(P* prob, int thN = 0, double swapRate = 0.23, int rounds = 200, int iterations = 10, double tMin = 0.0);

		// ---- v2: checkpoint ----
		// The state (temperatures, solution of every replica, best solution, rounds and
		// time already used) is written to 'file' at the end of start() and, if
		// everySeconds > 0, also periodically while it runs. 'writer' writes one solution.
		void setCheckpoint(const std::string& file, SolutionWriter writer, double everySeconds = 0.0){
			ckptFile = file; ckptWriter = std::move(writer); ckptEvery = everySeconds;
		}
		// Writes the state of the last start().
		void saveCheckpoint(const std::string& file, SolutionWriter writer) const {
			writeCheckpointFile(file, writer, finalTemps, finalSols, best, indexBestPTL, bestMs, roundsDone, elapsedMs);
		}
		// Reads a checkpoint: the next start() continues that execution (remaining PTL
		// and remaining time). Returns false if the file does not exist.
		bool loadCheckpoint(const std::string& file, SolutionReader reader);

		// ---- v2: results of the last start() ----
		unsigned long long getSeed() const { return usedSeed; }
		long long getEvaluations() const { return evaluations; }
		double getElapsedMs() const { return elapsedMs; }
		double getTimeToBestMs() const { return bestMs; }
		int getThreadsUsed() const { return threadsUsed; }
		std::vector<double> getInitialTemperatures() const { return std::vector<double>(allTemps.begin(), allTemps.end()); }
		const std::vector<double>& getTemperatures() const { return finalTemps; }            // after the dynamic updates
		const std::vector<double>& getMoveAcceptanceRates() const { return moveRates; }      // one per temperature
		const std::vector<double>& getSwapAcceptanceRates() const { return swapRates; }      // one per adjacent pair
		const std::vector<S>& getReplicas() const { return finalSols; }                      // final solution of each temperature
		const std::vector<PTImprovement>& getConvergence() const { return history; }         // every new best solution
		const std::vector<PTMoveStat>& getMoveStatistics() const { return moveStats; }       // one per registered move
		// Replicas that travelled from the first temperature to the last and back: the
		// measure of how well the ladder carries solutions between its two ends.
		long long getRoundTrips() const { return roundTrips; }
		// Sum of the exchange rejection rates (global communication barrier of Syed
		// et al.): about half of the number of replicas that the interval needs.
		double getCommunicationBarrier() const { double b = 0.0; for(double r : swapRates) b += 1.0 - r; return b; }
};


// PT initialization: temperature values, MCMC length and number of exchange proposals
template<typename S>
PT<S>::PT(float tempMin, float tempMax, int tempL, int MKL, int PTL, int tempD, int upType, int tempUp)
	: tempMin_(tempMin)
	, tempMax_(tempMax)
	, tempL_(tempL)
	, MKL_(MKL)
	, PTL_(PTL)
	, ut(upType)
	, tempUpdate(tempUp)
	, tempD_(tempD)
{
	if(tempL < 1) throw std::invalid_argument("PTAPI: the number of replicas must be >= 1");
	if(upType < 0 || upType > 6) throw std::invalid_argument("PTAPI: invalid temperature update type (0..6)");

	allTemps = makeLadder(tempMin, tempMax, tempL);
}

template<typename S>
PT<S>::~PT(){
}

template<typename S>
void PT<S>::setTemperatures(const std::vector<double>& temps){
	if(temps.empty()) throw std::invalid_argument("PTAPI: empty temperature list");
	allTemps.assign(temps.begin(), temps.end());
	tempL_ = (int)temps.size();
}

template<typename S>
template<typename P>
S PT<S>::start(int thN, P* userProb){
	// Couples the evaluation function of the user class (either form) with the API
	Problem<S>* prob = userProb;
	prob->template ptapiBind<P>();

	// Continuation of a saved execution (loadCheckpoint); used only once
	Resume res;
	std::swap(res, resume_);
	if(res.active){
		allTemps.assign(res.temps.begin(), res.temps.end());
		tempL_ = (int)res.temps.size();
	}
	const int firstRound = res.active ? res.rounds : 0;
	const double elapsedBefore = res.active ? res.elapsedMs : 0.0;

	const int n = (int)allTemps.size();

	// There are never more than n MC nodes ready at the same time
	unsigned nTh = (thN > 0) ? (unsigned)thN : std::thread::hardware_concurrency();
	if(nTh == 0) nTh = 1;
	if(nTh > (unsigned)n) nTh = (unsigned)n;
	threadsUsed = (int)nTh;

	usedSeed = hasSeed ? seed : (((unsigned long long)std::random_device{}() << 32) ^ std::random_device{}());
	double limit = (timeLimit > 0.0) ? timeLimit : (double)timePTEnd;
	bool nothingLeft = (firstRound >= PTL_);
	if(limit > 0.0 && res.active){
		limit -= elapsedBefore/1000.0;          // remaining time
		if(limit <= 0.0){ limit = 0.0; nothingLeft = true; }
	}
	const unsigned long long streamBase = (unsigned long long)firstRound << 34;   // a resumed run does not repeat random numbers

	std::vector<std::unique_ptr<NodeMCMC<S>>> nMCMC;
	std::vector<std::unique_ptr<NodeSwap<S>>> nSwap;
	std::vector<std::unique_ptr<NodeTempUp<S>>> nUpTemp;
	std::exception_ptr failure;

	{
		Consumer<S> consumer(nTh, PTL_, perPTLEnd);
		if(hasTarget) consumer.setTarget(target);
		if(onImprove) consumer.setOnImprove(onImprove);
		if(stagRounds > 0) consumer.setStagnationRounds(stagRounds);
		if(stopCondition) consumer.setCondition(stopCondition);
		if(res.active) consumer.resume(res.best, res.bestIndex, res.bestMs, elapsedBefore);
		// From here on requestStop() works (also from the callback, while the initial
		// solutions are built). The guard also clears it if an exception is thrown.
		struct RunGuard{
			PT<S>* pt;
			RunGuard(PT<S>* pt_, Consumer<S>* c): pt(pt_){ std::lock_guard<std::mutex> lock(pt->runMtx); pt->running = c; }
			~RunGuard(){ std::lock_guard<std::mutex> lock(pt->runMtx); pt->running = nullptr; }
		} runGuard(this, &consumer);
		auto deadline = std::chrono::steady_clock::now() + std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double>(limit));

		// One MC node for each temperature
		for(int i = 0; i < n; ++i){
			nMCMC.emplace_back(new NodeMCMC<S>(MKL_, PTL_, allTemps[i], prob, &consumer, makeRng(usedSeed, streamBase + (unsigned long long)i), firstRound));
			nMCMC[i]->setProfiling(profiling);
		}

		// Initial solution of each replica: saved (resume), given by the user or built by the problem
		const double construction0 = consumer.elapsedMs();
		auto initReplica = [&](int i){
			if(res.active) nMCMC[i]->init(&res.sols[i], i, n);
			else if(initialSols.empty()) nMCMC[i]->init(nullptr, i, n);
			else{
				S given = initialSols[i % initialSols.size()];
				prob->ptapiEvaluate(given);
				nMCMC[i]->init(&given, i, n);
			}
		};
		if(!parallelConstruction || nTh < 2 || res.active){
			for(int i = 0; i < n; ++i) initReplica(i);
		}else{
			// each replica has its own random engine, so the solutions do not depend on the threads
			std::atomic<int> nextReplica{0};
			std::mutex errMtx;
			std::exception_ptr err;
			std::vector<std::thread> builders;
			for(unsigned t = 0; t < nTh; ++t) builders.emplace_back([&]{
				for(;;){
					const int i = nextReplica.fetch_add(1, std::memory_order_relaxed);
					if(i >= n) break;
					try{ initReplica(i); }
					catch(...){
						std::lock_guard<std::mutex> lock(errMtx);
						if(!err) err = std::current_exception();
						nextReplica.store(n, std::memory_order_relaxed);
					}
				}
			});
			for(std::thread& t : builders) t.join();
			if(err) std::rethrow_exception(err);
		}
		const double constructionMs = consumer.elapsedMs() - construction0;
		nMCMC.front()->setFirstTemp();
		nMCMC.back()->setLastTemp();

		// One exchange node (and one temperature update) for each pair of adjacent temperatures
		for(int i = 0; i + 1 < n; ++i){
			NodeTempUp<S>* tu = nullptr;
			if(ut){
				nUpTemp.emplace_back(new NodeTempUp<S>(nMCMC[i+1].get(), nMCMC[i].get(), i ? nUpTemp[i-1].get() : nullptr, ut, tempL_, !v1Compat,
					pilot_.active ? std::numeric_limits<double>::infinity() : *std::max_element(allTemps.begin(), allTemps.end())));
				tu = nUpTemp.back().get();
				if(pilot_.active) tu->setSA(pilot_.target, pilot_.updates);
				if(i) nUpTemp[i-1]->setUpTempUp(tu);
			}
			const bool autoUpdate = (tempUpdate == PT_UPDATE_AUTO);
			nSwap.emplace_back(new NodeSwap<S>(nMCMC[i+1].get(), nMCMC[i].get(), tu, &consumer, autoUpdate ? autoFirst : tempUpdate,
				makeRng(usedSeed, streamBase + 0x100000000ull + (unsigned long long)i), v1Compat, i, swapScheme == PT_SWAP_EVEN_ODD, autoUpdate ? autoGrowth : 1.0));
			nSwap.back()->setProfiling(profiling);
		}

		// Edges. The next exchange node comes first in the successor list: the
		// worker follows the (cheap) exchange chain and hands the MC nodes over.
		for(int i = 0; i + 1 < n; ++i){
			if(i + 2 < n) Node::addEdge(nSwap[i].get(), nSwap[i+1].get());
			Node::addEdge(nSwap[i].get(), nMCMC[i].get());
			Node::addEdge(nSwap[i].get(), nMCMC[i+1].get());
			Node::addEdge(nMCMC[i].get(), nSwap[i].get());
			Node::addEdge(nMCMC[i+1].get(), nSwap[i].get());
		}
		if(n == 1) Node::addEdge(nMCMC[0].get(), nMCMC[0].get()); // single replica: plain Metropolis

		// execution begins
		const bool periodic = !ckptFile.empty() && ckptWriter && (ckptEvery > 0.0);
		if(periodic) consumer.nextCheckpointEpoch();   // the nodes copy their state at the end of the round
		if(!nothingLeft) for(int i = 0; i < n; ++i) consumer.execAsync(nMCMC[i].get());
		consumer.getPool().seal();
		{
			// The calling thread sleeps until the graph drains; it wakes up for the time
			// limit, for the stop by time without improvement and for the periodic checkpoint.
			typedef std::chrono::steady_clock Clock;
			const bool hasDeadline = (limit > 0.0);
			const bool hasStagTime = (stagSeconds > 0.0);
			const Clock::duration every = std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(periodic ? ckptEvery : 0.0));
			Clock::time_point next = Clock::now() + every;
			for(;;){
				if(!hasDeadline && !hasStagTime && !periodic){ consumer.getPool().waitDone(); break; }
				Clock::time_point wake = hasDeadline ? deadline : Clock::time_point::max();
				if(periodic && next < wake) wake = next;
				if(hasStagTime){
					const double left = std::max(0.0005, stagSeconds - consumer.msSinceBest()/1000.0);
					const Clock::time_point t = Clock::now() + std::chrono::duration_cast<Clock::duration>(std::chrono::duration<double>(left));
					if(t < wake) wake = t;
				}
				if(consumer.getPool().waitUntil(wake)) break;
				const Clock::time_point now = Clock::now();
				if(hasDeadline && now >= deadline){
					consumer.requestStop(PT_STOP_TIME);
					consumer.getPool().waitDone();
					break;
				}
				if(hasStagTime && consumer.msSinceBest() >= stagSeconds*1000.0){
					consumer.requestStop(PT_STOP_STAGNATION);
					consumer.getPool().waitDone();
					break;
				}
				if(!periodic || now < next) continue;
				// state of every replica as copied by its node (at most one interval old)
				std::vector<S> sols(n);
				std::vector<double> temps(n);
				int rounds = PTL_;
				bool complete = true;
				for(int i = 0; i < n && complete; ++i){
					int r = 0;
					complete = nMCMC[i]->getSnapshot(sols[i], temps[i], r);
					rounds = std::min(rounds, r);
				}
				if(complete){
					S b; int bIndex = 0; double bMs = 0.0;
					consumer.getBest(b, bIndex, bMs);
					writeCheckpointFile(ckptFile, ckptWriter, temps, sols, b, bIndex, bMs, rounds, consumer.elapsedMs());
				}
				consumer.nextCheckpointEpoch();
				next = Clock::now() + every;
			}
		}

		failure = consumer.getError();
		stopReason = consumer.stopReason();
		timing.wallMs = consumer.elapsedMs() - elapsedBefore;
		timing.constructionMs = constructionMs;
		timing.threads = (int)nTh;
		best = consumer.getBestSol();
		indexBestPTL = consumer.getIndexPT();
		bestMs = consumer.getBestMs();
		history = consumer.getHistory();
		elapsedMs = consumer.elapsedMs();
	} // the threads are joined here

	// an exception thrown by the problem code inside a worker thread
	if(failure) std::rethrow_exception(failure);

	// statistics
	endPTL = 0;
	roundsDone = PTL_;
	evaluations = 0;
	moveStats.assign(prob->getMoveCount(), PTMoveStat{0, 0});
	finalSols.clear(); finalTemps.clear(); moveRates.clear(); swapRates.clear();
	for(auto& node : nMCMC){
		endPTL = std::max(endPTL, node->getRound());
		roundsDone = std::min(roundsDone, node->getRound());
		for(size_t k = 0; k < moveStats.size(); ++k){
			moveStats[k].proposed += node->getMoveProposed()[k];
			moveStats[k].accepted += node->getMoveAccepted()[k];
		}
		evaluations += node->getEvals();
		finalSols.push_back(node->getSol());
		finalTemps.push_back(node->getTemp());
		moveRates.push_back(node->getAccept());
	}
	long long chainsNs = 0, swapNs = 0, updateNs = 0;
	for(auto& node : nMCMC) chainsNs += node->getRunNs();
	for(auto& node : nSwap){
		swapRates.push_back(node->getAcceptRate());
		swapNs += node->getRunNs();
		updateNs += node->getUpdateNs();
	}
	if(pilot_.active){
		pilot_.meanRho.clear();
		for(auto& tu : nUpTemp) pilot_.meanRho.push_back(tu->hasRho() ? tu->getMeanRho() : std::numeric_limits<double>::quiet_NaN());
		if(!nUpTemp.empty()) pilot_.updates = nUpTemp.front()->getSAUpdates();
	}
	timing.chainsMs = (double)chainsNs/1.0e6;
	timing.updatesMs = (double)updateNs/1.0e6;
	timing.exchangesMs = (double)(swapNs - updateNs)/1.0e6;
	timing.waitingMs = profiling ? std::max(0.0, (double)timing.threads*(timing.wallMs - timing.constructionMs) - (double)(chainsNs + swapNs)/1.0e6) : 0.0;
	roundTrips = nMCMC.front()->getRoundTrips();

	if(!ckptFile.empty() && ckptWriter) saveCheckpoint(ckptFile, ckptWriter);

	return best;
}

// Temperatures from a sample of the moves of the problem.
// A descent (only non-worsening moves are taken) is run for 'samples' steps from
// a constructed solution. The worsening moves seen in its second half, when the
// solution is close to a local optimum, give the function
//     A(T) = average of exp(-delta/T)
// which is the fraction of those moves that a replica at temperature T accepts.
// Tmin and Tmax are the solutions of A(Tmin) = acceptCold and A(Tmax) = acceptHot.
template<typename S>
template<typename P>
std::pair<double,double> PT<S>::calibrateTemperatures(P* userProb, int samples, double acceptHot, double acceptCold){
	if(samples < 20) throw std::invalid_argument("PTAPI: calibrateTemperatures needs at least 20 samples");
	if(!(acceptCold > 0.0) || !(acceptHot < 1.0) || !(acceptCold < acceptHot)) throw std::invalid_argument("PTAPI: calibrateTemperatures needs 0 < acceptCold < acceptHot < 1");

	Problem<S>* prob = userProb;
	prob->template ptapiBind<P>();

	ptapi::Rng gen = makeRng(hasSeed ? seed : (((unsigned long long)std::random_device{}() << 32) ^ std::random_device{}()), 0xCA11B8A7E0000000ull);
	struct RngGuard{
		RngGuard(ptapi::Rng* g){ ptapi::detail::currentRng() = g; }
		~RngGuard(){ ptapi::detail::currentRng() = nullptr; }
	} guard(&gen);

	std::vector<double> early, ups;
	S cur = prob->construction();
	S cand = cur;
	for(int i = 0; i < samples; ++i){
		prob->neighborEval(cur, cand);
		const double d = cand.evalSol - cur.evalSol;
		if(d > 0.0){ if(std::isfinite(d)) ((i < samples/2) ? early : ups).push_back(d); }
		else std::swap(cur, cand);
	}
	if(ups.empty()) ups.swap(early);
	if(ups.empty()) throw std::runtime_error("PTAPI: calibrateTemperatures found no worsening move in the sample");

	const double dMin = *std::min_element(ups.begin(), ups.end());
	const double dMax = *std::max_element(ups.begin(), ups.end());
	auto solve = [&ups, dMin, dMax](double goal){
		double lo = std::log(dMin/700.0), hi = std::log(dMax*1.0e6);   // A(lo) ~ 0, A(hi) ~ 1
		for(int it = 0; it < 100; ++it){
			const double mid = 0.5*(lo + hi), t = std::exp(mid);
			double a = 0.0;
			for(double d : ups) a += std::exp(-d/t);
			if(a/(double)ups.size() < goal) lo = mid; else hi = mid;
		}
		return std::exp(0.5*(lo + hi));
	};
	const double tMin = solve(acceptCold);
	const double tMax = solve(acceptHot);

	tempMin_ = (float)tMin;
	tempMax_ = (float)tMax;
	allTemps = makeLadder(tempMin_, tempMax_, tempL_);
	return std::make_pair((double)tempMin_, (double)tempMax_);
}

// Ladder with a target exchange rate between adjacent temperatures.
// The PT itself is used as a pilot, with the stochastic approximation of update
// type 5 (Miasojedow, Moulines and Vihola, 2013) aimed at 'swapRate': the ladder
// is T[i+1] = T[i]*exp(exp(rho[i])) and, every few exchange proposals, each pair
// moves its rho by  step*(average exchange probability - swapRate), with a step
// that decreases along the whole calibration. 'iterations' pilot runs of
// 'rounds' rounds are made, each one continuing from the replicas and the ladder
// of the previous one. The result is the mean of the rhos over the second half
// of the pilot runs (the measured rates are noisy, above all at temperatures
// where the replicas are trapped; the mean filters that noise).
// The first temperature is kept; the others, including the last, are results.
// If there are more replicas than the target rate needs, the ladder stops at the
// first temperature that accepts 90% of the moves (a random walk: above it the
// exchange rate no longer falls) and the remaining replicas stay close to it.
template<typename S>
template<typename P>
std::vector<double> PT<S>::calibrateLadder(P* userProb, int thN, double swapRate, int rounds, int iterations, double tMin){
	if(!(swapRate > 0.0) || !(swapRate < 1.0)) throw std::invalid_argument("PTAPI: calibrateLadder needs 0 < swapRate < 1");
	if(rounds < 10 || iterations < 1) throw std::invalid_argument("PTAPI: calibrateLadder needs rounds >= 10 and iterations >= 1");

	const int n = (int)allTemps.size();
	const double first = (tMin > 0.0) ? tMin : std::min(allTemps.front(), allTemps.back());
	if(!(first > 0.0) || !std::isfinite(first)) throw std::invalid_argument("PTAPI: calibrateLadder needs a positive first temperature");

	// initial ladder: the current one, in ascending order, moved to the first temperature
	std::vector<double> ladder(allTemps.begin(), allTemps.end());
	std::sort(ladder.begin(), ladder.end());
	for(int k = n - 1; k >= 0; --k) ladder[k] *= first/ladder[0];
	for(int k = 1; k < n; ++k) if(!(ladder[k] > ladder[k-1]*1.0001)) ladder[k] = ladder[k-1]*1.0001;
	if(n < 2){ allTemps.assign(ladder.begin(), ladder.end()); tempMin_ = tempMax_ = (float)first; return ladder; }

	// The pilot runs are plain executions: no stop criteria, callbacks or checkpoint.
	// The configuration of the user is put back at the end (also after an exception).
	struct Saved{
		PT<S>* pt;
		int PTL, ut, tempUpdate, timePTEnd; float perPTLEnd; double timeLimit; bool hasTarget, hasSeed; unsigned long long seed;
		std::vector<S> initialSols; std::function<void(const S&, int, double)> onImprove; std::string ckptFile;
		int stagRounds; double stagSeconds; std::function<bool(const PTStatus&)> stopCondition; bool profiling;
		explicit Saved(PT<S>* p): pt(p), PTL(p->PTL_), ut(p->ut), tempUpdate(p->tempUpdate), timePTEnd(p->timePTEnd), perPTLEnd(p->perPTLEnd), timeLimit(p->timeLimit),
			hasTarget(p->hasTarget), hasSeed(p->hasSeed), seed(p->seed), initialSols(std::move(p->initialSols)), onImprove(std::move(p->onImprove)), ckptFile(std::move(p->ckptFile)),
			stagRounds(p->stagRounds), stagSeconds(p->stagSeconds), stopCondition(std::move(p->stopCondition)), profiling(p->profiling) {
			p->stagRounds = 0; p->stagSeconds = 0.0; p->stopCondition = nullptr; p->profiling = false;
			p->timePTEnd = 0; p->perPTLEnd = 0.0f; p->timeLimit = 0.0; p->hasTarget = false;
			p->initialSols.clear(); p->onImprove = nullptr; p->ckptFile.clear(); p->resume_ = Resume();
		}
		~Saved(){
			pt->PTL_ = PTL; pt->ut = ut; pt->tempUpdate = tempUpdate; pt->timePTEnd = timePTEnd; pt->perPTLEnd = perPTLEnd; pt->timeLimit = timeLimit;
			pt->hasTarget = hasTarget; pt->hasSeed = hasSeed; pt->seed = seed;
			pt->stagRounds = stagRounds; pt->stagSeconds = stagSeconds; pt->stopCondition = std::move(stopCondition); pt->profiling = profiling;
			pt->initialSols = std::move(initialSols); pt->onImprove = std::move(onImprove); pt->ckptFile = std::move(ckptFile);
			pt->resume_ = Resume(); pt->pilot_ = Pilot();
		}
	} saved(this);

	// Beyond the first temperature that is already a random walk the ratio is 1.01.
	const double closeGap = std::log(1.01);
	auto packHotEnd = [&](std::vector<double>& gaps){
		for(int k = 1; k < n; ++k){
			if(moveRates[k] < 0.9) continue;
			for(int j = k; j < n - 1; ++j) gaps[j] = std::min(gaps[j], closeGap);
			break;
		}
	};
	auto fromGaps = [&](const std::vector<double>& gaps){
		for(int k = 1; k < n; ++k){
			ladder[k] = ladder[k-1]*std::exp(gaps[k-1]);
			if(!(ladder[k] > ladder[k-1]*1.0001) || !std::isfinite(ladder[k])) ladder[k] = ladder[k-1]*1.0001;
		}
	};

	const unsigned long long baseSeed = saved.hasSeed ? saved.seed : (((unsigned long long)std::random_device{}() << 32) ^ std::random_device{}());
	PTL_ = rounds;
	ut = 5;
	tempUpdate = std::min(PT_CALIBRATION_INTERVAL, rounds);
	pilot_ = Pilot();
	pilot_.active = true;
	pilot_.target = swapRate;

	std::vector<double> sum(n - 1, 0.0), gaps(n - 1);
	int used = 0;
	for(int it = 0; it < iterations; ++it){
		allTemps.assign(ladder.begin(), ladder.end());
		hasSeed = true;
		seed = baseSeed + 0x9E3779B97F4A7C15ull*(unsigned long long)(it + 1);
		if(it > 0){
			// continues from the replicas of the previous pilot run
			resume_.active = true;
			resume_.temps = ladder;
			resume_.sols = finalSols;
			resume_.best = best;
			resume_.bestIndex = 0; resume_.bestMs = 0.0; resume_.rounds = 0; resume_.elapsedMs = 0.0;
		}
		start(thN, userProb);

		// ladder of the next pilot run: where this one ended
		for(int k = 0; k + 1 < n; ++k){
			const double g = std::log(finalTemps[k+1]/finalTemps[k]);
			gaps[k] = (std::isfinite(g) && g > 0.0) ? g : std::log(ladder[k+1]/ladder[k]);
		}
		packHotEnd(gaps);
		fromGaps(gaps);

		// second half of the pilot runs: mean of rho
		if(it >= iterations/2){
			for(int k = 0; k + 1 < n; ++k){
				const double r = pilot_.meanRho[k];
				sum[k] += std::isfinite(r) ? r : std::log(gaps[k]);
			}
			++used;
		}
	}
	for(int k = 0; k + 1 < n; ++k) gaps[k] = std::exp(sum[k]/(double)used);
	packHotEnd(gaps);
	fromGaps(gaps);

	allTemps.assign(ladder.begin(), ladder.end());
	tempMin_ = (float)ladder.front();
	tempMax_ = (float)ladder.back();
	return ladder;
}

template<typename S>
void PT<S>::writeCheckpointFile(const std::string& file, const SolutionWriter& writer, const std::vector<double>& temps,
                                const std::vector<S>& sols, const S& bestSol, int bestIndex, double bestMsV, int rounds, double elapsed){
	if(!writer) throw std::invalid_argument("PTAPI: checkpoint without a function to write the solution");
	if(temps.empty() || temps.size() != sols.size()) throw std::logic_error("PTAPI: there is no execution to save");
	const std::string tmp = file + ".tmp";
	{
		std::ofstream out(tmp.c_str(), std::ios::trunc);
		if(!out) throw std::runtime_error("PTAPI: cannot write the checkpoint file " + tmp);
		out.precision(17);
		out << "PTAPI-CHECKPOINT 1\n";
		out << temps.size() << ' ' << rounds << ' ' << bestIndex << ' ' << elapsed << ' ' << bestMsV << '\n';
		for(double t : temps) out << t << ' ';
		out << '\n';
		out << bestSol.evalSol << '\n';
		writer(out, bestSol);
		out << '\n';
		for(const S& s : sols){
			out << s.evalSol << '\n';
			writer(out, s);
			out << '\n';
		}
		out << "PTAPI-END\n";
		out.flush();
		if(!out) throw std::runtime_error("PTAPI: error writing the checkpoint file " + tmp);
	}
	// the previous checkpoint is replaced only when the new one is complete
	if(std::rename(tmp.c_str(), file.c_str()) != 0) throw std::runtime_error("PTAPI: cannot replace the checkpoint file " + file);
}

template<typename S>
bool PT<S>::loadCheckpoint(const std::string& file, SolutionReader reader){
	if(!reader) throw std::invalid_argument("PTAPI: checkpoint without a function to read the solution");
	std::ifstream in(file.c_str());
	if(!in) return false;

	auto bad = [&file]{ return std::runtime_error("PTAPI: invalid checkpoint file " + file); };
	auto readValue = [&in, &bad]() -> double {
		std::string tok;
		if(!(in >> tok)) throw bad();
		char* end = nullptr;
		double v = std::strtod(tok.c_str(), &end);
		if(end == tok.c_str() || *end != '\0') throw bad();
		in.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		return v;
	};
	auto readSol = [&](S& s){
		const double v = readValue();
		reader(in, s);
		if(!in) throw bad();
		s.evalSol = v;
	};

	std::string magic;
	int version = 0;
	long long count = 0;
	Resume r;
	if(!(in >> magic >> version) || magic != "PTAPI-CHECKPOINT" || version != 1) throw bad();
	if(!(in >> count >> r.rounds >> r.bestIndex >> r.elapsedMs >> r.bestMs) || count < 1 || count > 100000000 || r.rounds < 0) throw bad();
	r.temps.resize((size_t)count);
	for(double& t : r.temps) if(!(in >> t) || !(t > 0.0)) throw bad();
	readSol(r.best);
	r.sols.resize((size_t)count);
	for(S& s : r.sols) readSol(s);
	if(!(in >> magic) || magic != "PTAPI-END") throw bad();

	r.active = true;
	resume_ = std::move(r);
	return true;
}

template<typename S>
std::deque<double> PT<S>::tempPG(float tempMin, float tempMax, int tempL){
	std::deque<double> temps(tempL);
	if(tempL == 1){ temps[0] = tempMin; return temps; }
	double b = tempMax/tempMin;
	double ex = 1.0/((double)(tempL - 1));
	double r =  std::pow(b,ex);
	for(unsigned int ind = 0; ind < temps.size();++ind){
		temps[ind] = tempMin * std::pow(r,ind);
	}
	return temps;
}

template<typename S>
std::deque<double> PT<S>::tempExp(float tempMin, float tempMax, int tempL){
	std::deque<double>temps(tempL);
	if(tempL == 1){ temps[0] = tempMin; return temps; }
	double tempMaxExp = std::log(tempMax);
	double tempMinExp = std::log(tempMin);
	for(unsigned int ind = 0; ind < temps.size();++ind){
		double currentTemp = tempMinExp + ((tempMaxExp - tempMinExp) * ((double)ind/(double)(tempL - 1)));
		temps[ind] = std::exp(currentTemp);
	}
	return temps;
}

template<typename S>
std::deque<double> PT<S>::tempLinear(float tempMin, float tempMax, int tempL){
	std::deque<double> temps(tempL);
	if(tempL == 1){ temps[0] = tempMin; return temps; }
	for(unsigned int ind = 0; ind < temps.size();++ind){
		temps[ind] = tempMin + ((tempMax - tempMin) * ((double)ind/(double)(tempL - 1)));
	}
	return temps;
}

template<typename S>
std::deque<double> PT<S>::tempILinear(float tempMin, float tempMax, int tempL){
	std::deque<double> temps(tempL);
	if(tempL == 1){ temps[0] = tempMin; return temps; }
	double bm = 1.0/tempMax;
	double b1 = 1.0/tempMin;
	for(unsigned int ind = 0; ind < temps.size();++ind){
		double bi = bm + ((b1-bm) * ((double)ind/(double)(tempL - 1)));
		temps[ind] = 1.0 / bi;
	}
	return temps;
}

template<typename S>
S PT<S>::getBestSol(){
 return best;
}

template<typename S>
int PT<S>::getBestPTLSolIndex(){
 return indexBestPTL;
}

template<typename S>
int PT<S>::getPTLEnd(){
 return endPTL;
}

template<typename S>
void PT<S>::setTimeEnd(int timeEnd){
 timePTEnd = timeEnd;
 timeLimit = 0.0;
}

template<typename S>
void PT<S>::setPTLStopEnd(float per){
 perPTLEnd = per;
}

// v1 exported the std namespace from its headers and existing user code may
// rely on it. Define PTAPI_NO_USING_NAMESPACE_STD to avoid it.
#ifndef PTAPI_NO_USING_NAMESPACE_STD
using namespace std;
#endif

#endif
