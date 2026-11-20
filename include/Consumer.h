#ifndef __CONSUMER_H__
#define __CONSUMER_H__

#include <atomic>
#include <chrono>
#include <mutex>
#include <vector>
#include <limits>
#include <functional>
#include <exception>
#include "Pool.h"

// One entry of the convergence history (a new global best solution).
struct PTImprovement{
	double timeMs;  // elapsed time since start()
	int ptl;        // exchange round (PTL index) in which it was found
	double value;   // objective value
};

// Why the last start() ended (PT::getStopReason).
enum PTStopReason{
	PT_STOP_PTL = 0,         // every replica completed the PTL rounds
	PT_STOP_TIME = 1,        // time limit
	PT_STOP_TARGET = 2,      // target value / lower bound reached
	PT_STOP_STAGNATION = 3,  // rounds or seconds without a new best solution
	PT_STOP_REQUEST = 4,     // requestStop() called by the user
	PT_STOP_CONDITION = 5,   // the function given to setStopCondition returned true
	PT_STOP_ERROR = 6        // exception thrown by the problem code
};

// What the function given to PT::setStopCondition receives at the end of each round.
struct PTStatus{
	int round;          // exchange rounds completed
	double elapsedMs;   // time since start()
	double bestValue;   // best objective value so far
	int bestRound;      // round in which it was found
	double bestMs;      // time at which it was found
};

// State shared by all the nodes of one execution: the thread pool, the stop
// flag and the global best solution.
template<typename S>
class Consumer{
	private:
		Pool pool;
		// read on every Markov chain step: kept away from the data the pool writes
		alignas(64) std::atomic<bool> stop{false};
		std::atomic<int> reason{PT_STOP_PTL};
		alignas(64) std::atomic<double> bestV{std::numeric_limits<double>::infinity()};
		std::atomic<int> indexPT{0};
		std::atomic<int> ckptEpoch{0};
		std::atomic<double> lastImproveMs{0.0};   // for the stop by seconds without improvement
		int stagRounds = 0;
		std::function<bool(const PTStatus&)> condition;
		alignas(64) std::mutex mtxSol;
		S bestSol;
		double bestMs = 0.0;
		std::vector<PTImprovement> history;
		std::function<void(const S&, int, double)> onImprove;
		bool hasTarget = false;
		double target = 0.0;
		int PTL = 0;
		float maxPercPTL = 0.0;
		std::chrono::steady_clock::time_point t0;
		std::exception_ptr error;
	public:
		Consumer(unsigned nTh, int PTL_, float maxPerc)
		: pool(nTh), PTL(PTL_), maxPercPTL(maxPerc), t0(std::chrono::steady_clock::now()) {}

		Pool& getPool(){ return pool; }

		// Periodic checkpoint: each MC node copies its state when it sees a new epoch.
		int checkpointEpoch() const { return ckptEpoch.load(std::memory_order_relaxed); }
		void nextCheckpointEpoch(){ ckptEpoch.fetch_add(1, std::memory_order_relaxed); }

		// A resumed execution continues the clock and the best solution of the saved one.
		void resume(const S& b, int index, double ms, double elapsedBeforeMs){
			std::lock_guard<std::mutex> lock(mtxSol);
			t0 -= std::chrono::duration_cast<std::chrono::steady_clock::duration>(std::chrono::duration<double, std::milli>(elapsedBeforeMs));
			bestSol = b;
			bestMs = ms;
			indexPT.store(index, std::memory_order_relaxed);
			bestV.store(b.evalSol, std::memory_order_relaxed);
			lastImproveMs.store(elapsedBeforeMs, std::memory_order_relaxed);   // the stagnation clock starts again
		}
		void getBest(S& b, int& index, double& ms){
			std::lock_guard<std::mutex> lock(mtxSol);
			b = bestSol; index = indexPT.load(std::memory_order_relaxed); ms = bestMs;
		}
		void execAsync(Node* item){ pool.push(item); }

		bool stopped() const { return stop.load(std::memory_order_relaxed); }
		// The first reason given is the one reported.
		void requestStop(int why = PT_STOP_REQUEST){
			int none = PT_STOP_PTL;
			reason.compare_exchange_strong(none, why, std::memory_order_relaxed);
			stop.store(true, std::memory_order_relaxed);
		}
		int stopReason() const { return reason.load(std::memory_order_relaxed); }

		// An exception left the problem code: stop everything; PT::start() rethrows it.
		void fail(std::exception_ptr e){
			{
				std::lock_guard<std::mutex> lock(mtxSol);
				if(!error) error = e;
			}
			requestStop(PT_STOP_ERROR);
		}
		std::exception_ptr getError(){ std::lock_guard<std::mutex> lock(mtxSol); return error; }

		void setTarget(double v){ hasTarget = true; target = v; }
		void setOnImprove(std::function<void(const S&, int, double)> f){ onImprove = std::move(f); }
		void setStagnationRounds(int r){ stagRounds = r; }
		void setCondition(std::function<bool(const PTStatus&)> f){ condition = std::move(f); }
		bool hasCondition() const { return (bool)condition; }
		double msSinceBest() const { return elapsedMs() - lastImproveMs.load(std::memory_order_relaxed); }

		// User-defined stop condition; called once per round (by the node of the first temperature).
		void checkCondition(int round){
			PTStatus st;
			st.round = round;
			st.elapsedMs = elapsedMs();
			{
				std::lock_guard<std::mutex> lock(mtxSol);
				st.bestValue = bestV.load(std::memory_order_relaxed);
				st.bestRound = indexPT.load(std::memory_order_relaxed);
				st.bestMs = bestMs;
			}
			if(condition(st)) requestStop(PT_STOP_CONDITION);
		}

		double elapsedMs() const {
			return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
		}

		// Called by a replica whenever it improves its current solution.
		// Lock-free unless the solution really is a new global best.
		void setBestSol(const S& b, int i){
			if(!(b.evalSol < bestV.load(std::memory_order_relaxed))) return;
			std::lock_guard<std::mutex> lock(mtxSol);
			if(!(b.evalSol < bestV.load(std::memory_order_relaxed))) return;
			bestSol = b;
			bestMs = elapsedMs();
			indexPT.store(i, std::memory_order_relaxed);
			bestV.store(b.evalSol, std::memory_order_relaxed);
			lastImproveMs.store(bestMs, std::memory_order_relaxed);
			history.push_back({bestMs, i, b.evalSol});
			if(onImprove) onImprove(bestSol, i, bestMs);
			if(hasTarget && b.evalSol <= target) requestStop(PT_STOP_TARGET);
		}

		// Stop criterion: no new best solution in the last (maxPercPTL * PTL) rounds,
		// or in the last 'stagRounds' rounds.
		void checkStagnation(int round){
			if((maxPercPTL > 0.0) && ((indexPT.load(std::memory_order_relaxed) + (PTL * maxPercPTL)) < round)) requestStop(PT_STOP_STAGNATION);
			if((stagRounds > 0) && ((indexPT.load(std::memory_order_relaxed) + stagRounds) < round)) requestStop(PT_STOP_STAGNATION);
		}

		S getBestSol(){ std::lock_guard<std::mutex> lock(mtxSol); return bestSol; }
		int getIndexPT() const { return indexPT.load(); }
		double getBestMs() const { return bestMs; }
		const std::vector<PTImprovement>& getHistory() const { return history; }
};

#endif
