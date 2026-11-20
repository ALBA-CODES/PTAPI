#ifndef __NODEMCMC_H__
#define __NODEMCMC_H__

#include <cmath>
#include <random>
#include <atomic>
#include <utility>
#include <vector>
#include <mutex>
#include <chrono>
#include "Problem.h"
#include "Consumer.h"

// MC node: one homogeneous Markov chain (MCL steps) at one temperature.
template<typename S>
class alignas(64) NodeMCMC: public Node{
	private:
		int MCL;
		int PTL;
		int round = 0;                 // completed executions of this node
		std::atomic<double> temp;
		Problem<S>* prob;
		Consumer<S>* pool;
		S sol;
		S neigh;
		long long accept = 0;          // accepted moves (whole run)
		long long evals = 0;           // evaluated neighbors (whole run)
		int label = 0;                 // round-trip label of the replica: 0 none, 1 up, 2 down
		long long roundTrips = 0;
		std::atomic<long long> nUp{0};
		std::atomic<long long> nDown{0};
		bool firstTemp = false;
		bool lastTemp = false;
		ptapi::Rng gen;
		std::uniform_real_distribution<double> dis{0.0, 1.0};
		// statistics of the moves registered with Problem::setMoveProbabilities
		int nMoves = 0;
		std::vector<long long> moveProposed, moveAccepted;
		// copy of the state for the periodic checkpoint (see Consumer::checkpointEpoch)
		int seenEpoch = 0;
		std::mutex snapMtx;
		bool hasSnap = false;
		S snapSol;
		double snapTemp = 0.0;
		int snapRound = 0;
		// time spent in this node (only measured with PT::setProfiling)
		bool prof = false;
		long long runNs = 0;

		void snapshot(int epoch){
			std::lock_guard<std::mutex> lock(snapMtx);
			snapSol = sol;
			snapTemp = temp.load(std::memory_order_relaxed);
			snapRound = round;
			hasSnap = true;
			seenEpoch = epoch;
		}

	public:
		NodeMCMC(int MCL_, int PTL_, double temp_, Problem<S>* prob_, Consumer<S>* pool_, ptapi::Rng gen_, int firstRound = 0)
		: MCL(MCL_), PTL(PTL_), round(firstRound), temp(temp_), prob(prob_), pool(pool_), gen(gen_) {
			nMoves = prob->getMoveCount();
			moveProposed.assign(nMoves, 0);
			moveAccepted.assign(nMoves, 0);
		}

		// Initial solution: built by the problem, or given by the user (already evaluated).
		// 'index' and 'count': position of the replica in the ladder (see Problem::constructionReplica)
		void init(const S* given, int index = 0, int count = 1){
			ptapi::detail::currentRng() = &gen;
			struct Clear{ ~Clear(){ ptapi::detail::currentRng() = nullptr; } } clear;   // also if the problem code throws
			if(given) sol = *given; else sol = prob->constructionReplica(index, count);
			ptapi::detail::currentRng() = nullptr;
			neigh = sol;
			pool->setBestSol(sol, round);
		}

		bool run() override {
			if(!prof) return runChain();
			const std::chrono::steady_clock::time_point p0 = std::chrono::steady_clock::now();
			const bool again = runChain();
			runNs += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - p0).count();
			return again;
		}

		void setProfiling(bool on){ prof = on; }
		long long getRunNs() const { return runNs; }

		bool runChain(){
			if(pool->stopped()) return false;

			ptapi::detail::currentRng() = &gen;
			const double t = temp.load(std::memory_order_relaxed);
			int done = 0;

			try{
				if(nMoves == 0){
					while(done < MCL){
						// neighbor + evaluation (problem dependent)
						prob->neighborEval(sol, neigh);
						++done;

						double delta = neigh.evalSol - sol.evalSol;

						// Metropolis acceptance
						if(delta <= 0.0 || dis(gen) <= std::exp(-(delta/t))){
							std::swap(sol, neigh);
							++accept;
							if(delta < 0.0) pool->setBestSol(sol, round);
						}

						if(pool->stopped()) break;
					}
				}else{
					// the same loop, also counting which registered move was used
					int& lastMove = ptapi::detail::lastMove();
					while(done < MCL){
						lastMove = -1;
						prob->neighborEval(sol, neigh);
						++done;
						const int mv = lastMove;
						const bool counted = (mv >= 0) && (mv < nMoves);
						if(counted) ++moveProposed[mv];

						double delta = neigh.evalSol - sol.evalSol;

						if(delta <= 0.0 || dis(gen) <= std::exp(-(delta/t))){
							std::swap(sol, neigh);
							++accept;
							if(counted) ++moveAccepted[mv];
							if(delta < 0.0) pool->setBestSol(sol, round);
						}

						if(pool->stopped()) break;
					}
				}
			}catch(...){
				// exception thrown by the problem code
				ptapi::detail::currentRng() = nullptr;
				pool->fail(std::current_exception());
				return false;
			}

			ptapi::detail::currentRng() = nullptr;
			evals += done;
			++round;
			pool->checkStagnation(round);
			if(firstTemp && pool->hasCondition()){
				try{ pool->checkCondition(round); }
				catch(...){ pool->fail(std::current_exception()); return false; }
			}
			const int epoch = pool->checkpointEpoch();
			if(epoch != seenEpoch) snapshot(epoch);

			return (round < PTL) && !pool->stopped();
		}

		double getTemp() const { return temp.load(std::memory_order_relaxed); }
		void setTemp(double t){ temp.store(t, std::memory_order_relaxed); }
		double getEnergia() const { return sol.evalSol; }
		const S& getSol() const { return sol; }
		int getRound() const { return round; }
		long long getEvals() const { return evals; }
		double getAccept() const { return evals ? (double)accept/(double)evals : 0.0; }

		void setFirstTemp(){ firstTemp = true; }
		void setLastTemp(){ lastTemp = true; }

		// Exchange of replicas (solution + round-trip label) between two temperatures.
		void swapWith(NodeMCMC<S>& other){
			std::swap(sol, other.sol);
			std::swap(label, other.label);
		}

		// Load round-trip labels
		void trySetLabels(){
			if(firstTemp){
				if(label == 2) ++roundTrips;   // came back from the last temperature
				label = 1;
			}
			else if(lastTemp) label = 2;
		}
		// Replicas that went from the first temperature to the last and came back
		// (counted at the node of the first temperature).
		long long getRoundTrips() const { return roundTrips; }

		void updateFlow(){
			if(label == 1) nUp.fetch_add(1, std::memory_order_relaxed);
			else if(label == 2) nDown.fetch_add(1, std::memory_order_relaxed);
		}

		double getFlow() const {
			long long u = nUp.load(std::memory_order_relaxed), d = nDown.load(std::memory_order_relaxed);
			if((u+d) == 0) return 0;
			return (double)u/(double)(u+d);
		}

		// true if a labelled replica has visited this temperature since the last reset
		bool hasFlow() const {
			return (nUp.load(std::memory_order_relaxed) + nDown.load(std::memory_order_relaxed)) > 0;
		}

		const std::vector<long long>& getMoveProposed() const { return moveProposed; }
		const std::vector<long long>& getMoveAccepted() const { return moveAccepted; }

		// Last copy taken for the checkpoint (false if there is none yet).
		bool getSnapshot(S& s, double& t, int& r){
			std::lock_guard<std::mutex> lock(snapMtx);
			if(!hasSnap) return false;
			s = snapSol; t = snapTemp; r = snapRound;
			return true;
		}

		void resetFlow(){
			nUp.store(0, std::memory_order_relaxed);
			nDown.store(0, std::memory_order_relaxed);
		}
};

#endif
