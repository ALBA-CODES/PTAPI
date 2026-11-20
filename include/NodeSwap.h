#ifndef __NODESWAP_H__
#define __NODESWAP_H__

#include <cmath>
#include <random>
#include <chrono>
#include <algorithm>
#include "NodeTempUp.h"
#include "ReplicaEst.h"

// SW node: exchange proposal between two adjacent temperatures.
template<typename S>
class NodeSwap: public Node{
	private:
		NodeMCMC<S>* left;
		NodeMCMC<S>* right;
		NodeTempUp<S>* tempUp;
		Consumer<S>* pool;
		int round = 0;
		long long attempts = 0;      // exchange proposals (= round, except in the even-odd scheme)
		int index = 0;               // position of the pair in the ladder
		bool evenOdd = false;        // true: the pair only proposes on alternate rounds
		int accept = 0;              // accepted since the last temperature update
		long long acceptTotal = 0;
		// Temperature updates: the next one happens when 'attempts' reaches 'nextUpdate';
		// 'window' is the number of proposals since the previous one. With growth > 1
		// the window grows after each update (automatic interval).
		long long nextUpdate = 0;
		long long window = 0;
		double growth = 1.0;
		bool v1Compat;
		// time spent in this node and in its temperature updates (PT::setProfiling)
		bool prof = false;
		long long runNs = 0;
		long long updateNs = 0;
		ReplicaEst replicaEstL;
		ReplicaEst replicaEstR;
		ptapi::Rng gen;
		std::uniform_real_distribution<double> dis{0.0, 1.0};
	public:
		// UPL_: proposals between two temperature updates (<= 0: never). growth_ > 1:
		// UPL_ is only the first interval and each of the following is growth_ times longer.
		NodeSwap(NodeMCMC<S>* left_, NodeMCMC<S>* right_, NodeTempUp<S>* tempUp_, Consumer<S>* pool_, int UPL_, ptapi::Rng gen_, bool v1Compat_ = false, int index_ = 0, bool evenOdd_ = false, double growth_ = 1.0)
		: left(left_), right(right_), tempUp(tempUp_), pool(pool_), index(index_), evenOdd(evenOdd_), nextUpdate(UPL_ > 0 ? UPL_ : -1), window(UPL_), growth(growth_), v1Compat(v1Compat_), gen(gen_) {}

		bool run() override {
			if(!prof) return runExchange();
			const std::chrono::steady_clock::time_point p0 = std::chrono::steady_clock::now();
			const bool again = runExchange();
			runNs += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - p0).count();
			return again;
		}

		void setProfiling(bool on){ prof = on; }
		long long getRunNs() const { return runNs; }        // includes getUpdateNs()
		long long getUpdateNs() const { return updateNs; }

		bool runExchange(){
			if(pool->stopped()) return false;
			++round;
			// Even-odd scheme: pairs (0,1), (2,3), ... propose on odd rounds and pairs
			// (1,2), (3,4), ... on even rounds. The node still runs, to keep the graph.
			if(tempUp) tempUp->propagate();
			if(evenOdd && (((index + round) & 1) == 0)) return true;
			++attempts;

			double temp1 = left->getTemp();
			double temp2 = right->getTemp();
			double energia1 = left->getEnergia();
			double energia2 = right->getEnergia();
			if(v1Compat){ energia1 = (double)(int)energia1; energia2 = (double)(int)energia2; }

			if(tempUp){
				replicaEstL.Push(energia1);
				replicaEstR.Push(energia2);
			}

			// Calculation of the exponent of exp
			double deltaBeta = ((1.0/temp1) - (1.0/temp2)) * (energia1-energia2);

			// When delta is greater than or equal to 0, it always changes
			const double prob = (deltaBeta >= 0) ? 1.0 : std::exp(deltaBeta);
			if(tempUp) tempUp->addProbability(prob);
			if(deltaBeta >= 0 || dis(gen) <= prob){
				swap();
				++accept;
				++acceptTotal;
			}

			if(tempUp && (attempts == nextUpdate)){
				std::chrono::steady_clock::time_point u0;
				if(prof) u0 = std::chrono::steady_clock::now();
				tempUp->setacceptRate((double)accept/(double)(v1Compat ? (long long)round : window));
				tempUp->setEnergy(replicaEstL.Mean(), replicaEstR.Mean());
				tempUp->run();
				accept = 0;
				replicaEstL.Clear();
				replicaEstR.Clear();
				if(growth > 1.0) window = std::max(window + 1, (long long)((double)window*growth + 0.5));
				nextUpdate += window;
				if(prof) updateNs += std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now() - u0).count();
			}

			return true;
		}

		void swap(){
			left->swapWith(*right);

			// try to set the labels up and down
			left->trySetLabels();
			right->trySetLabels();

			// update flow
			left->updateFlow();
			right->updateFlow();
		}

		double getAcceptRate() const { return attempts ? (double)acceptTotal/(double)attempts : 0.0; }
};

#endif
