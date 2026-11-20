#ifndef __NODETEMPUP_H__
#define __NODETEMPUP_H__

#include <cmath>
#include <algorithm>
#include <limits>
#include <vector>
#include <atomic>
#include <stdexcept>
#include "NodeMCMC.h"

// TU: dynamic update of the temperatures. It is executed inside the exchange
// node it belongs to (every 'tempUpdate' exchange proposals).
// 'left' is the replica with the HIGHER temperature of the pair, 'right' the lower.
//
// Since v2 the first and the last temperature of the ladder are ANCHORED: the
// updates only redistribute the temperatures between them (the v1 rules of
// types 1 and 2, which let the ladder drift, are kept behind
// PT::setV1Compatibility).
template<typename S>
class NodeTempUp{
	private:
		// Values of update type 3 that the pair of the first temperature writes for
		// the other pairs while their exchange nodes may be running.
		struct Relaxed{
			std::atomic<double> v{0.0};
			operator double() const { return v.load(std::memory_order_relaxed); }
			Relaxed& operator=(double x){ v.store(x, std::memory_order_relaxed); return *this; }
		};
		NodeMCMC<S>* left;
		NodeMCMC<S>* right;
		NodeTempUp<S>* downTempUp;
		NodeTempUp<S>* upTempUp = nullptr;
		int m = 0;
		int ut = 0;
		Relaxed c;
		Relaxed fl;
		Relaxed fr;
		double deltaE = 0;
		double energyL = 0;
		double energyR = 0;
		double acceptRate = 0;
		bool anchored = true;
		// update types 5 and 6: exchange probabilities of the pair since the last update
		Relaxed probSum;
		Relaxed probCount;
		// update type 5: log(log(T[i+1]/T[i])) of the pair, number of updates, upper limit
		double rho = 0.0;
		bool rhoReady = false;
		long long saUpdates = 0;
		double tCap = std::numeric_limits<double>::infinity();
		double saTarget = 0.234;
		// pilot runs of PT::calibrateLadder: mean of rho over the updates of this execution
		double rhoSum = 0.0;
		long long rhoCount = 0;

		static bool valid(double t){ return std::isfinite(t) && (t > 0.0); }
	public:
		NodeTempUp(NodeMCMC<S>* left_, NodeMCMC<S>* right_, NodeTempUp<S>* downTempUp_, int ut_, int tL, bool anchored_ = true, double tCap_ = std::numeric_limits<double>::infinity())
		: left(left_), right(right_), downTempUp(downTempUp_), m(tL), ut(ut_), anchored(anchored_), tCap(tCap_) {}

		void run(){
			switch (ut){
				case 1: if(anchored) tempUp23Anchored(); else tempUp23(); break;
				case 2: tempUpEnergy(); break;   // it never changes the first and the last temperature
				case 3: tempUpFB(); break;           // feedback-optimized as in the thesis
				case 4: tempUpFBAnchored(); break;   // feedback-optimized, quantiles of the density
				case 5: tempUpSA(); break;           // stochastic approximation towards 23.4%
				case 6: tempUpRejection(); break;    // equal rejection rates
				default: break;
			}
		}

		// Called by the exchange node on every proposal with the probability of the exchange.
		void addProbability(double p){
			if(ut >= 5){ probSum = (double)probSum + p; probCount = (double)probCount + 1.0; }
		}

		// Update type 5 (Miasojedow, Moulines and Vihola, 2013). The ladder is written
		// as T[i+1] = T[i]*exp(exp(rho[i])) and each pair moves its own rho with a
		// decreasing step:
		//     rho[i] += (n+1)^-0.6 * (average exchange probability - target)      target = 0.234
		// The first temperature is fixed and the others follow from the rhos, so the
		// last temperature is a result; it is limited to the highest temperature of
		// the initial ladder. Deterministic: each pair only writes its hotter temperature.
		void tempUpSA(){
			const double count = probCount;
			if(!rhoReady || !(count > 0.0)) return;
			const double h = ((double)probSum/count) - saTarget;
			rho += std::pow((double)(saUpdates + 1), -0.6)*h;
			rho = std::min(std::log(5.0), std::max(std::log(1.0e-4), rho));   // ratio between 1.0001 and 148
			++saUpdates;
			rhoSum += rho;
			++rhoCount;
			probSum = 0.0;
			probCount = 0.0;
		}
		// Type 5: the hotter temperature of the pair follows the colder one (every round).
		void propagate(){
			if(ut != 5) return;
			const double tR = right->getTemp();
			if(!rhoReady){
				const double tL = left->getTemp();
				if(!(tL > tR) || !(tR > 0.0)) return;       // needs an ascending ladder
				rho = std::log(std::log(tL/tR));
				rhoReady = true;
				return;
			}
			const double t = std::min(tCap, tR*std::exp(std::exp(rho)));
			if(valid(t)) left->setTemp(t);
		}

		// Update type 6 (Syed et al., 2019): equal rejection rates with fixed end
		// temperatures. With r[i] the average rejection probability of pair i, the
		// cumulative barrier L(k) = r[0] + ... + r[k-1] is interpolated linearly in
		// 1/T and the new temperatures are the points where it is k/(n-1) of the
		// total. Executed by the pair of the first temperature for the whole ladder;
		// half of the step is applied.
		void tempUpRejection(){
			if(downTempUp) return;
			std::vector<NodeMCMC<S>*> nodes;
			std::vector<NodeTempUp<S>*> pairs;
			nodes.push_back(right);
			for(NodeTempUp<S>* tu = this; tu; tu = tu->upTempUp){ nodes.push_back(tu->left); pairs.push_back(tu); }
			const int n = (int)nodes.size();
			if(n < 3) return;

			std::vector<double> beta(n), r(n-1);
			for(int k = 0; k < n; ++k){
				const double t = nodes[k]->getTemp();
				if(!valid(t)) return;
				beta[k] = 1.0/t;
			}
			double total = 0.0;
			const bool up = beta[1] < beta[0];
			for(int k = 0; k + 1 < n; ++k){
				if((beta[k+1] < beta[k]) != up || beta[k+1] == beta[k]) return;   // not monotonic
				const double count = pairs[k]->probCount;
				if(!(count > 0.0)) return;                                      // a pair without data
				r[k] = std::min(1.0, std::max(0.0, 1.0 - (double)pairs[k]->probSum/count));
				total += r[k];
			}
			if(!(total > 0.0)) return;                                           // nothing is rejected

			double cum = 0.0;
			int j = 0;
			for(int k = 1; k + 1 < n; ++k){
				const double goal = total * (double)k/(double)(n-1);
				while((j < n-2) && (cum + r[j] < goal)){ cum += r[j]; ++j; }
				const double frac = (r[j] > 0.0) ? std::min(1.0, (goal - cum)/r[j]) : 0.0;
				const double newBeta = beta[j] + frac*(beta[j+1] - beta[j]);
				const double t = 1.0/(0.5*(beta[k] + newBeta));
				if(valid(t)) nodes[k]->setTemp(t);
			}
			for(NodeTempUp<S>* tu : pairs){ tu->probSum = 0.0; tu->probCount = 0.0; }
		}

		// Target exchange rate of 23% with fixed end temperatures.
		// The pair changes its own gap (in 1/T) towards the gap that gives 23%:
		// narrower when the measured rate is below 20%, wider when it is above 26%.
		// The change is shared by the two temperatures of the pair, except the
		// first and the last of the ladder, which never move. Where 23% cannot be
		// reached by every pair (too few or too many replicas for the interval),
		// the pulls of adjacent pairs balance and the deviation is spread.
		// It only reads data that the dataflow graph orders (deterministic for an
		// update interval >= 2).
		void tempUp23Anchored(){
			const bool moveL = (upTempUp != nullptr), moveR = (downTempUp != nullptr);
			if(!moveL && !moveR) return;                          // two replicas: both anchored
			if(acceptRate >= 0.20 && acceptRate <= 0.26) return;
			const double dE = energyL - energyR;
			const double tL = left->getTemp(), tR = right->getTemp();
			if(!(dE > 0.0) || !std::isfinite(dE) || !(tL > tR) || !(tR > 0.0)) return;

			const double gap = (1.0/tR) - (1.0/tL);
			const double want = -std::log(0.23)/dE;
			double delta = want - gap;                            // > 0: widen
			if((acceptRate < 0.20 && delta >= 0.0) || (acceptRate > 0.26 && delta <= 0.0)) return;
			delta *= 0.5;                                         // damping (the estimates are noisy)
			const double share = (moveL && moveR) ? 0.5 : 1.0;

			if(moveR){
				double b = (1.0/tR) + share*delta;
				const double below = downTempUp->getDownTemp();   // next colder temperature
				if(below > 0.0 && below < tR) b = std::min(b, 0.5*((1.0/tR) + (1.0/below)));
				if(valid(1.0/b)) right->setTemp(1.0/b);
			}
			if(moveL){
				double b = (1.0/tL) - share*delta;
				const double above = upTempUp->getLeftTemp();     // next hotter temperature
				if(above > tL) b = std::max(b, 0.5*((1.0/tL) + (1.0/above)));
				if(valid(1.0/b)) left->setTemp(1.0/b);
			}
		}

		// v1: target exchange rate of 23% (the ladder is not anchored)
		void tempUp23(){
			if(acceptRate < 0.2){
				double k = (energyR - energyL)/right->getTemp();
				double newT = (energyL - energyR) / (std::log(0.23) - k);
				if(valid(newT)) left->setTemp(newT);
			}
		}

		// Equal exchange probability between adjacent temperatures (Hukushima)
		void tempUpEnergy(){
			deltaE = energyL - energyR;

			if(downTempUp){
				double k1 = downTempUp->getDeltaE();
				double tMenos = downTempUp->getDownTemp();
				double eq = ((k1/tMenos) + (deltaE/left->getTemp()));
				double newT = ((deltaE+k1)/eq);
				// v2: the ladder must stay ordered (v1 accepted any positive value)
				if(valid(newT) && (!anchored || ((newT > tMenos) && (newT < left->getTemp())))) right->setTemp(newT);
			}
		}

		// Update type 4: feedback-optimized (Katzgraber et al., 2006), computed for
		// the whole ladder at once, with fixed end temperatures.
		// f(T) is the fraction of the replicas that visited temperature T coming from
		// the first temperature. The density of temperatures is made proportional to
		// sqrt((df/dT)/deltaT): the new temperatures are the points that divide the
		// integral of this density, between the first and the last temperature, in
		// equal parts. Executed by the pair of the first temperature for the whole
		// ladder, when every temperature has flow data; half of the step is applied.
		void tempUpFBAnchored(){
			if(downTempUp) return;
			std::vector<NodeMCMC<S>*> nodes;
			nodes.push_back(right);
			for(NodeTempUp<S>* tu = this; tu; tu = tu->upTempUp) nodes.push_back(tu->left);
			const int n = (int)nodes.size();
			if(n < 3) return;

			std::vector<double> T(n), f(n), w(n-1);
			for(int k = 0; k < n; ++k){
				if(!nodes[k]->hasFlow()) return;
				T[k] = nodes[k]->getTemp();
				f[k] = nodes[k]->getFlow();
			}
			double total = 0.0;
			const bool up = T[1] > T[0];
			for(int k = 0; k + 1 < n; ++k){
				if((T[k+1] > T[k]) != up || T[k+1] == T[k]) return;   // not monotonic
				w[k] = std::sqrt(std::fabs(f[k] - f[k+1]));         // integral of the density in the interval
				total += w[k];
			}
			if(!(total > 0.0) || !std::isfinite(total)) return;

			double cum = 0.0;
			int j = 0;
			for(int k = 1; k + 1 < n; ++k){
				const double goal = total * (double)k/(double)(n-1);
				while((j < n-2) && (cum + w[j] < goal)){ cum += w[j]; ++j; }
				const double frac = (w[j] > 0.0) ? std::min(1.0, (goal - cum)/w[j]) : 0.0;
				const double newT = T[j] + frac*(T[j+1] - T[j]);
				const double t = 0.5*(T[k] + newT);
				if(valid(t)) nodes[k]->setTemp(t);
			}
			for(int k = 0; k < n; ++k) nodes[k]->resetFlow();
		}

		// Feedback-optimized as described in the thesis (update type 3)
		void tempUpFB(){
			double upper = upTempUp ? upTempUp->getLeftTemp() : std::numeric_limits<double>::infinity();

			if(!downTempUp){
				double density = this->getProbDens();
				double dfdt = ((fl - fr)/(left->getTemp() - right->getTemp()));
				if((dfdt != 0) && (density != 0) && std::isfinite(density)){
					c = 1.0/density;
					double eq = std::pow((1.0/(double)(m*c)),2)/std::fabs(dfdt);
					double newT = eq + right->getTemp();
					if((newT > right->getTemp()) && (newT < upper)) left->setTemp(newT);
				}
			}else if(upTempUp){
				c = downTempUp->getC();
				double dfdt = ((fl - fr)/(left->getTemp() - right->getTemp()));
				if((dfdt != 0) && (c != 0)){
					double eq = std::fabs(std::pow((1.0/(double)(m*c)),2)/dfdt);
					double newT = eq + right->getTemp();
					if((newT > right->getTemp()) && (newT < upper)) left->setTemp(newT);
				}
			}
		}

		double getProbDens(){
			fr = right->getFlow();
			fl = left->getFlow();
			if((fr == 0) || (fl == 0)) return 0;

			right->resetFlow();
			double den = left->getTemp() - right->getTemp();
			double dfdt = std::fabs((fl - fr)/den);
			double atual = std::sqrt(dfdt/den)*den;
			if(upTempUp) atual += upTempUp->getProbDens();
			return atual;
		}

		// Used by PT::calibrateLadder: the target rate of type 5 and the number of
		// updates already made (the step keeps decreasing from one pilot run to the next).
		void setSA(double target, long long updatesDone){ saTarget = target; saUpdates = updatesDone; }
		long long getSAUpdates() const { return saUpdates; }
		// mean of rho over the updates of this execution (the last value if there was none)
		double getMeanRho() const { return rhoCount ? rhoSum/(double)rhoCount : rho; }
		bool hasRho() const { return rhoReady; }

		void setUpTempUp(NodeTempUp<S>* utu){ upTempUp = utu; }
		void setacceptRate(double aR){ acceptRate = aR; }
		void setEnergy(double l, double r){ energyL = l; energyR = r; }
		double getDeltaE() const { return deltaE; }
		double getDownTemp() const { return right->getTemp(); }
		double getLeftTemp() const { return left->getTemp(); }
		double getC() const { return (double)c; }
};

#endif
