#ifndef __PROBLEM_H__
#define __PROBLEM_H__

#include <vector>
#include <climits>
#include <random>
#include <cstdint>
#include <utility>
#include <stdexcept>
#include <type_traits>

// ---------------------------------------------------------------------------
// PTAPI v2 - problem interface
// ---------------------------------------------------------------------------

struct solution{
  double evalSol = INT_MAX;
  // Kept only for source compatibility with v1. Since v2 the round-trip labels
  // are stored by the API itself; the user no longer needs to copy them.
  bool Nup = false, Ndown = false;
};

namespace ptapi{
	typedef std::mt19937_64 Rng;

	namespace detail{
		inline Rng*& currentRng(){
			thread_local Rng* r = nullptr;
			return r;
		}
		// Index of the move returned by the last Problem::chooseMove() on this thread.
		inline int& lastMove(){
			thread_local int m = -1;
			return m;
		}
	}

	// Random engine of the replica that is currently running on this thread.
	// Using it inside construction()/neighbor()/evaluate() instead of creating a
	// std::random_device + engine on every call is much faster and makes the whole
	// run reproducible through PT::setSeed().
	inline Rng& rng(){
		Rng* r = detail::currentRng();
		if(r) return *r;
		thread_local Rng fallback{std::random_device{}()};
		return fallback;
	}

	// Uniform integer in [lo, hi]
	inline int randInt(int lo, int hi){
		return std::uniform_int_distribution<int>(lo, hi)(rng());
	}

	// Uniform real in [0, 1)
	inline double rand01(){
		return std::uniform_real_distribution<double>(0.0, 1.0)(rng());
	}
}

template<typename S>
class Problem{
	private:
		// Coupling with the user's evaluate(), bound by PT::start() (see ptapiBind).
		void (*ptapiStep_)(Problem<S>*, const S&, S&) = nullptr;
		void (*ptapiEval_)(Problem<S>*, S&) = nullptr;
		std::vector<double> moveCum_;   // cumulative probabilities of the registered moves

		template<typename P, typename = void>
		struct HasEvaluate: std::false_type{};
		template<typename P>
		struct HasEvaluate<P, decltype(void(std::declval<P&>().evaluate(std::declval<S&>())))>: std::true_type{};

		// Evaluates 'n' with the user's function, whatever its form, and stores it in 'out'.
		template<typename P>
		static void evalInto(P* p, S& n, S& out){
			if constexpr (std::is_arithmetic<typename std::decay<decltype(p->evaluate(n))>::type>::value){
				double v = p->evaluate(n);           // double evaluate(S) or double evaluate(const S&)
				n.evalSol = v;
				out = std::move(n);
			}else{
				out = p->evaluate(std::move(n));     // S evaluate(S)
			}
		}

		template<typename P>
		static void classicStep(Problem<S>* self, const S& current, S& candidate){
			P* p = static_cast<P*>(self);
			S n = p->neighbor(current);
			evalInto<P>(p, n, candidate);
		}

		template<typename P>
		static void classicEval(Problem<S>* self, S& s){
			S n = std::move(s);
			evalInto<P>(static_cast<P*>(self), n, s);
		}
	public:
		Problem(){}
		virtual ~Problem(){}
		void setBestSol(S){} // kept for v1 compatibility (no effect)

		// --- required ---
		virtual S construction()=0;   // must return an evaluated solution
		virtual S neighbor(S sol)=0;

		// --- optional (v2) ---
		// Initial solution of one replica: 'index' goes from 0 (first temperature) to
		// 'count' - 1 (last temperature). The default calls construction(). Override it
		// to give each replica its own constructive heuristic (e.g. a greedy solution
		// for the cold replicas and random ones for the hot replicas). It must return
		// an evaluated solution. With PT::setParallelConstruction(true) it is called
		// from several threads at the same time.
		virtual S constructionReplica(int index, int count){ (void)index; (void)count; return construction(); }

		// The evaluation function is declared by the user class in one of the forms
		//     S      evaluate(S sol);          (may also change the solution)
		//     double evaluate(S sol);          (v1)
		//     double evaluate(const S& sol);   (no copy of the solution)
		// It is not needed if neighborEval() is overridden and no initial
		// solutions are given to PT.

		// --- optional (v2) ---
		// One Markov chain step: writes into 'candidate' an EVALUATED neighbor of
		// 'current'. The default calls neighbor() + evaluate(). Override it to avoid
		// copies (reuse the memory already held by 'candidate') and/or to compute
		// candidate.evalSol incrementally from current.evalSol (delta evaluation).
		virtual void neighborEval(const S& current, S& candidate){
			if(!ptapiStep_) throw std::logic_error("PTAPI: the problem class has no evaluate(); pass a pointer of the problem class itself to PT::start() or override neighborEval()");
			ptapiStep_(this, current, candidate);
		}

		// Several kinds of move (v2). Register the probability of each one, usually
		// in the constructor of the problem class:
		//     setMoveProbabilities({0.6, 0.3, 0.1});
		// and select the move inside neighbor() / neighborEval():
		//     switch(chooseMove()){ case 0: ...; case 1: ...; case 2: ...; }
		// The values are weights (they are normalized). PT::getMoveStatistics()
		// reports how often each move was proposed and accepted.
		void setMoveProbabilities(const std::vector<double>& weights){
			double sum = 0.0;
			for(double w : weights){
				if(!(w >= 0.0)) throw std::invalid_argument("PTAPI: move probabilities must be >= 0");
				sum += w;
			}
			if(weights.empty() || !(sum > 0.0)) throw std::invalid_argument("PTAPI: at least one move must have a positive probability");
			moveCum_.clear();
			double acc = 0.0;
			for(double w : weights){ acc += w/sum; moveCum_.push_back(acc); }
			moveCum_.back() = 1.0;
		}
		int getMoveCount() const { return (int)moveCum_.size(); }
		std::vector<double> getMoveProbabilities() const {
			std::vector<double> p;
			for(size_t i = 0; i < moveCum_.size(); ++i) p.push_back(moveCum_[i] - (i ? moveCum_[i-1] : 0.0));
			return p;
		}

		// Draws the index of a move according to the registered probabilities.
		int chooseMove(){
			if(moveCum_.empty()) throw std::logic_error("PTAPI: chooseMove() called before setMoveProbabilities()");
			const double u = ptapi::rand01();
			int m = 0;
			const int last = (int)moveCum_.size() - 1;
			while(m < last && u >= moveCum_[m]) ++m;
			ptapi::detail::lastMove() = m;
			return m;
		}

		// --- used by PT ---
		template<typename P>
		void ptapiBind(){
			if constexpr (HasEvaluate<P>::value){
				ptapiStep_ = &classicStep<P>;
				ptapiEval_ = &classicEval<P>;
			}
		}

		void ptapiEvaluate(S& s){
			if(!ptapiEval_) throw std::logic_error("PTAPI: the problem class has no evaluate()");
			ptapiEval_(this, s);
		}
};

#endif
