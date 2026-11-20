/*
 * PTAPI - unit tests.
 *
 * No external framework: only the C++ standard library, like the API itself.
 *
 *   make test            (in this folder)
 *   ./test_ptapi         runs everything
 *   ./test_ptapi name    runs the tests whose name contains "name"
 *
 * The tests only use small synthetic problems with three different encodings
 * (permutation, binary and continuous): the API must work for any optimization
 * problem by changing only the problem class.
 *
 * To add a test: write TEST(name){ ... CHECK(condition); ... }.
 */
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <thread>
#include <cmath>
#include <atomic>
#include <chrono>
#include <string>
#include <vector>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <sys/resource.h>
#include "../include/PT.h"
#include "../include/ExecTime.h"
#include "../include/Experiment.h"
#include "../include/Comparison.h"
#include "../include/Options.h"
#include <set>
#include <mutex>
#include <filesystem>

// ---------------------------------------------------------------------------
// minimal test framework
// ---------------------------------------------------------------------------
struct TestCase{ const char* name; void (*fn)(); };
static std::vector<TestCase>& registry(){ static std::vector<TestCase> r; return r; }
static int gFailures = 0;
struct Registrar{ Registrar(const char* n, void (*f)()){ registry().push_back({n, f}); } };

#define TEST(name) static void test_##name(); static Registrar reg_##name(#name, test_##name); static void test_##name()
#define CHECK(cond) do{ if(!(cond)){ ++gFailures; std::printf("    FAILED %s:%d  %s\n", __FILE__, __LINE__, #cond); } }while(0)
#define CHECK_NEAR(a, b, tol) do{ double va = (a), vb = (b); if(!(std::fabs(va - vb) <= (tol))){ ++gFailures; std::printf("    FAILED %s:%d  %s = %g, expected %g (tol %g)\n", __FILE__, __LINE__, #a, va, vb, (double)(tol)); } }while(0)
#define CHECK_THROWS(expr, type) do{ bool thrown = false; try{ expr; }catch(const type&){ thrown = true; }catch(...){} if(!thrown){ ++gFailures; std::printf("    FAILED %s:%d  %s did not throw %s\n", __FILE__, __LINE__, #expr, #type); } }while(0)

static double cpuSeconds(){
	struct rusage r; getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_utime.tv_usec*1e-6 + r.ru_stime.tv_sec + r.ru_stime.tv_usec*1e-6;
}

// ---------------------------------------------------------------------------
// Problem 1 (permutation): put a permutation in order. f = sum |p[i] - i|, optimum 0.
// Written in the four ways the API accepts; with the same seed they must give
// exactly the same run.
// ---------------------------------------------------------------------------
struct solPerm: public solution{ std::vector<int> sol; };

static double permCost(const std::vector<int>& p){
	double c = 0;
	for(int i = 0; i < (int)p.size(); i++) c += std::abs(p[i] - i);
	return c;
}

static bool isPermutation(std::vector<int> p){
	std::sort(p.begin(), p.end());
	for(int i = 0; i < (int)p.size(); i++) if(p[i] != i) return false;
	return true;
}

// common part: construction and move
template<typename Derived>
class PermBase: public Problem<solPerm>{
	public:
		int n;
		std::atomic<double> minSeen{1e300};
		std::atomic<long long> calls{0};
		explicit PermBase(int n_): n(n_){}

		void note(double v){
			calls.fetch_add(1, std::memory_order_relaxed);
			double c = minSeen.load();
			while(v < c && !minSeen.compare_exchange_weak(c, v)){}
		}

		solPerm construction(){
			solPerm s;
			s.sol.resize(n);
			std::iota(s.sol.begin(), s.sol.end(), 0);
			std::shuffle(s.sol.begin(), s.sol.end(), ptapi::rng());
			s.evalSol = permCost(s.sol);
			note(s.evalSol);
			return s;
		}

		solPerm neighbor(solPerm s){
			int a = ptapi::randInt(0, n-1), b = ptapi::randInt(0, n-1);
			std::swap(s.sol[a], s.sol[b]);
			return s;
		}
};

// form used by the "source" folder: S evaluate(S)
class PermS: public PermBase<PermS>{
	public:
		using PermBase::PermBase;
		solPerm evaluate(solPerm s){ s.evalSol = permCost(s.sol); note(s.evalSol); return s; }
};

// form of the published v1: double evaluate(S)
class PermD: public PermBase<PermD>{
	public:
		using PermBase::PermBase;
		double evaluate(solPerm s){ double v = permCost(s.sol); note(v); return v; }
};

// v2: double evaluate(const S&)
class PermC: public PermBase<PermC>{
	public:
		using PermBase::PermBase;
		double evaluate(const solPerm& s){ double v = permCost(s.sol); note(v); return v; }
};

// v2: neighborEval with delta evaluation (and no evaluate() at all)
class PermFast: public PermBase<PermFast>{
	public:
		using PermBase::PermBase;
		void neighborEval(const solPerm& cur, solPerm& cand){
			int a = ptapi::randInt(0, n-1), b = ptapi::randInt(0, n-1);
			cand.sol = cur.sol;
			std::swap(cand.sol[a], cand.sol[b]);
			cand.evalSol = cur.evalSol - std::abs(cur.sol[a] - a) - std::abs(cur.sol[b] - b)
			                           + std::abs(cur.sol[b] - a) + std::abs(cur.sol[a] - b);
			note(cand.evalSol);
		}
};

// neighborEval overridden in a class derived from one that has evaluate()
class PermSFast: public PermS{
	public:
		using PermS::PermS;
		std::atomic<long long> fastCalls{0};
		void neighborEval(const solPerm& cur, solPerm& cand){
			fastCalls.fetch_add(1, std::memory_order_relaxed);
			cand = evaluate(neighbor(cur));
		}
};

// ---------------------------------------------------------------------------
// Problem 2 (binary): OneMax as a minimization. f = number of zeros, optimum 0.
// ---------------------------------------------------------------------------
struct solBin: public solution{ std::vector<char> bits; };

class OneMax: public Problem<solBin>{
	public:
		int n;
		explicit OneMax(int n_): n(n_){}
		solBin construction(){
			solBin s;
			s.bits.resize(n);
			for(char& b : s.bits) b = (char)ptapi::randInt(0, 1);
			s.evalSol = evaluate(s);
			return s;
		}
		solBin neighbor(solBin s){
			int i = ptapi::randInt(0, n-1);
			s.bits[i] = !s.bits[i];
			return s;
		}
		double evaluate(const solBin& s){
			double z = 0;
			for(char b : s.bits) z += (b == 0);
			return z;
		}
};

// ---------------------------------------------------------------------------
// Problem 3 (continuous): sphere function in [-0.5, 0.5]^d. f = sum x^2 (always < 1).
// ---------------------------------------------------------------------------
struct solReal: public solution{ std::vector<double> x; };

class Sphere: public Problem<solReal>{
	public:
		int d;
		explicit Sphere(int d_): d(d_){}
		solReal construction(){
			solReal s;
			s.x.resize(d);
			for(double& v : s.x) v = ptapi::rand01() - 0.5;
			s.evalSol = evaluate(s);
			return s;
		}
		solReal neighbor(solReal s){
			int i = ptapi::randInt(0, d-1);
			s.x[i] += 0.1*(ptapi::rand01() - 0.5);
			if(s.x[i] > 0.5) s.x[i] = 0.5;
			if(s.x[i] < -0.5) s.x[i] = -0.5;
			return s;
		}
		double evaluate(solReal s){
			double f = 0;
			for(double v : s.x) f += v*v;
			return f;
		}
};

// ---------------------------------------------------------------------------
// Problem 4: six states on a ring, used to check the sampled distribution.
// ---------------------------------------------------------------------------
struct solState: public solution{ int state = 0; };
static const double kRingEnergy[6] = {0.0, 1.0, 2.0, 3.0, 1.0, 0.5};

class Ring: public Problem<solState>{
	public:
		solState construction(){
			solState s;
			s.state = ptapi::randInt(0, 5);
			s.evalSol = kRingEnergy[s.state];
			return s;
		}
		solState neighbor(solState s){
			s.state = (s.state + (ptapi::randInt(0, 1) ? 1 : 5)) % 6; // symmetric proposal
			return s;
		}
		double evaluate(const solState& s){ return kRingEnergy[s.state]; }
};

// Problem that fails after some evaluations
class Faulty: public PermBase<Faulty>{
	public:
		using PermBase::PermBase;
		double evaluate(const solPerm& s){
			if(calls.fetch_add(1) > 500) throw std::runtime_error("evaluation failed");
			return permCost(s.sol);
		}
};

// Problem with an evaluation that costs some time
class Slow: public PermBase<Slow>{
	public:
		using PermBase::PermBase;
		double evaluate(const solPerm& s){
			volatile double w = 0;
			for(int i = 0; i < 20000; i++) w = w + std::sqrt((double)i);
			return permCost(s.sol);
		}
};

// Two kinds of move chosen with registered probabilities (75% / 25%)
class PermMoves: public PermBase<PermMoves>{
	public:
		explicit PermMoves(int n_): PermBase(n_){ setMoveProbabilities({3.0, 1.0}); }
		solPerm neighbor(solPerm s){
			int a = ptapi::randInt(0, n-1), b = ptapi::randInt(0, n-1);
			switch(chooseMove()){
				case 0: std::swap(s.sol[a], s.sol[b]); break;
				default: { int v = s.sol[a]; s.sol.erase(s.sol.begin()+a); s.sol.insert(s.sol.begin()+b, v); }
			}
			return s;
		}
		double evaluate(const solPerm& s){ return permCost(s.sol); }
};

// functions that write and read one solution (checkpoint)
static void writePerm(std::ostream& out, const solPerm& s){ out << s.sol.size(); for(int v : s.sol) out << ' ' << v; }
static void readPerm(std::istream& in, solPerm& s){ size_t k = 0; in >> k; if(k > 100000) k = 0; s.sol.resize(k); for(int& v : s.sol) in >> v; }

// ---------------------------------------------------------------------------
// tests
// ---------------------------------------------------------------------------

TEST(temperature_distributions){
	for(int dist = 1; dist <= 4; dist++){
		PT<solPerm> algo(0.5f, 8.0f, 9, 10, 10, dist, 0, 50);
		std::vector<double> t = algo.getInitialTemperatures();
		CHECK(t.size() == 9);
		// Distribution 2 (inverse linear) goes from the maximum to the minimum
		// temperature, as in v1 ("the index is inverted"); the others are ascending.
		if(dist == 2) std::reverse(t.begin(), t.end());
		CHECK_NEAR(t.front(), 0.5, 1e-5);
		CHECK_NEAR(t.back(), 8.0, 1e-4);
		for(size_t i = 1; i < t.size(); i++) CHECK(t[i] > t[i-1]);
	}
	// linear: constant difference; geometric (3 and 4): constant ratio
	std::vector<double> lin = PT<solPerm>(1.0f, 9.0f, 5, 1, 1, 1, 0, 1).getInitialTemperatures();
	for(size_t i = 1; i < lin.size(); i++) CHECK_NEAR(lin[i] - lin[i-1], 2.0, 1e-6);
	std::vector<double> inv = PT<solPerm>(1.0f, 9.0f, 5, 1, 1, 2, 0, 1).getInitialTemperatures();
	CHECK(inv.front() > inv.back());
	for(size_t i = 2; i < inv.size(); i++) CHECK_NEAR(1.0/inv[i] - 1.0/inv[i-1], 1.0/inv[i-1] - 1.0/inv[i-2], 1e-9);   // constant step in 1/T
	for(int dist = 3; dist <= 4; dist++){
		std::vector<double> geo = PT<solPerm>(1.0f, 16.0f, 5, 1, 1, dist, 0, 1).getInitialTemperatures();
		for(size_t i = 1; i < geo.size(); i++) CHECK_NEAR(geo[i]/geo[i-1], 2.0, 1e-6);
	}
	// one replica: the minimum temperature
	CHECK_NEAR(PT<solPerm>(0.7f, 9.0f, 1, 1, 1, 1, 0, 1).getInitialTemperatures()[0], 0.7, 1e-6);
}

TEST(invalid_parameters_throw){
	CHECK_THROWS(PT<solPerm>(1, 2, 4, 10, 10, 0, 0, 50), std::invalid_argument);
	CHECK_THROWS(PT<solPerm>(1, 2, 4, 10, 10, 5, 0, 50), std::invalid_argument);
	CHECK_THROWS(PT<solPerm>(1, 2, 4, 10, 10, 1, 7, 50), std::invalid_argument);
	CHECK_THROWS(PT<solPerm>(1, 2, 4, 10, 10, 1, -1, 50), std::invalid_argument);
	CHECK_THROWS(PT<solPerm>(1, 2, 0, 10, 10, 1, 0, 50), std::invalid_argument);
	PT<solPerm> algo(1, 2, 4, 10, 10, 1, 0, 50);
	CHECK_THROWS(algo.setTemperatures(std::vector<double>()), std::invalid_argument);
}

TEST(basic_run_and_statistics){
	const int rep = 8, MCL = 50, PTL = 200;
	PermS prob(30);
	PT<solPerm> algo(0.5f, 6.0f, rep, MCL, PTL, 4, 0, 50);
	algo.setSeed(1);
	solPerm best = algo.start(4, &prob);

	CHECK(isPermutation(best.sol));
	CHECK_NEAR(best.evalSol, permCost(best.sol), 1e-9);       // value consistent with the solution
	CHECK_NEAR(best.evalSol, prob.minSeen.load(), 1e-9);      // the best evaluated solution is returned
	CHECK_NEAR(algo.getBestSol().evalSol, best.evalSol, 0);
	CHECK(algo.getEvaluations() == (long long)rep*MCL*PTL);
	CHECK(prob.calls.load() == (long long)rep*MCL*PTL + rep);  // + initial solutions
	CHECK(algo.getPTLEnd() == PTL);
	CHECK(algo.getBestPTLSolIndex() >= 0 && algo.getBestPTLSolIndex() < PTL);
	CHECK(algo.getThreadsUsed() == 4);
	CHECK(algo.getSeed() == 1);
	CHECK(algo.getElapsedMs() > 0 && algo.getTimeToBestMs() <= algo.getElapsedMs());
	CHECK(algo.getReplicas().size() == (size_t)rep);
	for(const solPerm& r : algo.getReplicas()){ CHECK(isPermutation(r.sol)); CHECK_NEAR(r.evalSol, permCost(r.sol), 1e-9); CHECK(r.evalSol >= best.evalSol); }
	CHECK(algo.getTemperatures().size() == (size_t)rep);
	CHECK(algo.getMoveAcceptanceRates().size() == (size_t)rep);
	CHECK(algo.getSwapAcceptanceRates().size() == (size_t)rep - 1);
	for(double r : algo.getMoveAcceptanceRates()) CHECK(r >= 0.0 && r <= 1.0);
	for(double r : algo.getSwapAcceptanceRates()) CHECK(r >= 0.0 && r <= 1.0);
	// hotter replicas accept more moves
	CHECK(algo.getMoveAcceptanceRates().back() > algo.getMoveAcceptanceRates().front());
}

TEST(finds_optimum_permutation){
	PermS prob(12);
	PT<solPerm> algo(0.3f, 4.0f, 8, 100, 300, 4, 0, 50);
	algo.setSeed(3);
	CHECK_NEAR(algo.start(4, &prob).evalSol, 0.0, 0);
}

TEST(finds_optimum_binary){
	OneMax prob(60);
	PT<solBin> algo(0.2f, 2.0f, 8, 100, 200, 4, 0, 50);
	algo.setSeed(3);
	solBin best = algo.start(4, &prob);
	CHECK_NEAR(best.evalSol, 0.0, 0);
	CHECK(std::count(best.bits.begin(), best.bits.end(), 1) == 60);
}

TEST(finds_optimum_continuous){
	Sphere prob(5);
	PT<solReal> algo(0.0005f, 0.1f, 8, 100, 300, 4, 0, 50);
	algo.setSeed(3);
	solReal best = algo.start(4, &prob);
	CHECK(best.evalSol < 0.002);
	CHECK_NEAR(best.evalSol, prob.evaluate(best), 1e-12);
}

TEST(deterministic_with_seed_for_any_thread_count){
	for(int upType : {0, 1, 2, 5}){
		std::vector<int> refSol; double refVal = -1; int refPTL = -1; std::vector<double> refTemps; std::vector<solPerm> refRep;
		for(int threads : {1, 3, 8}){
			PermS prob(30);
			PT<solPerm> algo(0.5f, 6.0f, 8, 40, 300, 4, upType, 50);
			algo.setSeed(77);
			solPerm best = algo.start(threads, &prob);
			if(refVal < 0){ refSol = best.sol; refVal = best.evalSol; refPTL = algo.getBestPTLSolIndex(); refTemps = algo.getTemperatures(); refRep = algo.getReplicas(); }
			CHECK(best.evalSol == refVal);
			CHECK(algo.getBestPTLSolIndex() == refPTL);
			CHECK(algo.getTemperatures() == refTemps);
			for(size_t i = 0; i < refRep.size(); i++) CHECK(algo.getReplicas()[i].sol == refRep[i].sol);
		}
	}
}

TEST(seed_changes_the_run_and_start_can_be_repeated){
	PermS prob(30);
	PT<solPerm> algo(0.5f, 6.0f, 6, 20, 50, 4, 0, 50);
	algo.setSeed(5);
	solPerm a = algo.start(3, &prob);
	std::vector<solPerm> ra = algo.getReplicas();
	solPerm b = algo.start(3, &prob);                 // same object, second run
	CHECK(a.sol == b.sol && a.evalSol == b.evalSol);
	algo.setSeed(6);
	algo.start(3, &prob);
	bool different = false;
	for(size_t i = 0; i < ra.size(); i++) different = different || (ra[i].sol != algo.getReplicas()[i].sol);
	CHECK(different);

	// without setSeed every run uses a new seed, which can be read back
	PT<solPerm> free1(0.5f, 6.0f, 4, 10, 10, 4, 0, 50);
	free1.start(2, &prob); unsigned long long s1 = free1.getSeed();
	free1.start(2, &prob); unsigned long long s2 = free1.getSeed();
	CHECK(s1 != s2);
}

template<typename P>
static void runPerm(P& prob, double& value, int& ptl, std::vector<int>& sol, int upType = 0){
	PT<solPerm> algo(0.5f, 6.0f, 8, 40, 200, 4, upType, 50);
	algo.setSeed(2024);
	solPerm best = algo.start(4, &prob);
	value = best.evalSol; ptl = algo.getBestPTLSolIndex(); sol = best.sol;
}

TEST(all_evaluate_forms_give_the_same_run){
	double v0, v; int p0, p; std::vector<int> s0, s;
	PermS ps(30);    runPerm(ps, v0, p0, s0);
	PermD pd(30);    runPerm(pd, v, p, s);  CHECK(v == v0 && p == p0 && s == s0);  // double evaluate(S)
	PermC pc(30);    runPerm(pc, v, p, s);  CHECK(v == v0 && p == p0 && s == s0);  // double evaluate(const S&)
	PermFast pf(30); runPerm(pf, v, p, s);  CHECK(v == v0 && p == p0 && s == s0);  // neighborEval + delta evaluation
	CHECK_NEAR(pf.minSeen.load(), ps.minSeen.load(), 1e-9);
}

TEST(neighborEval_is_used_through_a_base_class_pointer){
	PermSFast fast(20);
	PermS* base = &fast;
	PT<solPerm> algo(0.5f, 6.0f, 4, 25, 40, 4, 0, 50);
	algo.setSeed(1);
	algo.start(2, base);
	CHECK(fast.fastCalls.load() == algo.getEvaluations());
	CHECK(algo.getEvaluations() == 4LL*25*40);
}

TEST(problem_pointer_without_evaluate_reports_an_error){
	// start() receives a Problem<S>*: the API cannot see evaluate() and the class
	// does not override neighborEval(). The error must be reported, not crash.
	Problem<solPerm>* generic = new PermS(10);
	PT<solPerm> algo(0.5f, 6.0f, 4, 10, 10, 4, 0, 50);
	CHECK_THROWS(algo.start(2, generic), std::logic_error);
	delete generic;
}

TEST(exception_in_problem_code_reaches_the_caller){
	Faulty prob(20);
	PT<solPerm> algo(0.5f, 6.0f, 6, 50, 1000, 4, 0, 50);
	CHECK_THROWS(algo.start(4, &prob), std::runtime_error);
	// the object is still usable
	PermS ok(20);
	algo.setSeed(1);
	CHECK(isPermutation(algo.start(4, &ok).sol));
}

TEST(stop_by_time_limit){
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 50, 2000000000, 4, 0, 50);
	algo.setTimeLimit(0.25);
	ExecTime et;
	solPerm best = algo.start(4, &prob);
	int ms = et.getTimeMs();
	CHECK(ms >= 240 && ms < 1500);
	CHECK(algo.getPTLEnd() > 0 && algo.getPTLEnd() < 2000000000);
	CHECK_NEAR(best.evalSol, prob.minSeen.load(), 1e-9);   // nothing is lost when stopping by time
	CHECK(isPermutation(best.sol));

	PT<solPerm> v1style(0.5f, 6.0f, 4, 50, 2000000000, 4, 0, 50);
	v1style.setTimeEnd(1);                                 // v1 function, whole seconds
	ExecTime et2;
	v1style.start(2, &prob);
	ms = et2.getTimeMs();
	CHECK(ms >= 990 && ms < 2500);
}

TEST(stop_by_target_value){
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 50, 1000000, 4, 0, 50);
	algo.setSeed(9);
	algo.setTargetValue(100.0);
	solPerm best = algo.start(4, &prob);
	CHECK(best.evalSol <= 100.0);
	CHECK(algo.getPTLEnd() < 1000000);
}

TEST(stop_by_stagnation){
	PermS prob(12);
	const int PTL = 200000;
	PT<solPerm> algo(0.3f, 4.0f, 6, 30, PTL, 4, 0, 50);
	algo.setSeed(4);
	algo.setPTLStopEnd(0.01f);                       // 2000 rounds without a new best
	algo.start(4, &prob);
	CHECK(algo.getPTLEnd() < PTL);
	CHECK(algo.getPTLEnd() > algo.getBestPTLSolIndex() + 2000);
	CHECK(algo.getPTLEnd() <= algo.getBestPTLSolIndex() + 2000 + 50);
}

TEST(stop_requested_by_callback_and_by_another_thread){
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 50, 2000000000, 4, 0, 50);
	int improvements = 0;
	algo.setOnImprove([&](const solPerm& s, int ptl, double ms){
		CHECK(isPermutation(s.sol)); CHECK(ptl >= 0); CHECK(ms >= 0);
		if(++improvements >= 5) algo.requestStop();
	});
	algo.start(4, &prob);
	CHECK(improvements >= 5);
	CHECK((size_t)improvements == algo.getConvergence().size());

	PT<solPerm> other(0.5f, 6.0f, 8, 50, 2000000000, 4, 0, 50);
	std::thread stopper([&]{ std::this_thread::sleep_for(std::chrono::milliseconds(150)); other.requestStop(); });
	ExecTime et;
	other.start(4, &prob);
	stopper.join();
	CHECK(et.getTimeMs() >= 140 && et.getTimeMs() < 1500);
}

TEST(stop_requested_during_construction){
	// The first improvement is reported while the initial solutions are being
	// built; a stop requested there must not be lost.
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 50, 2000000000, 4, 0, 50);
	algo.setOnImprove([&](const solPerm&, int, double){ algo.requestStop(); });
	solPerm best = algo.start(4, &prob);
	CHECK(isPermutation(best.sol));
	CHECK(algo.getPTLEnd() == 0);
}

TEST(convergence_history){
	PermS prob(30);
	PT<solPerm> algo(0.5f, 6.0f, 8, 40, 200, 4, 0, 50);
	algo.setSeed(8);
	solPerm best = algo.start(4, &prob);
	const std::vector<PTImprovement>& h = algo.getConvergence();
	CHECK(!h.empty());
	for(size_t i = 1; i < h.size(); i++){ CHECK(h[i].value < h[i-1].value); CHECK(h[i].timeMs >= h[i-1].timeMs); }
	CHECK(h.back().value == best.evalSol);
	CHECK(h.back().ptl == algo.getBestPTLSolIndex());
	CHECK_NEAR(h.back().timeMs, algo.getTimeToBestMs(), 1e-9);
}

TEST(edge_cases_replicas_threads_and_lengths){
	PermS prob(20);
	{	// one replica: plain Metropolis, no exchanges
		PT<solPerm> algo(0.5f, 6.0f, 1, 30, 100, 4, 0, 50);
		solPerm best = algo.start(8, &prob);
		CHECK(isPermutation(best.sol));
		CHECK(algo.getEvaluations() == 30LL*100);
		CHECK(algo.getSwapAcceptanceRates().empty());
		CHECK(algo.getThreadsUsed() == 1);
	}
	{	// two replicas, more threads than replicas
		PT<solPerm> algo(0.5f, 6.0f, 2, 30, 100, 4, 3, 10);
		algo.start(16, &prob);
		CHECK(algo.getEvaluations() == 2LL*30*100);
		CHECK(algo.getThreadsUsed() == 2);
		CHECK(algo.getSwapAcceptanceRates().size() == 1);
	}
	{	// more replicas than threads; thN = 0 means all cores
		PT<solPerm> algo(0.5f, 6.0f, 33, 5, 60, 4, 0, 50);
		algo.start(2, &prob);
		CHECK(algo.getEvaluations() == 33LL*5*60);
		CHECK(algo.getPTLEnd() == 60);
		algo.start(0, &prob);
		CHECK(algo.getThreadsUsed() >= 1 && algo.getThreadsUsed() <= 33);
	}
	{	// PTL = 0: the best initial solution
		PT<solPerm> algo(0.5f, 6.0f, 4, 30, 0, 4, 0, 50);
		solPerm best = algo.start(2, &prob);
		CHECK(isPermutation(best.sol));
		CHECK(algo.getEvaluations() == 0 && algo.getPTLEnd() == 0);
	}
	{	// MCL = 1
		PT<solPerm> algo(0.5f, 6.0f, 4, 1, 500, 4, 0, 50);
		algo.start(4, &prob);
		CHECK(algo.getEvaluations() == 4LL*500);
	}
}

TEST(dynamic_temperature_updates){
	PermS prob(30);
	{
		PT<solPerm> algo(0.5f, 6.0f, 8, 40, 400, 4, 0, 50);
		algo.start(4, &prob);
		CHECK(algo.getTemperatures() == algo.getInitialTemperatures());   // type 0: fixed
	}
	for(int upType = 1; upType <= 6; upType++){
		PT<solPerm> algo(0.5f, 6.0f, 8, 40, 400, 4, upType, 50);
		algo.setSeed(11);
		solPerm best = algo.start(4, &prob);
		CHECK(isPermutation(best.sol));
		CHECK(algo.getEvaluations() == 8LL*40*400);
		CHECK(algo.getTemperatures().size() == 8);
		bool changed = false;
		for(size_t i = 0; i < 8; i++){
			double t = algo.getTemperatures()[i];
			CHECK(std::isfinite(t) && t > 0.0);
			changed = changed || (t != algo.getInitialTemperatures()[i]);
		}
		CHECK(changed);
	}
}

TEST(user_defined_temperatures_and_initial_solutions){
	PermS prob(20);
	PT<solPerm> algo(0.5f, 6.0f, 8, 30, 50, 4, 0, 50);
	algo.setTemperatures({1.0, 2.0, 4.0});
	algo.start(4, &prob);
	CHECK(algo.getInitialTemperatures() == std::vector<double>({1.0, 2.0, 4.0}));
	CHECK(algo.getReplicas().size() == 3);

	// warm start with the optimum, not evaluated: the API evaluates it
	solPerm opt;
	opt.sol.resize(20);
	std::iota(opt.sol.begin(), opt.sol.end(), 0);
	opt.evalSol = 12345;
	algo.setInitialSolutions({opt});
	solPerm best = algo.start(4, &prob);
	CHECK_NEAR(best.evalSol, 0.0, 0);
	CHECK(algo.getBestPTLSolIndex() == 0);

	// initial solutions also work when the class uses the v1 form of evaluate
	PermD probD(20);
	PT<solPerm> algoD(0.5f, 6.0f, 4, 30, 50, 4, 0, 50);
	algoD.setInitialSolutions({opt});
	CHECK_NEAR(algoD.start(4, &probD).evalSol, 0.0, 0);
}

TEST(real_valued_objective_and_v1_compatibility){
	// All the energies are below 1. v1 truncated them to int in the exchange test,
	// so every exchange was accepted. v2 must use the real values.
	Sphere prob(3);
	PT<solReal> algo(0.001f, 0.2f, 6, 20, 500, 4, 0, 50);
	algo.setSeed(5);
	algo.start(4, &prob);
	bool rejected = false;
	for(double r : algo.getSwapAcceptanceRates()) rejected = rejected || (r < 0.95);
	CHECK(rejected);

	algo.setV1Compatibility(true);
	algo.start(4, &prob);
	for(double r : algo.getSwapAcceptanceRates()) CHECK_NEAR(r, 1.0, 0);

	// the compatibility mode also runs with the dynamic updates
	PermS perm(30);
	PT<solPerm> legacy(0.5f, 6.0f, 8, 40, 400, 4, 1, 50);
	legacy.setV1Compatibility(true);
	CHECK(isPermutation(legacy.start(4, &perm).sol));
}

TEST(sampled_distribution_is_boltzmann){
	// Statistical check of the acceptance and exchange rules: after many steps the
	// state at each temperature must follow exp(-E/T)/Z.
	// Checked for the two exchange schemes.
	for(int scheme : {PT_SWAP_SEQUENTIAL, PT_SWAP_EVEN_ODD}){
	const std::vector<double> temps = {0.5, 1.0, 2.0};
	const int runs = 4000;
	std::vector<std::vector<int>> count(temps.size(), std::vector<int>(6, 0));
	Ring prob;
	PT<solState> algo(1, 2, 3, 4, 60, 1, 0, 50);
	algo.setTemperatures(temps);
	algo.setSwapScheme(scheme);
	for(int r = 0; r < runs; r++){
		algo.setSeed(1000 + r);
		algo.start(1, &prob);
		for(size_t k = 0; k < temps.size(); k++) count[k][algo.getReplicas()[k].state]++;
	}
	for(size_t k = 0; k < temps.size(); k++){
		double Z = 0;
		for(double e : kRingEnergy) Z += std::exp(-e/temps[k]);
		for(int s = 0; s < 6; s++){
			double p = std::exp(-kRingEnergy[s]/temps[k])/Z;
			double tol = 4.0*std::sqrt(p*(1-p)/runs) + 0.003;
			CHECK_NEAR((double)count[k][s]/runs, p, tol);
		}
	}
	}
}

TEST(exchange_schemes_and_round_trips){
	// small problem and hot ladder: the replicas travel between the two ends
	PermS prob(8);
	for(int scheme : {PT_SWAP_SEQUENTIAL, PT_SWAP_EVEN_ODD}){
		PT<solPerm> algo(1.0f, 12.0f, 6, 10, 2000, 4, 0, 50);
		algo.setSwapScheme(scheme);
		algo.setSeed(8);
		solPerm one = algo.start(1, &prob);
		long long trips = algo.getRoundTrips();
		std::vector<double> rates = algo.getSwapAcceptanceRates();
		CHECK(trips > 0);
		CHECK(algo.getEvaluations() == 6LL*10*2000);
		double barrier = 0.0;
		for(double r : rates){ CHECK(r > 0.0 && r <= 1.0); barrier += 1.0 - r; }
		CHECK_NEAR(algo.getCommunicationBarrier(), barrier, 1e-12);
		CHECK(barrier >= 0.0 && barrier < 5.0);
		// same seed, any number of threads
		solPerm many = algo.start(8, &prob);
		CHECK(many.sol == one.sol && many.evalSol == one.evalSol);
		CHECK(algo.getRoundTrips() == trips);
		CHECK(algo.getSwapAcceptanceRates() == rates);
		// the dynamic updates also work with the scheme
		PT<solPerm> dyn(1.0f, 12.0f, 6, 10, 600, 4, 1, 10);
		dyn.setSwapScheme(scheme);
		dyn.start(4, &prob);
		CHECK(dyn.getTemperatures().front() == dyn.getInitialTemperatures().front());
		CHECK(dyn.getTemperatures().back() == dyn.getInitialTemperatures().back());
	}
	PT<solPerm> bad(1.0f, 12.0f, 6, 10, 10, 4, 0, 50);
	CHECK_THROWS(bad.setSwapScheme(7), std::invalid_argument);
}

TEST(dynamic_updates_keep_the_end_temperatures){
	PermS prob(30);
	for(int upType = 1; upType <= 6; upType++){
		for(int interval : {2, 20}){
			PT<solPerm> algo(0.5f, 6.0f, 8, 40, 600, 4, upType, interval);
			algo.setSeed(21);
			algo.start(4, &prob);
			std::vector<double> ini = algo.getInitialTemperatures(), t = algo.getTemperatures();
			CHECK(t.front() == ini.front());
			// type 5 leaves the last temperature free, up to the initial one
			if(upType == 5) CHECK(t.back() <= ini.back()); else CHECK(t.back() == ini.back());
			// type 3 is the rule of the thesis: it keeps the ends, the order is not checked
			if(upType == 5) for(size_t i = 1; i < t.size(); i++) CHECK(t[i] >= t[i-1]);
			else if(upType != 3) for(size_t i = 1; i < t.size(); i++) CHECK(t[i] > t[i-1]);
		}
	}
	// two replicas: both temperatures are anchored
	PT<solPerm> two(0.5f, 6.0f, 2, 40, 300, 4, 1, 10);
	two.start(2, &prob);
	CHECK(two.getTemperatures() == two.getInitialTemperatures());
	// the v1 rule (not anchored) is still available
	PT<solPerm> old(0.5f, 6.0f, 8, 40, 600, 4, 1, 20);
	old.setV1Compatibility(true);
	CHECK(isPermutation(old.start(4, &prob).sol));
}

TEST(calibrate_temperatures){
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 40, 200, 4, 0, 50);
	algo.setSeed(5);
	std::pair<double,double> t = algo.calibrateTemperatures(&prob);
	CHECK(t.first > 0.0 && t.second > t.first);
	std::vector<double> ladder = algo.getInitialTemperatures();
	CHECK(ladder.size() == 8);
	CHECK_NEAR(ladder.front(), t.first, 1e-5*t.first);
	CHECK_NEAR(ladder.back(), t.second, 1e-5*t.second);

	// same seed, same result; lower acceptance targets give lower temperatures
	PT<solPerm> other(0.5f, 6.0f, 8, 40, 200, 4, 0, 50);
	other.setSeed(5);
	CHECK(other.calibrateTemperatures(&prob) == t);
	std::pair<double,double> colder = other.calibrateTemperatures(&prob, 5000, 0.001, 0.00001);
	CHECK(colder.first < t.first && colder.second < t.second);

	CHECK_THROWS(algo.calibrateTemperatures(&prob, 5), std::invalid_argument);
	CHECK_THROWS(algo.calibrateTemperatures(&prob, 1000, 0.5, 0.6), std::invalid_argument);
	CHECK_THROWS(algo.calibrateTemperatures(&prob, 1000, 1.0, 0.1), std::invalid_argument);

	// the run uses the calibrated ladder
	solPerm best = algo.start(4, &prob);
	CHECK(isPermutation(best.sol));
	CHECK(algo.getTemperatures() == ladder);
	CHECK(algo.getMoveAcceptanceRates().front() < algo.getMoveAcceptanceRates().back());

	// a problem that only has neighborEval
	PermFast fast(40);
	PT<solPerm> f(0.5f, 6.0f, 4, 40, 50, 4, 0, 50);
	std::pair<double,double> tf = f.calibrateTemperatures(&fast, 1000);
	CHECK(tf.first > 0.0 && tf.second > tf.first);
}

TEST(calibrate_ladder_for_exchange_rate){
	PermS prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 40, 300, 4, 0, 50);
	algo.setSeed(5);
	int callbacks = 0;
	algo.setOnImprove([&](const solPerm&, int, double){ ++callbacks; });
	std::vector<double> ladder = algo.calibrateLadder(&prob, 4);
	CHECK(callbacks == 0);                       // the pilot runs are not reported
	CHECK(ladder.size() == 8);
	CHECK_NEAR(ladder.front(), 0.5, 1e-9);        // the first temperature is kept
	for(size_t i = 1; i < ladder.size(); i++) CHECK(ladder[i] > ladder[i-1]);
	CHECK(algo.getInitialTemperatures() == ladder);

	// same seed, same ladder
	PT<solPerm> other(0.5f, 6.0f, 8, 40, 300, 4, 0, 50);
	other.setSeed(5);
	CHECK(other.calibrateLadder(&prob, 2) == ladder);

	// the configuration of the user is back, and the exchange rates are around the target
	solPerm best = algo.start(4, &prob);
	CHECK(isPermutation(best.sol));
	CHECK(callbacks > 0);
	CHECK(algo.getPTLEnd() == 300);
	CHECK(algo.getEvaluations() == 8LL*40*300);
	CHECK(algo.getTemperatures() == ladder);
	double mean = 0.0;
	for(double r : algo.getSwapAcceptanceRates()) mean += r/7.0;
	std::printf("    exchange rates with the calibrated ladder: mean %.3f (target 0.23)\n", mean);
	CHECK(mean > 0.13 && mean < 0.36);

	// every pair, not only the mean: a longer run with the calibrated ladder
	PT<solPerm> longRun(0.5f, 6.0f, 8, 40, 3000, 4, 0, 50);
	longRun.setTemperatures(ladder);
	longRun.setSeed(6);
	longRun.start(4, &prob);
	double worst = 0.0;
	for(double r : longRun.getSwapAcceptanceRates()) worst = std::max(worst, std::fabs(r - 0.23));
	std::printf("    largest deviation of a pair from the target in 3000 rounds: %.3f\n", worst);
	CHECK(worst < 0.10);

	// a higher target rate needs closer temperatures
	std::vector<double> closer = other.calibrateLadder(&prob, 4, 0.60);
	CHECK(closer.back() < ladder.back());
	// first temperature given in the call; problem that only has neighborEval
	PermFast fast(40);
	PT<solPerm> f(0.5f, 6.0f, 4, 40, 50, 4, 2, 10);
	std::vector<double> lf = f.calibrateLadder(&fast, 2, 0.23, 50, 3, 1.0);
	CHECK(lf.size() == 4);
	CHECK_NEAR(lf.front(), 1.0, 1e-9);
	CHECK(isPermutation(f.start(2, &fast).sol));

	CHECK_THROWS(algo.calibrateLadder(&prob, 4, 0.0), std::invalid_argument);
	CHECK_THROWS(algo.calibrateLadder(&prob, 4, 1.0), std::invalid_argument);
	CHECK_THROWS(algo.calibrateLadder(&prob, 4, 0.23, 5), std::invalid_argument);
}

TEST(several_moves_with_probabilities){
	PermMoves prob(40);
	CHECK(prob.getMoveCount() == 2);
	CHECK_NEAR(prob.getMoveProbabilities()[0], 0.75, 1e-12);
	PT<solPerm> algo(0.5f, 6.0f, 8, 50, 400, 4, 0, 50);
	algo.setSeed(9);
	solPerm best = algo.start(4, &prob);
	CHECK(isPermutation(best.sol));
	const std::vector<PTMoveStat>& st = algo.getMoveStatistics();
	CHECK(st.size() == 2);
	if(st.size() == 2){
		CHECK(st[0].proposed + st[1].proposed == algo.getEvaluations());
		CHECK_NEAR((double)st[0].proposed/(double)algo.getEvaluations(), 0.75, 0.01);
		CHECK(st[0].accepted > 0 && st[0].accepted <= st[0].proposed);
		CHECK(st[1].accepted > 0 && st[1].accepted <= st[1].proposed);
		double rate = 0.0;
		for(double r : algo.getMoveAcceptanceRates()) rate += r/8.0;
		CHECK_NEAR((double)(st[0].accepted + st[1].accepted)/(double)algo.getEvaluations(), rate, 1e-9);
	}

	CHECK_THROWS(prob.setMoveProbabilities({}), std::invalid_argument);
	CHECK_THROWS(prob.setMoveProbabilities({0.0, 0.0}), std::invalid_argument);
	CHECK_THROWS(prob.setMoveProbabilities({1.0, -1.0}), std::invalid_argument);

	// a move with probability zero is never chosen
	prob.setMoveProbabilities({0.0, 1.0});
	algo.start(4, &prob);
	CHECK(algo.getMoveStatistics()[0].proposed == 0);
	CHECK(algo.getMoveStatistics()[1].proposed == algo.getEvaluations());

	// problems without registered moves: no statistics, chooseMove is an error
	PermS plain(20);
	PT<solPerm> p(0.5f, 6.0f, 4, 20, 20, 4, 0, 50);
	p.start(2, &plain);
	CHECK(p.getMoveStatistics().empty());
	CHECK_THROWS(plain.chooseMove(), std::logic_error);
}

TEST(checkpoint_and_resume){
	const std::string file = "test_ptapi_checkpoint.tmp";
	std::remove(file.c_str());
	PermS prob(60);

	// first part: stopped by a target value, state saved at the end
	PT<solPerm> first(0.5f, 6.0f, 8, 50, 400, 4, 2, 20);
	CHECK(!first.loadCheckpoint(file, readPerm));   // no file
	first.setSeed(3);
	first.setTargetValue(150.0);
	first.setCheckpoint(file, writePerm);
	solPerm part = first.start(4, &prob);
	CHECK(part.evalSol <= 150.0);
	CHECK(first.getPTLEnd() > 0 && first.getPTLEnd() < 400);

	// second part: continues with the same replicas, temperatures, best and round
	PT<solPerm> second(0.5f, 6.0f, 8, 50, 400, 4, 2, 20);
	CHECK(second.loadCheckpoint(file, readPerm));
	solPerm end = second.start(4, &prob);
	CHECK(isPermutation(end.sol));
	CHECK_NEAR(end.evalSol, permCost(end.sol), 1e-9);
	CHECK(end.evalSol <= part.evalSol);
	CHECK(second.getPTLEnd() == 400);
	CHECK(second.getInitialTemperatures() == first.getTemperatures());
	CHECK(second.getEvaluations() >= 8LL*50*(400 - first.getPTLEnd()));
	CHECK(second.getEvaluations() < 8LL*50*400);
	CHECK(second.getElapsedMs() >= first.getElapsedMs());
	CHECK(second.getBestPTLSolIndex() >= first.getBestPTLSolIndex());

	// a finished execution has nothing left to do
	second.saveCheckpoint(file, writePerm);
	PT<solPerm> third(0.5f, 6.0f, 8, 50, 400, 4, 2, 20);
	CHECK(third.loadCheckpoint(file, readPerm));
	solPerm same = third.start(4, &prob);
	CHECK(same.sol == end.sol && same.evalSol == end.evalSol);
	CHECK(third.getEvaluations() == 0);
	CHECK(third.getPTLEnd() == 400);
	// the checkpoint is used once: the next start() is a new execution
	third.start(4, &prob);
	CHECK(third.getEvaluations() == 8LL*50*400);

	// damaged files are rejected
	std::string content;
	{ std::ifstream in(file.c_str()); std::stringstream ss; ss << in.rdbuf(); content = ss.str(); }
	{ std::ofstream out(file.c_str(), std::ios::trunc); out << content.substr(0, content.size()/2); }
	CHECK_THROWS(third.loadCheckpoint(file, readPerm), std::runtime_error);
	{ std::ofstream out(file.c_str(), std::ios::trunc); out << "this is not a checkpoint\n"; }
	CHECK_THROWS(third.loadCheckpoint(file, readPerm), std::runtime_error);
	PT<solPerm> none(0.5f, 6.0f, 8, 50, 400, 4, 0, 20);
	CHECK_THROWS(none.saveCheckpoint(file, writePerm), std::logic_error);   // nothing was executed
	std::remove(file.c_str());

	// periodic checkpoint while it runs, and time limit of a resumed execution
	Slow slow(40);
	PT<solPerm> timed(0.5f, 6.0f, 8, 20, 2000000000, 4, 0, 50);
	timed.setTimeLimit(0.8);
	timed.setCheckpoint(file, writePerm, 0.1);
	bool seenWhileRunning = false;
	std::thread watcher([&]{
		std::this_thread::sleep_for(std::chrono::milliseconds(600));
		PT<solPerm> reader(0.5f, 6.0f, 8, 20, 2000000000, 4, 0, 50);
		seenWhileRunning = reader.loadCheckpoint(file, readPerm);
	});
	solPerm tb = timed.start(4, &slow);
	watcher.join();
	CHECK(seenWhileRunning);
	PT<solPerm> rest(0.5f, 6.0f, 8, 20, 2000000000, 4, 0, 50);
	rest.setTimeLimit(0.8);
	CHECK(rest.loadCheckpoint(file, readPerm));
	solPerm rb = rest.start(4, &slow);
	CHECK(rest.getEvaluations() == 0);          // the time was already used
	CHECK(rb.evalSol == tb.evalSol);
	rest.setTimeLimit(1.1);
	CHECK(rest.loadCheckpoint(file, readPerm));
	ExecTime et;
	rb = rest.start(4, &slow);
	CHECK(rest.getEvaluations() > 0);           // only the remaining time (about 0.3 s)
	CHECK(et.getTimeMs() < 700);
	CHECK(rb.evalSol <= tb.evalSol);
	std::remove(file.c_str());
	std::remove((file + ".tmp").c_str());
}

// ---------------------------------------------------------------------------
// stop criteria added after the first release candidate
// ---------------------------------------------------------------------------
TEST(stop_by_rounds_without_improvement){
	PermS prob(12);
	const int PTL = 200000;
	PT<solPerm> algo(0.3f, 4.0f, 6, 30, PTL, 4, 0, 50);
	algo.setSeed(4);
	algo.setStagnationRounds(1500);
	algo.start(4, &prob);
	CHECK(algo.getPTLEnd() < PTL);
	CHECK(algo.getPTLEnd() > algo.getBestPTLSolIndex() + 1500);
	CHECK(algo.getPTLEnd() <= algo.getBestPTLSolIndex() + 1500 + 50);
	CHECK(algo.getStopReason() == PT_STOP_STAGNATION);
}

TEST(stop_by_time_without_improvement){
	PermS prob(10);
	PT<solPerm> algo(0.3f, 4.0f, 6, 30, 2000000000, 4, 0, 50);
	algo.setStagnationTime(0.25);
	algo.setTimeLimit(20.0);                          // safety net
	solPerm best = algo.start(4, &prob);
	CHECK(best.evalSol == 0.0);                       // found long before the stop
	CHECK(algo.getStopReason() == PT_STOP_STAGNATION);
	CHECK(algo.getElapsedMs() >= algo.getTimeToBestMs() + 250.0 - 1.0);
	CHECK(algo.getElapsedMs() <= algo.getTimeToBestMs() + 250.0 + 400.0);
}

TEST(stop_reason_is_reported){
	PermS prob(40);
	{
		PT<solPerm> algo(0.5f, 6.0f, 6, 30, 50, 4, 0, 50);
		algo.start(4, &prob);
		CHECK(algo.getStopReason() == PT_STOP_PTL);
		CHECK(std::string(algo.getStopReasonName()) == "PTL");
	}
	{
		PT<solPerm> algo(0.5f, 6.0f, 6, 30, 2000000000, 4, 0, 50);
		algo.setTimeLimit(0.1);
		algo.start(4, &prob);
		CHECK(algo.getStopReason() == PT_STOP_TIME);
	}
	{
		PT<solPerm> algo(0.5f, 6.0f, 6, 30, 2000000000, 4, 0, 50);
		algo.setTargetValue(300.0);
		algo.setTimeLimit(20.0);
		algo.start(4, &prob);
		CHECK(algo.getStopReason() == PT_STOP_TARGET);
	}
	{
		PT<solPerm> algo(0.5f, 6.0f, 6, 30, 2000000000, 4, 0, 50);
		std::thread other([&]{ std::this_thread::sleep_for(std::chrono::milliseconds(80)); algo.requestStop(); });
		algo.start(4, &prob);
		other.join();
		CHECK(algo.getStopReason() == PT_STOP_REQUEST);
	}
	{
		Faulty bad(20);
		PT<solPerm> algo(0.5f, 6.0f, 6, 30, 1000, 4, 0, 50);
		CHECK_THROWS(algo.start(4, &bad), std::runtime_error);
		CHECK(algo.getStopReason() == PT_STOP_ERROR);
	}
}

TEST(stop_condition_and_lower_bound){
	PermS prob(40);
	// any rule written by the user, checked once per round
	PT<solPerm> algo(0.5f, 6.0f, 6, 30, 100000, 4, 0, 50);
	algo.setSeed(2);
	std::atomic<int> calls{0};
	std::atomic<bool> sane{true};
	algo.setStopCondition([&](const PTStatus& st){
		calls.fetch_add(1);
		if(st.round < 1 || st.elapsedMs < 0.0 || st.bestRound > st.round || st.bestMs > st.elapsedMs + 1.0 || !(st.bestValue >= 0.0)) sane.store(false);
		return st.round >= 37;
	});
	solPerm best = algo.start(4, &prob);
	CHECK(algo.getStopReason() == PT_STOP_CONDITION);
	CHECK(calls.load() == 37);
	CHECK(sane.load());
	CHECK(algo.getPTLEnd() >= 37 && algo.getPTLEnd() <= 38);
	CHECK(isPermutation(best.sol));

	// lower bound with a tolerance: stops at bound + gap*|bound|
	PT<solPerm> lb(0.5f, 6.0f, 6, 30, 2000000000, 4, 0, 50);
	lb.setLowerBound(100.0, 0.5);                     // 150
	lb.setTimeLimit(20.0);
	best = lb.start(4, &prob);
	CHECK(best.evalSol <= 150.0);
	CHECK(lb.getStopReason() == PT_STOP_TARGET);

	// an exception thrown by the condition reaches the caller
	PT<solPerm> bad(0.5f, 6.0f, 6, 30, 1000, 4, 0, 50);
	bad.setStopCondition([](const PTStatus& st) -> bool { if(st.round >= 5) throw std::runtime_error("condition failed"); return false; });
	CHECK_THROWS(bad.start(4, &prob), std::runtime_error);
}

// ---------------------------------------------------------------------------
// time of each component
// ---------------------------------------------------------------------------
TEST(time_of_each_component){
	Slow prob(30);
	PT<solPerm> algo(0.5f, 6.0f, 4, 20, 40, 4, 1, 5);
	algo.setSeed(6);
	solPerm plain = algo.start(4, &prob);
	const PTTiming off = algo.getTiming();
	CHECK(off.wallMs > 0.0 && off.constructionMs >= 0.0 && off.constructionMs <= off.wallMs);
	CHECK(off.chainsMs == 0.0 && off.exchangesMs == 0.0 && off.updatesMs == 0.0 && off.waitingMs == 0.0);   // not measured
	CHECK(off.threads == 4);

	algo.setProfiling(true);
	solPerm measured = algo.start(4, &prob);
	const PTTiming t = algo.getTiming();
	CHECK(measured.sol == plain.sol);                 // measuring does not change the run
	CHECK(t.chainsMs > 0.0 && t.exchangesMs > 0.0 && t.updatesMs > 0.0 && t.waitingMs >= 0.0);
	// the parts add up to the time of the threads after the construction
	const double threadMs = t.threads*(t.wallMs - t.constructionMs);
	CHECK(t.chainsMs + t.exchangesMs + t.updatesMs <= threadMs*1.02 + 1.0);
	CHECK_NEAR(t.chainsMs + t.exchangesMs + t.updatesMs + t.waitingMs, threadMs, threadMs*0.02 + 1.0);
	// an expensive evaluation: nearly all the busy time is in the chains
	CHECK(t.chainsMs > 20.0*(t.exchangesMs + t.updatesMs));
	std::printf("    wall %.1f ms, construction %.2f ms, chains %.1f, exchanges %.3f, updates %.3f, waiting %.1f (thread-ms, %d threads)\n",
		t.wallMs, t.constructionMs, t.chainsMs, t.exchangesMs, t.updatesMs, t.waitingMs, t.threads);
}

// ---------------------------------------------------------------------------
// initial solutions
// ---------------------------------------------------------------------------
// A constructive heuristic for the first replica, random solutions for the others
class PermHeuristic: public PermBase<PermHeuristic>{
	public:
		std::mutex m;
		std::set<int> seen;
		int countSeen = -1;
		int sleepMs = 0;
		int failAt = -1;
		using PermBase::PermBase;
		double evaluate(const solPerm& s){ return permCost(s.sol); }
		solPerm constructionReplica(int index, int count) override {
			{ std::lock_guard<std::mutex> lock(m); seen.insert(index); countSeen = count; }
			if(sleepMs) std::this_thread::sleep_for(std::chrono::milliseconds(sleepMs));
			if(index == failAt) throw std::runtime_error("construction failed");
			if(index != 0) return construction();
			solPerm s;
			s.sol.resize(n);
			std::iota(s.sol.begin(), s.sol.end(), 0);   // the optimum
			s.evalSol = 0.0;
			return s;
		}
};

TEST(construction_per_replica_and_in_parallel){
	// the problem is told which replica it is building
	PermHeuristic prob(40);
	PT<solPerm> algo(0.5f, 6.0f, 8, 30, 20, 4, 0, 50);
	algo.setSeed(9);
	solPerm best = algo.start(4, &prob);
	CHECK(prob.seen.size() == 8 && *prob.seen.begin() == 0 && *prob.seen.rbegin() == 7);
	CHECK(prob.countSeen == 8);
	CHECK(best.evalSol == 0.0 && algo.getBestPTLSolIndex() == 0);

	// built by several threads: same solutions (each replica has its own random engine)
	PermS plain(40);
	PT<solPerm> seq(0.5f, 6.0f, 8, 30, 100, 4, 0, 50), par(0.5f, 6.0f, 8, 30, 100, 4, 0, 50);
	seq.setSeed(31); par.setSeed(31);
	par.setParallelConstruction(true);
	solPerm a = seq.start(4, &plain), b = par.start(4, &plain);
	CHECK(a.sol == b.sol && a.evalSol == b.evalSol);
	CHECK(seq.getReplicas().size() == par.getReplicas().size());
	for(size_t i = 0; i < seq.getReplicas().size(); ++i) CHECK(seq.getReplicas()[i].sol == par.getReplicas()[i].sol);

	// and it takes less time when the construction is expensive
	PermHeuristic slow(40);
	slow.sleepMs = 30;
	PT<solPerm> s1(0.5f, 6.0f, 8, 30, 5, 4, 0, 50), s2(0.5f, 6.0f, 8, 30, 5, 4, 0, 50);
	s2.setParallelConstruction(true);
	s1.start(8, &slow);
	s2.start(8, &slow);
	CHECK(s1.getTiming().constructionMs >= 8*30.0 - 1.0);
	CHECK(s2.getTiming().constructionMs < 0.5*s1.getTiming().constructionMs);
	std::printf("    construction of 8 replicas (30 ms each): sequential %.0f ms, parallel %.0f ms\n", s1.getTiming().constructionMs, s2.getTiming().constructionMs);

	// an exception in the construction reaches the caller in both modes
	PermHeuristic bad(40);
	bad.failAt = 5;
	CHECK_THROWS(s1.start(8, &bad), std::runtime_error);
	CHECK_THROWS(s2.start(8, &bad), std::runtime_error);
}

// ---------------------------------------------------------------------------
// automatic interval between temperature updates
// ---------------------------------------------------------------------------
TEST(automatic_update_interval){
	PermS prob(40);
	for(int upType = 1; upType <= 6; ++upType){
		PT<solPerm> algo(0.5f, 8.0f, 8, 30, 600, 1, upType, PT_UPDATE_AUTO);   // linear ladder: a bad one
		algo.setSeed(12);
		solPerm best = algo.start(4, &prob);
		CHECK(isPermutation(best.sol));
		const std::vector<double> ini = algo.getInitialTemperatures();
		const std::vector<double>& fin = algo.getTemperatures();
		CHECK(fin.size() == ini.size());
		bool changed = false;
		for(size_t i = 0; i < fin.size(); ++i){
			CHECK(std::isfinite(fin[i]) && fin[i] > 0.0);
			if(std::fabs(fin[i] - ini[i]) > 1e-9) changed = true;
		}
		// the updates did happen (the rule of type 3 only accepts a new temperature
		// between its neighbours, which may not occur in such a short run)
		if(upType != 3) CHECK(changed);
		CHECK_NEAR(fin.front(), ini.front(), 1e-9);
		if(upType != 3 && upType != 5) CHECK_NEAR(fin.back(), ini.back(), 1e-9);
		if(upType == 5) CHECK(fin.back() <= ini.back() + 1e-9);
	}
	// still reproducible for any number of threads (types that are with a fixed interval)
	for(int upType : {1, 2, 5}){
		std::vector<int> ref;
		for(int th : {1, 3, 8}){
			PT<solPerm> algo(0.5f, 8.0f, 8, 30, 300, 4, upType, PT_UPDATE_AUTO);
			algo.setSeed(77);
			solPerm best = algo.start(th, &prob);
			if(ref.empty()) ref = best.sol;
			CHECK(best.sol == ref);
		}
	}
	// the first interval and the growth can be chosen
	PT<solPerm> algo(0.5f, 8.0f, 8, 30, 300, 1, 1, PT_UPDATE_AUTO);
	algo.setAutoUpdate(4, 1.5);
	algo.setSeed(1);
	CHECK(isPermutation(algo.start(4, &prob).sol));
	CHECK_THROWS(algo.setAutoUpdate(1, 2.0), std::invalid_argument);
	CHECK_THROWS(algo.setAutoUpdate(10, 0.5), std::invalid_argument);
}

// ---------------------------------------------------------------------------
// experiments and command line
// ---------------------------------------------------------------------------
TEST(experiment_runs_and_reports){
	const std::string dir = "test_ptapi_experiment.tmp";
	std::filesystem::remove_all(dir);
	std::atomic<int> calls{0};
	std::vector<unsigned long long> seeds;
	auto runner = [&](const std::string& instance, unsigned long long seed){
		calls.fetch_add(1);
		seeds.push_back(seed);
		PermS prob(std::stoi(instance));              // the "instance" is the size
		PT<solPerm> algo(0.5f, 6.0f, 6, 20, 15, 4, 0, 50);
		algo.setSeed(seed);
		solPerm best = algo.start(4, &prob);
		PTRunResult r = ptResult(algo, best);
		std::ostringstream text; writePerm(text, best);
		r.solution = text.str();
		return r;
	};

	PTExperiment exp(dir);
	exp.setVerbose(false);
	exp.addInstance("30", "small");
	exp.setReference(1.0e9);                          // every run "reaches" it
	exp.addInstance("60", "large");
	exp.setRuns(5);
	exp.setFirstSeed(100);
	exp.run(runner);
	CHECK(calls.load() == 10);
	CHECK(seeds.size() == 10 && seeds[0] == 100 && seeds[4] == 104 && seeds[5] == 100);
	CHECK_THROWS(exp.addInstance("30", "small"), std::invalid_argument);

	// statistics recomputed from the individual results
	const std::vector<PTInstanceSummary>& sum = exp.getSummary();
	CHECK(sum.size() == 2 && sum[0].name == "small" && sum[1].name == "large");
	for(size_t i = 0; i < sum.size(); ++i){
		std::vector<double> v;
		for(const PTRunResult& r : exp.getRuns(i)) v.push_back(r.value);
		CHECK(v.size() == 5 && sum[i].runs == 5);
		std::sort(v.begin(), v.end());
		double mean = 0, ss = 0;
		for(double x : v) mean += x/5.0;
		for(double x : v) ss += (x - mean)*(x - mean);
		CHECK(sum[i].best == v[0] && sum[i].worst == v[4] && sum[i].median == v[2]);
		CHECK_NEAR(sum[i].mean, mean, 1e-9);
		CHECK_NEAR(sum[i].sd, std::sqrt(ss/4.0), 1e-9);
		CHECK(sum[i].meanSeconds > 0.0 && sum[i].meanEvaluations == 6.0*20*15);
	}
	CHECK(sum[0].hasReference && sum[0].hits == 5 && sum[0].gapBest < 0.0);
	CHECK(!sum[1].hasReference);

	// files
	for(const char* f : {"/runs.csv", "/summary.csv", "/summary.txt", "/summary.tex", "/small/run_1.txt", "/small/run_5.txt", "/large/run_5.txt"})
		CHECK(std::filesystem::exists(dir + f));
	{
		std::ifstream in((dir + "/runs.csv").c_str());
		std::string line; int lines = 0;
		while(std::getline(in, line)) ++lines;
		CHECK(lines == 11);                             // header + 10 executions
	}

	// an interrupted experiment continues: nothing that is on disk runs again
	std::filesystem::remove(dir + "/large/run_4.txt");
	PTExperiment again(dir);
	again.setVerbose(false);
	again.addInstance("30", "small");
	again.addInstance("60", "large");
	again.setRuns(5);
	again.setFirstSeed(100);
	calls = 0; seeds.clear();
	again.run(runner);
	CHECK(calls.load() == 1 && seeds.size() == 1 && seeds[0] == 103);
	for(size_t i = 0; i < 2; ++i){
		CHECK(again.getSummary()[i].best == sum[i].best && again.getSummary()[i].median == sum[i].median);
		CHECK_NEAR(again.getSummary()[i].mean, sum[i].mean, 1e-9);
		for(size_t k = 0; k < 5; ++k){
			CHECK(again.getRuns(i)[k].value == exp.getRuns(i)[k].value);
			CHECK(again.getRuns(i)[k].seed == exp.getRuns(i)[k].seed);
			CHECK(again.getRuns(i)[k].solution == exp.getRuns(i)[k].solution);
		}
	}
	// unless everything is asked to run again
	again.setResume(false);
	calls = 0;
	again.run(runner);
	CHECK(calls.load() == 10);

	// columns of runs.csv chosen by the user, values of the user, and the same after a resume
	std::filesystem::remove_all(dir);
	auto richRunner = [&](const std::string& instance, unsigned long long seed){
		PTRunResult r = runner(instance, seed);
		r.add("size", std::stod(instance));
		r.add("note", "a,b");                           // needs quoting
		return r;
	};
	auto firstLines = [&](std::string& header, std::string& row){
		std::ifstream in((dir + "/runs.csv").c_str());
		std::getline(in, header); std::getline(in, row);
	};
	for(const std::string& c : PTExperiment::csvColumnNames()) CHECK(!c.empty());
	std::string header1, row1, header2, row2;
	for(int pass = 0; pass < 2; ++pass){
		PTExperiment custom(dir);
		custom.setVerbose(false);
		custom.addInstance("30", "small");
		custom.setReference(10.0);
		custom.setRuns(2);
		custom.setCsvColumns({"instance", "seed", "value", "rounds", "stop", "improvements", "replicas", "temp_first", "temp_last", "barrier", "swap_rate_mean", "gap_percent", "size"});
		calls = 0;
		custom.run(richRunner);
		CHECK(calls.load() == (pass == 0 ? 2 : 0));      // second pass: read from the files
		firstLines(pass == 0 ? header1 : header2, pass == 0 ? row1 : row2);
	}
	CHECK(header1 == "instance,seed,value,rounds,stop,improvements,replicas,temp_first,temp_last,barrier,swap_rate_mean,gap_percent,size,note");
	CHECK(header1 == header2 && row1 == row2);          // nothing is lost when the run is read back
	CHECK(row1.rfind("small,1,", 0) == 0);              // instance and seed
	CHECK(row1.find(",15,PTL,") != std::string::npos);  // rounds and stop reason
	CHECK(row1.find(",6,0.5,6,") != std::string::npos); // replicas and the ends of the ladder
	CHECK(row1.size() > 10 && row1.substr(row1.size() - 9) == ",30,\"a,b\"");
	{
		PTExperiment custom(dir);
		CHECK_THROWS(custom.setCsvColumns({}), std::invalid_argument);
	}

	PTExperiment empty(dir);
	CHECK_THROWS(empty.run(runner), std::logic_error);
	CHECK_THROWS(empty.setRuns(0), std::invalid_argument);
	std::filesystem::remove_all(dir);
}

TEST(experiment_curves_time_to_target_and_parts){
	const std::string dir = "test_ptapi_curves.tmp";
	std::filesystem::remove_all(dir);
	std::atomic<int> calls{0};
	auto runner = [&](const std::string& instance, unsigned long long seed){
		calls.fetch_add(1);
		PermMoves prob(std::stoi(instance));          // two registered moves
		PT<solPerm> algo(0.5f, 6.0f, 6, 20, 40, 4, 0, 50);
		algo.setSeed(seed);
		solPerm best = algo.start(4, &prob);
		return ptResult(algo, best);
	};
	auto describe = [&](PTExperiment& e){
		e.setVerbose(false);
		e.setRuns(5);
		e.addInstance("30", "small");
		e.setReference(1.0e9);                        // every execution reaches it, with its first solution
		e.addInstance("60", "large");
		e.setCsvColumns({"instance", "run", "value", "temperatures_initial", "moves_proposed", "moves_accepted", "moves_rate", "seconds_to_reference"});
	};
	auto lines = [](const std::string& file){
		std::vector<std::vector<double>> rows;
		std::ifstream in(file.c_str());
		std::string line;
		std::getline(in, line);                       // header
		while(std::getline(in, line)){
			for(char& c : line) if(c == ',') c = ' ';
			std::istringstream is(line);
			std::vector<double> v; double x;
			while(is >> x) v.push_back(x);
			rows.push_back(v);
		}
		return rows;
	};

	// two parts: each one runs only its executions, none is repeated or forgotten
	PTExperiment part0(dir), part1(dir), whole(dir);
	describe(part0); part0.setShard(0, 2); part0.run(runner);
	CHECK(calls.load() == 5);
	CHECK(part0.getSummary()[0].runs + part0.getSummary()[1].runs == 5);
	describe(part1); part1.setShard(1, 2); part1.run(runner);
	CHECK(calls.load() == 10);
	CHECK(part1.getSummary()[0].runs == 5 && part1.getSummary()[1].runs == 5);   // it also finds those of the other part
	CHECK_THROWS(part0.setShard(2, 2), std::invalid_argument);

	// report() runs nothing and has everything, convergence curves included
	describe(whole); whole.report();
	CHECK(calls.load() == 10);
	const std::vector<PTInstanceSummary>& sum = whole.getSummary();
	CHECK(sum[0].runs == 5 && sum[1].runs == 5);
	for(size_t i = 0; i < 2; ++i) for(size_t k = 0; k < 5; ++k){
		const PTRunResult& r = whole.getRuns(i)[k];
		CHECK(!r.convergence.empty() && (int)r.convergence.size() == r.improvements);
		CHECK(r.convergence.back().value == r.value);
		CHECK(r.initialTemperatures.size() == 6 && r.initialTemperatures.front() == 0.5);
		CHECK(r.moves.size() == 2);
		if(r.moves.size() == 2) CHECK(r.moves[0].proposed + r.moves[1].proposed == 6LL*20*40 && r.moves[0].accepted > 0 && r.moves[0].accepted <= r.moves[0].proposed);
	}
	CHECK(sum[0].hits == 5 && sum[0].meanSecondsToReference >= 0.0 && sum[0].meanSecondsToReference <= sum[0].meanSecondsToBest);

	// convergence of an instance: by round and by time, never worse as it advances
	for(const char* f : {"/large/convergence.csv", "/large/convergence_time.csv", "/small/convergence.csv"}){
		const std::vector<std::vector<double>> rows = lines(dir + f);
		CHECK(rows.size() >= 5);
		for(size_t i = 0; i < rows.size(); ++i){
			CHECK(rows[i].size() == 6);
			const std::vector<double>& v = rows[i];                 // x, runs, best, median, mean, worst
			CHECK(v[1] >= 1 && v[1] <= 5 && v[2] <= v[3] && v[3] <= v[5] && v[2] <= v[4] && v[4] <= v[5]);
			if(i && rows[i-1][1] == 5){ CHECK(v[0] > rows[i-1][0]); CHECK(v[2] <= rows[i-1][2] && v[5] <= rows[i-1][5] && v[3] <= rows[i-1][3]); }
		}
		const std::vector<double>& last = rows.back();
		const PTInstanceSummary& m = (std::string(f).find("large") != std::string::npos) ? sum[1] : sum[0];
		CHECK(last[1] == 5 && last[2] == m.best && last[5] == m.worst && last[3] == m.median);
		CHECK_NEAR(last[4], m.mean, 1e-6*std::fabs(m.mean));
	}
	// scripts of the figures: one convergence figure per instance, time to target where there is a reference
	for(const char* f : {"/figures.gp", "/figures.tex"}){
		std::ifstream in((dir + f).c_str());
		const std::string text((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
		CHECK(text.find("small/convergence.csv") != std::string::npos && text.find("large/convergence.csv") != std::string::npos);
		CHECK(text.find("small/time_to_target.csv") != std::string::npos && text.find("large/time_to_target.csv") == std::string::npos);
	}
	// time to target: only for the instance with a reference; ascending times and probabilities
	CHECK(!std::filesystem::exists(dir + "/large/time_to_target.csv"));
	{
		const std::vector<std::vector<double>> rows = lines(dir + "/small/time_to_target.csv");
		CHECK(rows.size() == 5);
		for(size_t i = 0; i < rows.size(); ++i){
			CHECK(rows[i].size() == 3);
			CHECK_NEAR(rows[i][1], ((double)i + 0.5)/5.0, 1e-9);
			if(i) CHECK(rows[i][0] >= rows[i-1][0]);
		}
	}
	{	// the new columns
		std::ifstream in((dir + "/runs.csv").c_str());
		std::string header, row;
		std::getline(in, header); std::getline(in, row);
		CHECK(header == "instance,run,value,temperatures_initial,moves_proposed,moves_accepted,moves_rate,seconds_to_reference");
		CHECK(row.rfind("small,1,", 0) == 0 && row.find(",0.5 ") != std::string::npos);
		const PTRunResult& r = whole.getRuns(0)[0];
		if(r.moves.size() == 2) CHECK(row.find("," + std::to_string(r.moves[0].proposed) + " " + std::to_string(r.moves[1].proposed) + "," + std::to_string(r.moves[0].accepted) + " ") != std::string::npos);
	}
	std::filesystem::remove_all(dir);
}

TEST(statistical_tests){
	// Mann-Whitney U, exact: 3 against 3 without overlap -> 2/C(6,3)
	CHECK_NEAR(ptstat::mannWhitney({1, 2, 3}, {4, 5, 6}), 0.1, 1e-12);
	CHECK_NEAR(ptstat::mannWhitney({4, 5, 6}, {1, 2, 3}), 0.1, 1e-12);
	CHECK_NEAR(ptstat::mannWhitney({1, 2, 3}, {1, 2, 3}), 1.0, 1e-12);
	CHECK_NEAR(ptstat::mannWhitney({5, 5, 5}, {5, 5}), 1.0, 1e-12);                 // every value tied
	CHECK(std::isnan(ptstat::mannWhitney({}, {1, 2})));
	{	// with ties, against the enumeration of every division of the values
		const std::vector<double> a{1, 2, 2, 3}, b{2, 3, 4, 4, 5};
		std::vector<double> all(a); all.insert(all.end(), b.begin(), b.end());
		const std::vector<double> r = ptstat::ranks(all);
		CHECK_NEAR(r[1], 3.0, 1e-12);                                               // the three 2s share ranks 2, 3, 4
		double obs = 0; for(size_t i = 0; i < a.size(); ++i) obs += r[i];
		const double mean = 4*10/2.0;
		int total = 0, tail = 0;
		for(int mask = 0; mask < (1 << 9); ++mask){
			int bits = 0; double sum = 0;
			for(int i = 0; i < 9; ++i) if(mask >> i & 1){ ++bits; sum += r[i]; }
			if(bits != 4) continue;
			++total;
			if(std::fabs(sum - mean) >= std::fabs(obs - mean) - 1e-9) ++tail;
		}
		CHECK(total == 126);
		CHECK_NEAR(ptstat::mannWhitney(a, b), (double)tail/total, 1e-12);
	}
	{	// large samples: the normal approximation separates shifted samples and not equal ones
		std::vector<double> a, b, c;
		std::mt19937 g(7);
		std::uniform_real_distribution<double> u(0, 1);
		for(int i = 0; i < 100; ++i){ a.push_back(u(g)); b.push_back(u(g) + 0.5); c.push_back(u(g)); }
		CHECK(ptstat::mannWhitney(a, b) < 1e-6);
		CHECK(ptstat::mannWhitney(a, c) > 0.01);
	}

	// Wilcoxon signed-rank, exact: 6 pairs, all in the same direction -> 2/2^6
	CHECK_NEAR(ptstat::wilcoxonSignedRank({10, 12, 9, 14, 11, 13}, {11, 14, 12, 18, 16, 19}), 2.0/64.0, 1e-12);
	CHECK_NEAR(ptstat::wilcoxonSignedRank({1, 2, 3}, {1, 2, 3}), 1.0, 1e-12);       // no difference
	CHECK_NEAR(ptstat::wilcoxonSignedRank({1, 5, 3, 3}, {2, 3, 3, 6}), 0.75, 1e-12); // |d| = 1, 2, (0), 3: W+ = 2 of 6
	CHECK_THROWS(ptstat::wilcoxonSignedRank({1, 2}, {1}), std::invalid_argument);

	// chi-square tail and Friedman (values of scipy.stats)
	CHECK_NEAR(ptstat::chiSquareUpper(3.841458820694124, 1), 0.05, 1e-10);
	CHECK_NEAR(ptstat::chiSquareUpper(7.5, 3), 0.0575584519726364, 1e-10);
	std::vector<double> mr;
	CHECK_NEAR(ptstat::friedman({{1, 2, 3}, {2, 3, 1}, {1, 3, 2}, {1, 2, 3}, {1, 3, 2}, {2, 2, 3}}, &mr), 0.07049779029112152, 1e-10);
	CHECK(mr.size() == 3);
	CHECK_NEAR(mr[0], (1 + 2 + 1 + 1 + 1 + 1.5)/6.0, 1e-12);
	CHECK(std::isnan(ptstat::friedman({{1, 2, 3}})));

	// Holm: 0.01, 0.04, 0.03 -> 0.03, 0.06, 0.06
	const std::vector<double> h = ptstat::holm({0.01, 0.04, 0.03});
	CHECK_NEAR(h[0], 0.03, 1e-12); CHECK_NEAR(h[1], 0.06, 1e-12); CHECK_NEAR(h[2], 0.06, 1e-12);
}

TEST(comparison_of_experiments){
	const std::string dir = "test_ptapi_comparison.tmp";
	std::filesystem::remove_all(dir);
	// executions with known values: 'b' is 'a' with 10 more in every instance but the
	// last, where it is equal; 'c' is a little better than 'a' everywhere
	auto make = [&](const std::string& name, double shift, double shiftLast){
		PTExperiment e(dir + "/" + name);
		e.setVerbose(false);
		e.setRuns(8);
		for(int i = 1; i <= 7; ++i) e.addInstance(std::to_string(i), "inst" + std::to_string(i));
		e.run([=](const std::string& instance, unsigned long long seed){
			const int i = std::stoi(instance);
			PTRunResult r;
			r.seed = seed;
			r.value = 100.0*i + (double)((seed*7 + (unsigned long long)i*3) % 5) + (i == 7 ? shiftLast : shift);
			return r;
		});
		return e;
	};
	PTExperiment a = make("a", 0.0, 0.0), b = make("b", 10.0, 0.0), c = make("c", -0.25, -0.25);

	PTComparison cmp;
	CHECK_THROWS(cmp.compute(), std::logic_error);       // needs two configurations
	cmp.add("a", a);
	cmp.add("b", b);
	cmp.add("c", c);
	CHECK_THROWS(cmp.add("a", a), std::invalid_argument);
	CHECK(cmp.getInstances().size() == 7);
	const std::vector<std::vector<PTComparisonCell>>& cells = cmp.getCells();
	const std::vector<PTComparisonTotal>& tot = cmp.getTotals();
	for(size_t i = 0; i < 7; ++i){
		CHECK(cells[i][0].runs == 8 && cells[i][0].sign == 0 && std::isnan(cells[i][0].pValue));
		CHECK_NEAR(cells[i][1].mean - cells[i][0].mean, i == 6 ? 0.0 : 10.0, 1e-9);
		CHECK(cells[i][1].sign == (i == 6 ? 0 : -1));                      // worse than the control, except where equal
		if(i < 6) CHECK(cells[i][1].pValue < 0.001);
		else CHECK_NEAR(cells[i][1].pValue, 1.0, 1e-12);
		CHECK(cells[i][2].pValue >= 0.0 && cells[i][2].pValue <= 1.0);
	}
	CHECK(tot[0].name == "a" && tot[1].name == "b" && tot[2].name == "c");
	CHECK(tot[1].wins == 0 && tot[1].ties == 1 && tot[1].losses == 6);
	CHECK_NEAR(tot[1].pValue, 2.0/64.0, 1e-12);          // 6 pairs differ, all in the same direction
	CHECK_NEAR(tot[2].pValue, 2.0/128.0, 1e-12);         // 7 pairs
	CHECK_NEAR(tot[2].pHolm, 4.0/128.0, 1e-12);          // the smallest p of two comparisons is doubled
	CHECK_NEAR(tot[1].pHolm, 4.0/128.0, 1e-12);          // and the other is not below it
	CHECK_NEAR(tot[2].meanRank, 1.0, 1e-12);             // c has the best mean everywhere
	CHECK_NEAR(tot[0].meanRank, (6*2 + 2.5)/7.0, 1e-12); // a: second, tied with b in the last instance
	CHECK(tot[2].bestMean == 7 && tot[0].bestMean == 0);
	CHECK(tot[2].meanGap == 0.0 && tot[1].meanGap > tot[0].meanGap && tot[0].meanGap > 0.0);
	CHECK(cmp.getFriedmanP() < 0.01);
	CHECK(cmp.table().find("Friedman") != std::string::npos);

	cmp.write(dir);
	for(const char* f : {"/comparison.txt", "/comparison.csv", "/comparison_total.csv", "/comparison.tex"}) CHECK(std::filesystem::exists(dir + f));
	{
		std::ifstream in((dir + "/comparison.csv").c_str());
		std::string line; int n = 0;
		while(std::getline(in, line)) ++n;
		CHECK(n == 1 + 7*3);
	}

	// two configurations: no Friedman test; an experiment read back from its folder
	PTExperiment again(dir + "/b");
	again.setVerbose(false);
	again.setRuns(8);
	for(int i = 1; i <= 7; ++i) again.addInstance(std::to_string(i), "inst" + std::to_string(i));
	again.report();
	PTComparison two;
	two.add("a", a);
	two.add("b again", again);
	CHECK(std::isnan(two.getFriedmanP()));
	CHECK(two.getTotals()[1].losses == 6);
	CHECK_NEAR(two.getTotals()[1].pValue, 2.0/64.0, 1e-12);
	std::filesystem::remove_all(dir);
}

TEST(options_from_the_command_line){
	const char* args[] = {"prog", "40", "--TEMP_INIT", "0.5", "--TEMP_FIM", "6", "--N_REPLICAS", "6", "--MCL", "30", "--PTL", "80",
		"--TEMP_DIST", "3", "--TYPE_UPDATE", "5", "--TEMP_UPDATE", "0", "--THREAD_USED", "3", "--SEED", "17", "--STAGNATION_ROUNDS", "1000",
		"--ALPHA", "0.25", "--NAME", "abc", "--IRACE", "2", "--IRACE_COST", "value+time", "--IRACE_TIME_WEIGHT", "100"};
	const int argc = (int)(sizeof(args)/sizeof(args[0]));
	PTOptions opt(argc, const_cast<char**>(args));
	CHECK(opt.instance == "40");
	CHECK(opt.tempMin == 0.5f && opt.tempMax == 6.0f && opt.replicas == 6 && opt.MCL == 30 && opt.PTL == 80);
	CHECK(opt.tempDist == 3 && opt.updateType == 5 && opt.updateInterval == PT_UPDATE_AUTO && opt.threads == 3 && opt.seed == 17);
	CHECK(opt.stagnationRounds == 1000 && opt.tuner == 2);
	// options of the problem
	CHECK(opt.has("--ALPHA") && !opt.has("--BETA"));
	CHECK(opt.get("--ALPHA", 1.0) == 0.25 && opt.get("--BETA", 3) == 3 && opt.get("--NAME", "x") == "abc");

	PermS prob(std::stoi(opt.instance));
	solPerm best = opt.run<solPerm>(&prob);
	CHECK(isPermutation(best.sol));
	CHECK(opt.result.value == best.evalSol && opt.result.seed == 17 && opt.result.rounds == 80 && opt.result.threads == 3);
	CHECK(opt.result.evaluations == 6LL*30*80 && opt.result.temperatures.size() == 6);
	CHECK_NEAR(opt.tunerCostValue(), best.evalSol + 100.0*opt.result.elapsedMs/1000.0, 1e-9);
	// the same options give the same run
	PTOptions twin(argc, const_cast<char**>(args));
	CHECK(twin.run<solPerm>(&prob).sol == best.sol);

	// cost = time to reach a target; a run that misses it is penalized
	const char* timed[] = {"prog", "40", "--N_REPLICAS", "6", "--MCL", "30", "--PTL", "2000000000", "--TEMP_INIT", "0.5", "--TEMP_FIM", "6",
		"--TARGET", "300", "--TIME_LIMIT", "5", "--IRACE", "1", "--IRACE_COST", "time"};
	PTOptions hit((int)(sizeof(timed)/sizeof(timed[0])), const_cast<char**>(timed));
	hit.run<solPerm>(&prob);
	CHECK(hit.result.stopReason == PT_STOP_TARGET);
	CHECK(hit.tunerCostValue() > 0.0 && hit.tunerCostValue() < 5.0);
	hit.target = -1.0;                                // as if the target had been missed
	CHECK(hit.tunerCostValue() == 50.0);

	const char* bad1[] = {"prog", "40", "--MCL", "many"};
	CHECK_THROWS(PTOptions(4, const_cast<char**>(bad1)), std::invalid_argument);
	const char* bad2[] = {"prog", "40", "--MCL"};
	CHECK_THROWS(PTOptions(3, const_cast<char**>(bad2)), std::invalid_argument);
	const char* bad3[] = {"prog", "40", "--IRACE_COST", "time"};
	CHECK_THROWS(PTOptions(4, const_cast<char**>(bad3)), std::invalid_argument);
	const char* bad4[] = {"prog", "40", "41"};
	CHECK_THROWS(PTOptions(3, const_cast<char**>(bad4)), std::invalid_argument);
}

TEST(idle_threads_do_not_burn_cpu){
	// 2 replicas and 8 threads requested: only 2 threads may be busy.
	Slow prob(20);
	PT<solPerm> algo(0.5f, 6.0f, 2, 20, 2000000000, 4, 0, 50);
	algo.setTimeLimit(0.5);
	double c0 = cpuSeconds();
	ExecTime et;
	algo.start(8, &prob);
	double wall = et.getTimeMs()/1000.0, cpu = cpuSeconds() - c0;
	CHECK(cpu/wall < 2.5);
	std::printf("    cpu/wall = %.2f (2 replicas, 8 threads requested)\n", cpu/wall);
}

TEST(throughput_report){
	// Not a pass/fail test: prints the speed of the API itself (cheap problem) so
	// that versions can be compared. See also the bench folder.
	const int rep = 8, MCL = 500, PTL = 400;
	for(int fast = 0; fast <= 1; fast++){
		PermS classic(100); PermFast delta(100);
		PT<solPerm> algo(0.5f, 6.0f, rep, MCL, PTL, 4, 0, 50);
		algo.setSeed(1);
		if(fast) algo.start(rep, &delta); else algo.start(rep, &classic);
		std::printf("    %-34s %6.1f million steps/s\n", fast ? "neighborEval + delta evaluation:" : "neighbor + evaluate:", algo.getEvaluations()/algo.getElapsedMs()/1000.0);
	}
	CHECK(true);
}

// defined in test_second_unit.cpp: the headers must be usable from several .cpp files
double runFromSecondUnit();
TEST(headers_can_be_included_in_several_files){
	CHECK(runFromSecondUnit() >= 0.0);
}

int main(int argc, char** argv){
	int run = 0;
	for(const TestCase& t : registry()){
		if(argc > 1 && !std::strstr(t.name, argv[1])) continue;
		int before = gFailures;
		std::printf("[ RUN  ] %s\n", t.name);
		std::fflush(stdout);
		try{ t.fn(); }
		catch(const std::exception& e){ ++gFailures; std::printf("    FAILED unexpected exception: %s\n", e.what()); }
		catch(...){ ++gFailures; std::printf("    FAILED unexpected exception\n"); }
		std::printf("[ %s ] %s\n", gFailures == before ? " OK " : "FAIL", t.name);
		++run;
	}
	std::printf("\n%d test(s), %d failure(s) - PTAPI %d.%d\n", run, gFailures, PTAPI_VERSION_MAJOR, PTAPI_VERSION_MINOR);
	return gFailures ? 1 : 0;
}
