/*
 * PTAPI v2 - the forms of the evaluation function and neighborEval().
 *
 * The problem class declares evaluate() in any of three forms; the API finds
 * out which one at compile time. neighborEval() is optional: it replaces
 * neighbor() + evaluate() by one step without copies, in which the cost can be
 * computed from the cost of the current solution (delta evaluation).
 *
 * The four classes below make the same moves, so with the same seed they
 * reach the same solution; only the time changes.
 *
 *   g++ exampleEvaluate.cpp -std=c++17 -O3 -pthread -o exampleEvaluate
 */
#include "tsp.h"

// Form 1: S evaluate(S) - returns the evaluated solution (it may also change it)
class TSPReturnsSolution: public Problem<solTSP>{
	protected:
		TSP base;
	public:
		TSPReturnsSolution(int n): base(n){}
		solTSP construction(){ return base.construction(); }
		solTSP neighbor(solTSP s){ return base.neighbor(s); }
		solTSP evaluate(solTSP s){ s.evalSol = base.evaluate(s); return s; }
};

// Form 2: double evaluate(S) - the form of the examples published with v1
class TSPReturnsValue: public Problem<solTSP>{
	protected:
		TSP base;
	public:
		TSPReturnsValue(int n): base(n){}
		solTSP construction(){ return base.construction(); }
		solTSP neighbor(solTSP s){ return base.neighbor(s); }
		double evaluate(solTSP s){ return base.evaluate(s); }
};

// Form 3: double evaluate(const S&) - no copy of the solution
class TSPConstReference: public Problem<solTSP>{
	protected:
		TSP base;
	public:
		TSPConstReference(int n): base(n){}
		solTSP construction(){ return base.construction(); }
		solTSP neighbor(solTSP s){ return base.neighbor(s); }
		double evaluate(const solTSP& s){ return base.evaluate(s); }
};

// neighborEval(): the class TSP of tsp.h overrides it (copy-free step with
// delta evaluation). evaluate() is then only used by construction().

template<typename P>
static solTSP run(const char* name, P* prob){
	PT<solTSP> algo(1.0f, 60.0f, 12, 300, 500, 4, 0, 50);
	algo.setSeed(7);
	// the pointer must be of the problem class itself (P*, not Problem<solTSP>*):
	// this is how the API sees which evaluate() the class has
	solTSP s = algo.start(0, prob);
	printf("%-32s best %6.0f   %7.1f ms   %5.1f million evaluations/s\n", name, s.evalSol, algo.getElapsedMs(),
		algo.getEvaluations()/algo.getElapsedMs()/1000.0);
	return s;
}

int main(){
	const int n = 200;
	TSPReturnsSolution p1(n);
	TSPReturnsValue p2(n);
	TSPConstReference p3(n);
	TSP p4(n);

	solTSP a = run("S evaluate(S)", &p1);
	solTSP b = run("double evaluate(S)", &p2);
	solTSP c = run("double evaluate(const S&)", &p3);
	solTSP d = run("neighborEval() with delta", &p4);

	bool ok = (a.sol == b.sol) && (b.sol == c.sol) && (c.sol == d.sol);
	printf("%s\n", ok ? "the four forms reached the same solution" : "ERROR: the forms reached different solutions");
	return ok ? 0 : 1;
}
