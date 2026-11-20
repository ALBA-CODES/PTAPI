// Second translation unit of the tests: including the API headers in more than
// one .cpp file must compile and link.
#include "../include/PT.h"
#include "../include/ExecTime.h"
#include "../include/Experiment.h"
#include "../include/Comparison.h"
#include "../include/Options.h"

namespace{
	struct solInt: public solution{ int v = 0; };

	class Parabola: public Problem<solInt>{
		public:
			solInt construction(){ solInt s; s.v = ptapi::randInt(-100, 100); s.evalSol = evaluate(s); return s; }
			solInt neighbor(solInt s){ s.v += ptapi::randInt(0, 1) ? 1 : -1; return s; }
			double evaluate(const solInt& s){ return (double)s.v*s.v; }
	};
}

double runFromSecondUnit(){
	Parabola prob;
	PT<solInt> algo(0.5f, 20.0f, 4, 50, 100, 4, 0, 50);
	algo.setSeed(1);
	return algo.start(2, &prob).evalSol;
}
