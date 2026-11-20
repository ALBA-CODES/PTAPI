/*
 * PTAPI v2 - problem shared by the examples of this folder.
 *
 * Random Euclidean TSP with the 2-opt move. Each example*.cpp shows one group
 * of functions of the API; the problem is here so that the examples only hold
 * what is specific to the function they present.
 */
#ifndef __EXAMPLE_TSP_H__
#define __EXAMPLE_TSP_H__

#include <cstdio>
#include <cmath>
#include <vector>
#include <algorithm>
#include <iostream>
#include "../../include/PT.h"

// Step 1: solution encoding
struct solTSP: public solution{
	std::vector<int> sol;   // order of the cities
};

// Step 2: problem dependent components
class TSP: public Problem<solTSP>{
	protected:
		int n;
		std::vector<double> d;
		double dist(int a, int b) const { return d[a*n + b]; }
	public:
		// 'instanceSeed' defines the cities, not the search
		TSP(int n_, unsigned instanceSeed = 12345): n(n_), d(n_*n_){
			std::mt19937 g(instanceSeed);
			std::uniform_real_distribution<double> u(0, 1000);
			std::vector<double> x(n), y(n);
			for(int i = 0; i < n; i++){ x[i] = u(g); y[i] = u(g); }
			for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) d[i*n+j] = std::round(std::hypot(x[i]-x[j], y[i]-y[j]));
		}

		int size() const { return n; }

		solTSP construction(){
			solTSP s;
			for(int i = 0; i < n; i++) s.sol.push_back(i);
			std::shuffle(s.sol.begin(), s.sol.end(), ptapi::rng());   // engine of the replica
			s.evalSol = evaluate(s);
			return s;
		}

		solTSP neighbor(solTSP s){
			int a, b;
			do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);
			if(a > b) std::swap(a, b);
			std::reverse(s.sol.begin()+a, s.sol.begin()+b+1);
			return s;
		}

		double evaluate(const solTSP& s){
			double c = 0;
			for(int i = 0; i < n; i++) c += dist(s.sol[i], s.sol[(i+1)%n]);
			return c;
		}

		// One step without copies and with delta evaluation (optional, see exampleEvaluate.cpp)
		void neighborEval(const solTSP& cur, solTSP& cand){
			int a, b;
			do{ a = ptapi::randInt(0, n-1); b = ptapi::randInt(0, n-1); }while(a == b);
			if(a > b) std::swap(a, b);

			cand.sol = cur.sol;   // reuses the memory of 'cand'
			std::reverse(cand.sol.begin()+a, cand.sol.begin()+b+1);

			if(a == 0 && b == n-1){ cand.evalSol = cur.evalSol; return; }
			int pa = cur.sol[(a+n-1)%n], ca = cur.sol[a], cb = cur.sol[b], nb = cur.sol[(b+1)%n];
			cand.evalSol = cur.evalSol - dist(pa,ca) - dist(cb,nb) + dist(pa,cb) + dist(ca,nb);
		}

		// true if 's' is a tour with the cost it reports
		bool valid(const solTSP& s){
			std::vector<int> c = s.sol;
			std::sort(c.begin(), c.end());
			if((int)c.size() != n) return false;
			for(int i = 0; i < n; i++) if(c[i] != i) return false;
			return std::fabs(evaluate(s) - s.evalSol) < 1e-6;
		}
};

inline void printRow(const char* name, const std::vector<double>& v, const char* fmt = " %7.2f"){
	printf("%s", name);
	for(double x : v) printf(fmt, x);
	printf("\n");
}

#endif
