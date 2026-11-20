/*
 * Benchmark used to compare version 1 (v1/include) with version 2 (include).
 * It only uses the interface common to both versions. See the Makefile of this
 * folder (make bench / make bench-sd): it builds this file once per version, each
 * one next to its own include folder (SD.h includes "../include/Problem.h").
 *
 *   ./bench <instance | tsp:N> <replicas> <MCL> <PTL> <threads> [timeEnd s] [updateType] [updateInterval] [tempMin] [tempMax]
 *
 * Output: wall time, CPU time (user+sys), evaluations, time spent inside the problem
 * code (SD only) and a check that the returned solution is the best one evaluated.
 */
#include <cstdlib>
#include <atomic>
#include <chrono>
#include <random>
#include <cstdio>
#include <sys/resource.h>
#include "SD/SD.h"
#include "include/PT.h"

static inline long long nowNs(){
	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
static std::atomic<long long> gProbNs{0}, gEvals{0};
static std::atomic<double> gMinEval{1e300};
struct TL { long long ns=0, ev=0; ~TL(){ gProbNs+=ns; gEvals+=ev; } };
static thread_local TL tl;
static void noteMin(double v){ double c=gMinEval.load(); while(v<c && !gMinEval.compare_exchange_weak(c,v)){} }

// ---------- cheap synthetic problem: random Euclidean TSP, 2-opt, full O(n) evaluation
struct solTSP: public solution { std::vector<int> sol; };
class TSP: public Problem<solTSP>{
	int n; std::vector<double> d;
	static std::mt19937& rng(){ thread_local std::mt19937 g{std::random_device{}()}; return g; }
public:
	TSP(int n_, unsigned seed):n(n_),d(n_*n_){
		std::mt19937 g(seed); std::uniform_real_distribution<double> u(0,1000);
		std::vector<double> x(n),y(n); for(int i=0;i<n;i++){x[i]=u(g);y[i]=u(g);}
		for(int i=0;i<n;i++)for(int j=0;j<n;j++) d[i*n+j]=std::round(std::hypot(x[i]-x[j],y[i]-y[j]));
	}
	solTSP construction(){ solTSP s; for(int i=0;i<n;i++) s.sol.push_back(i); std::shuffle(s.sol.begin(),s.sol.end(),rng()); s=evaluate(s); s.Nup=false; s.Ndown=false; return s; }
	solTSP neighbor(solTSP sol){
		std::uniform_int_distribution<int> u(0,n-1); int a,b; do{a=u(rng());b=u(rng());}while(a==b); if(a>b)std::swap(a,b);
		std::reverse(sol.sol.begin()+a, sol.sol.begin()+b+1); return sol; }
	solTSP evaluate(solTSP sol){ double c=0; for(int i=0;i<n;i++) c+=d[sol.sol[i]*n+sol.sol[(i+1)%n]]; sol.evalSol=c; ++tl.ev; noteMin(c); return sol; }
};

// ---------- timing decorator for the real problem
class TimedSD: public Problem<solSD>{
	SD* p;
public:
	TimedSD(SD* p_):p(p_){}
	solSD construction(){ solSD s=p->construction(); noteMin(s.evalSol); return s; }
	solSD neighbor(solSD s){ long long t=nowNs(); solSD r=p->neighbor(std::move(s)); tl.ns+=nowNs()-t; return r; }
	solSD evaluate(solSD s){ long long t=nowNs(); solSD r=p->evaluate(std::move(s)); tl.ns+=nowNs()-t; ++tl.ev; noteMin(r.evalSol); return r; }
};

int main(int argc, char** argv){
	std::string prob = argv[1];           // "tsp:<n>" or instance path
	int rep=atoi(argv[2]), mcl=atoi(argv[3]), ptl=atoi(argv[4]), th=atoi(argv[5]);
	int tEnd = argc>6? atoi(argv[6]):0; int ut = argc>7? atoi(argv[7]):0; int tup = argc>8? atoi(argv[8]):50;
	float tmin = argc>9? atof(argv[9]):0.3f, tmax = argc>10? atof(argv[10]):3.0f;
	double best; 
	long long t0; struct rusage r0,r1;
	if(prob.rfind("tsp:",0)==0){
		TSP* p = new TSP(atoi(prob.c_str()+4), 12345);
		PT<solTSP> algo(tmin,tmax,rep,mcl,ptl,3,ut,tup);
		if(tEnd) algo.setTimeEnd(tEnd);
		getrusage(RUSAGE_SELF,&r0); t0=nowNs();
		solTSP s = algo.start(th,p); best=s.evalSol;
	}else{
		SD* sd = new SD(prob,3,2); TimedSD* p = new TimedSD(sd);
		PT<solSD> algo(tmin,tmax,rep,mcl,ptl,3,ut,tup);
		if(tEnd) algo.setTimeEnd(tEnd);
		getrusage(RUSAGE_SELF,&r0); t0=nowNs();
		solSD s = algo.start(th,p); best=s.evalSol;
	}
	long long t1=nowNs(); getrusage(RUSAGE_SELF,&r1);
	auto tv=[](timeval a){return a.tv_sec+a.tv_usec*1e-6;};
	double wall=(t1-t0)*1e-9, cpu=tv(r1.ru_utime)+tv(r1.ru_stime)-tv(r0.ru_utime)-tv(r0.ru_stime);
	// worker thread_local counters are flushed at thread exit (after start() joins)
	printf("wall=%.3f cpu=%.2f evals=%lld evals/s=%.0f probCpu=%.2f best=%.0f minEvaluated=%.0f %s\n", wall, cpu, gEvals.load(), gEvals.load()/wall, gProbNs.load()*1e-9, best, gMinEval.load(), best>gMinEval.load()? "LOST-BEST":"");
	return 0;
}
