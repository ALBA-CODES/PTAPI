#ifndef __PTAPI_EXPERIMENT_H__
#define __PTAPI_EXPERIMENT_H__

// ---------------------------------------------------------------------------
// PTAPI v2 - management of experiments (optional header)
//
//   PTRunResult     everything the API knows about one execution
//   ptResult()      fills a PTRunResult from a PT object after start()
//   PTExperiment    runs R executions for each instance and writes one file
//                   per execution plus the summaries (CSV, text and LaTeX), the
//                   convergence curves and the times to the target
//
// Comparison.h (also optional) compares several experiments, with the usual
// statistical tests.
//
// The header knows nothing about the problem: the user gives a function that
// receives (instance, seed), runs the PT and returns a PTRunResult.
// ---------------------------------------------------------------------------

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <functional>
#include <utility>
#include <stdexcept>
#include <filesystem>
#include <limits>
#include <cmath>
#include <cstdio>
#include "PT.h"

// One execution.
struct PTRunResult{
	double value = std::numeric_limits<double>::quiet_NaN();   // best objective value
	unsigned long long seed = 0;
	double elapsedMs = 0.0;
	double timeToBestMs = 0.0;
	int rounds = 0;               // exchange rounds executed
	int bestRound = 0;            // round of the best solution
	long long evaluations = 0;
	int threads = 0;
	int stopReason = PT_STOP_PTL;
	long long roundTrips = 0;
	PTTiming timing{0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0};
	std::vector<PTImprovement> convergence;
	std::vector<double> temperatures;   // final ladder
	std::vector<double> initialTemperatures;
	std::vector<PTMoveStat> moves;      // proposed and accepted, one per registered move
	std::vector<double> swapRates;
	std::vector<double> moveRates;
	int improvements = 0;         // number of new best solutions (size of the convergence curve)
	std::string solution;         // optional: the best solution as text, written by the user
	// Optional: values of the user (name, value), e.g. measures of the problem. Each
	// name becomes a column of runs.csv and is kept in the file of the execution.
	std::vector<std::pair<std::string, std::string>> extra;
	void add(const std::string& name, const std::string& v){ extra.emplace_back(name, v); }
	void add(const std::string& name, double v){ std::ostringstream o; o << std::setprecision(10) << v; extra.emplace_back(name, o.str()); }
};

// Collects the results of the last start() of 'pt'.
template<typename S>
PTRunResult ptResult(PT<S>& pt, const S& best){
	PTRunResult r;
	r.value = best.evalSol;
	r.seed = pt.getSeed();
	r.elapsedMs = pt.getElapsedMs();
	r.timeToBestMs = pt.getTimeToBestMs();
	r.rounds = pt.getPTLEnd();
	r.bestRound = pt.getBestPTLSolIndex();
	r.evaluations = pt.getEvaluations();
	r.threads = pt.getThreadsUsed();
	r.stopReason = pt.getStopReason();
	r.roundTrips = pt.getRoundTrips();
	r.timing = pt.getTiming();
	r.convergence = pt.getConvergence();
	r.improvements = (int)r.convergence.size();
	r.temperatures = pt.getTemperatures();
	r.initialTemperatures = pt.getInitialTemperatures();
	r.moves = pt.getMoveStatistics();
	r.swapRates = pt.getSwapAcceptanceRates();
	r.moveRates = pt.getMoveAcceptanceRates();
	return r;
}

// Statistics of the executions of one instance (minimization).
struct PTInstanceSummary{
	std::string name;
	int runs = 0;
	double best = 0.0, mean = 0.0, median = 0.0, sd = 0.0, worst = 0.0;
	double meanSeconds = 0.0;             // mean time of one execution
	double meanSecondsToBest = 0.0;       // mean time to the best solution of the execution
	double meanEvaluations = 0.0;
	bool hasReference = false;            // best known value / lower bound given by the user
	double reference = 0.0;
	double gapBest = 0.0, gapMean = 0.0;  // 100*(x - reference)/|reference|
	int hits = 0;                         // executions that reached the reference
	double meanSecondsToReference = 0.0;  // mean time to reach it, over those executions
};

class PTExperiment{
	public:
		// (instance path or name, seed) -> result of one execution
		typedef std::function<PTRunResult(const std::string&, unsigned long long)> Runner;

	private:
		struct Instance{
			std::string name, path;
			bool hasReference = false;
			double reference = 0.0;
			std::vector<PTRunResult> runs;
		};
		std::string dir;
		std::vector<Instance> instances;
		std::vector<PTInstanceSummary> summary;
		int nRuns = 10;
		unsigned long long firstSeed = 1;
		bool resume = true;
		bool verbose = true;
		char sep = ',';
		bool saveConvergence = true;
		int shardIndex = 0, shardCount = 1;
		int curvePoints = 100;
		std::vector<std::string> columns{"instance", "run", "seed", "value", "seconds", "seconds_to_best", "rounds", "best_round",
		                                 "evaluations", "threads", "stop", "round_trips"};

		static std::vector<double> numbers(const std::string& text){
			std::vector<double> v;
			std::istringstream in(text);
			double x;
			while(in >> x) v.push_back(x);
			return v;
		}
		// CSV field: quoted if it contains the separator, quotes or line breaks
		std::string field(const std::string& v) const {
			if(v.find(sep) == std::string::npos && v.find('"') == std::string::npos && v.find('\n') == std::string::npos) return v;
			std::string q = "\"";
			for(char c : v){ if(c == '"') q += '"'; q += c; }
			return q + "\"";
		}
		// Value of one built-in column of runs.csv for execution k of an instance.
		static bool builtin(const std::string& col, const Instance& in, size_t k, const PTRunResult& r, std::string& out){
			const double sec = r.elapsedMs/1000.0;
			auto rate = [&r](int what){
				if(r.swapRates.empty()) return std::string();
				double mn = r.swapRates[0], mx = r.swapRates[0], sum = 0.0;
				for(double x : r.swapRates){ mn = std::min(mn, x); mx = std::max(mx, x); sum += x; }
				return num(what == 0 ? sum/(double)r.swapRates.size() : (what == 1 ? mn : (what == 2 ? mx : (double)r.swapRates.size() - sum)));
			};
			if(col == "instance") out = in.name;
			else if(col == "path") out = in.path;
			else if(col == "run") out = std::to_string(k + 1);
			else if(col == "seed") out = std::to_string(r.seed);
			else if(col == "value") out = num(r.value, 17);
			else if(col == "seconds") out = num(sec);
			else if(col == "seconds_to_best") out = num(r.timeToBestMs/1000.0);
			else if(col == "rounds") out = std::to_string(r.rounds);
			else if(col == "best_round") out = std::to_string(r.bestRound);
			else if(col == "evaluations") out = std::to_string(r.evaluations);
			else if(col == "evaluations_per_second") out = (sec > 0.0) ? num((double)r.evaluations/sec) : std::string();
			else if(col == "threads") out = std::to_string(r.threads);
			else if(col == "stop") out = PT<solution>::stopReasonName(r.stopReason);
			else if(col == "round_trips") out = std::to_string(r.roundTrips);
			else if(col == "improvements") out = std::to_string(r.improvements);
			else if(col == "replicas") out = std::to_string(r.temperatures.size());
			else if(col == "temp_first") out = r.temperatures.empty() ? std::string() : num(r.temperatures.front());
			else if(col == "temp_last") out = r.temperatures.empty() ? std::string() : num(r.temperatures.back());
			else if(col == "temperatures") out = list(r.temperatures);
			else if(col == "temperatures_initial") out = list(r.initialTemperatures);
			else if(col == "moves_proposed" || col == "moves_accepted" || col == "moves_rate"){
				std::vector<double> v;
				for(const PTMoveStat& m : r.moves) v.push_back(col == "moves_proposed" ? (double)m.proposed : (col == "moves_accepted" ? (double)m.accepted : (m.proposed ? (double)m.accepted/(double)m.proposed : 0.0)));
				out = list(v, col == "moves_rate" ? 6 : 17);
			}
			else if(col == "seconds_to_reference"){ const double t = in.hasReference ? secondsToReference(in, r) : -1.0; out = (t >= 0.0) ? num(t) : std::string(); }
			else if(col == "swap_rate_mean") out = rate(0);
			else if(col == "swap_rate_min") out = rate(1);
			else if(col == "swap_rate_max") out = rate(2);
			else if(col == "barrier") out = rate(3);
			else if(col == "swap_rates") out = list(r.swapRates);
			else if(col == "move_rates") out = list(r.moveRates);
			else if(col == "construction_seconds") out = num(r.timing.constructionMs/1000.0);
			else if(col == "chains_ms") out = num(r.timing.chainsMs);
			else if(col == "exchanges_ms") out = num(r.timing.exchangesMs);
			else if(col == "updates_ms") out = num(r.timing.updatesMs);
			else if(col == "waiting_ms") out = num(r.timing.waitingMs);
			else if(col == "reference") out = in.hasReference ? num(in.reference, 17) : std::string();
			else if(col == "gap_percent") out = in.hasReference ? num(100.0*(r.value - in.reference)/((in.reference != 0.0) ? std::fabs(in.reference) : 1.0)) : std::string();
			else if(col == "solution") out = r.solution;
			else return false;
			return true;
		}

		static std::string num(double v, int precision = 10){
			std::ostringstream o;
			o << std::setprecision(precision) << v;
			return o.str();
		}
		static std::string fixed(double v, int decimals){
			std::ostringstream o;
			o << std::fixed << std::setprecision(decimals) << v;
			return o.str();
		}
		template<typename T>
		static std::string list(const std::vector<T>& v, int precision = 10){
			std::ostringstream o;
			o << std::setprecision(precision);
			for(size_t i = 0; i < v.size(); ++i) o << (i ? " " : "") << v[i];
			return o.str();
		}
		static std::string baseName(const std::string& path){
			std::string s = std::filesystem::path(path).stem().string();
			return s.empty() ? path : s;
		}
		// name usable as a folder
		static std::string safe(std::string s){
			for(char& c : s) if(c == '/' || c == '\\' || c == ':' || c == ' ') c = '_';
			return s;
		}
		std::string runFile(const Instance& in, int k) const {
			return dir + "/" + safe(in.name) + "/run_" + std::to_string(k + 1) + ".txt";
		}
		static int reasonCode(const std::string& name){
			for(int r = PT_STOP_PTL; r <= PT_STOP_ERROR; ++r) if(name == PT<solution>::stopReasonName(r)) return r;
			return PT_STOP_PTL;
		}

		void writeRun(const Instance& in, int k, const PTRunResult& r) const {
			const std::string file = runFile(in, k), tmp = file + ".tmp";
			{
				std::ofstream o(tmp.c_str(), std::ios::trunc);
				if(!o) throw std::runtime_error("PTAPI: cannot write " + tmp);
				o << "instance: " << in.name << "\n";
				o << "path: " << in.path << "\n";
				o << "run: " << (k + 1) << "\n";
				o << "seed: " << r.seed << "\n";
				o << "value: " << num(r.value, 17) << "\n";
				o << "elapsed_ms: " << num(r.elapsedMs) << "\n";
				o << "time_to_best_ms: " << num(r.timeToBestMs) << "\n";
				o << "rounds: " << r.rounds << "\n";
				o << "best_round: " << r.bestRound << "\n";
				o << "evaluations: " << r.evaluations << "\n";
				o << "threads: " << r.threads << "\n";
				o << "stop: " << PT<solution>::stopReasonName(r.stopReason) << "\n";
				o << "round_trips: " << r.roundTrips << "\n";
				o << "improvements: " << r.improvements << "\n";
				o << "construction_ms: " << num(r.timing.constructionMs) << "\n";
				if(r.timing.chainsMs > 0.0){
					o << "chains_ms: " << num(r.timing.chainsMs) << "\n";
					o << "exchanges_ms: " << num(r.timing.exchangesMs) << "\n";
					o << "updates_ms: " << num(r.timing.updatesMs) << "\n";
					o << "waiting_ms: " << num(r.timing.waitingMs) << "\n";
				}
				o << "temperatures: " << list(r.temperatures) << "\n";
				o << "temperatures_initial: " << list(r.initialTemperatures) << "\n";
				if(!r.moves.empty()){
					std::vector<long long> a, b;
					for(const PTMoveStat& m : r.moves){ a.push_back(m.proposed); b.push_back(m.accepted); }
					o << "moves_proposed: " << list(a) << "\n";
					o << "moves_accepted: " << list(b) << "\n";
				}
				o << "swap_rates: " << list(r.swapRates) << "\n";
				o << "move_rates: " << list(r.moveRates) << "\n";
				for(const std::pair<std::string, std::string>& e : r.extra){
					std::string one = e.second;
					for(char& c : one) if(c == '\n' || c == '\r') c = ' ';
					o << "extra." << e.first << ": " << one << "\n";
				}
				if(!r.solution.empty()){
					std::string one = r.solution;
					for(char& c : one) if(c == '\n' || c == '\r') c = ' ';
					o << "solution: " << one << "\n";
				}
				if(saveConvergence){
					o << "convergence: " << r.convergence.size() << "   (ms round value)\n";
					for(const PTImprovement& p : r.convergence) o << num(p.timeMs) << " " << p.ptl << " " << num(p.value, 17) << "\n";
				}
				o << "complete: 1\n";
				o.flush();
				if(!o) throw std::runtime_error("PTAPI: error writing " + tmp);
			}
			std::filesystem::rename(tmp, file);
		}

		// Reads back an execution saved by writeRun. False if the file does not exist
		// or is not complete.
		bool readRun(const Instance& in, int k, PTRunResult& r) const {
			std::ifstream f(runFile(in, k).c_str());
			if(!f) return false;
			bool complete = false;
			std::string line;
			std::vector<double> proposed, accepted;
			while(std::getline(f, line)){
				const size_t p = line.find(": ");
				if(p == std::string::npos) continue;
				const std::string key = line.substr(0, p), val = line.substr(p + 2);
				try{
					if(key == "seed") r.seed = std::stoull(val);
					else if(key == "value") r.value = std::stod(val);
					else if(key == "elapsed_ms") r.elapsedMs = std::stod(val);
					else if(key == "time_to_best_ms") r.timeToBestMs = std::stod(val);
					else if(key == "rounds") r.rounds = std::stoi(val);
					else if(key == "best_round") r.bestRound = std::stoi(val);
					else if(key == "evaluations") r.evaluations = std::stoll(val);
					else if(key == "threads") r.threads = std::stoi(val);
					else if(key == "stop") r.stopReason = reasonCode(val);
					else if(key == "round_trips") r.roundTrips = std::stoll(val);
					else if(key == "improvements") r.improvements = std::stoi(val);
					else if(key == "construction_ms") r.timing.constructionMs = std::stod(val);
					else if(key == "chains_ms") r.timing.chainsMs = std::stod(val);
					else if(key == "exchanges_ms") r.timing.exchangesMs = std::stod(val);
					else if(key == "updates_ms") r.timing.updatesMs = std::stod(val);
					else if(key == "waiting_ms") r.timing.waitingMs = std::stod(val);
					else if(key == "temperatures") r.temperatures = numbers(val);
					else if(key == "temperatures_initial") r.initialTemperatures = numbers(val);
					else if(key == "moves_proposed") proposed = numbers(val);
					else if(key == "moves_accepted") accepted = numbers(val);
					else if(key == "convergence"){
						const int count = std::stoi(val);
						for(int i = 0; i < count && std::getline(f, line); ++i){
							const std::vector<double> x = numbers(line);
							if(x.size() != 3) return false;
							r.convergence.push_back(PTImprovement{x[0], (int)x[1], x[2]});
						}
					}
					else if(key == "swap_rates") r.swapRates = numbers(val);
					else if(key == "move_rates") r.moveRates = numbers(val);
					else if(key.rfind("extra.", 0) == 0) r.extra.emplace_back(key.substr(6), val);
					else if(key == "solution") r.solution = val;
					else if(key == "complete") complete = (val == "1");
				}catch(...){ return false; }
			}
			for(size_t i = 0; i < proposed.size() && i < accepted.size(); ++i) r.moves.push_back(PTMoveStat{(long long)proposed[i], (long long)accepted[i]});
			return complete;
		}

		static double tolerance(double reference){ return 1e-9*std::max(1.0, std::fabs(reference)); }
		// Seconds the execution took to reach the reference of the instance; < 0 if it
		// did not. From the convergence curve; without it, the time of the best solution.
		static double secondsToReference(const Instance& in, const PTRunResult& r){
			const double limit = in.reference + tolerance(in.reference);
			if(!std::isfinite(r.value) || !(r.value <= limit)) return -1.0;
			for(const PTImprovement& p : r.convergence) if(p.value <= limit) return p.timeMs/1000.0;
			return r.timeToBestMs/1000.0;
		}

		// Best value the execution had at position x of the curve (x in rounds or in
		// seconds). False before its first solution. The curve is the list of the new
		// best solutions in the order they were found; the replicas are not all in
		// the same round, so the rounds of the list are not always ascending.
		static bool valueAt(const PTRunResult& r, bool byTime, double x, double& v){
			bool any = false;
			for(const PTImprovement& p : r.convergence){
				if(byTime){ if(p.timeMs/1000.0 > x) break; }
				else if((double)p.ptl > x) continue;
				v = p.value; any = true;
			}
			return any;
		}
		// Convergence of an instance over its executions: for each point of a grid
		// (closer points at the beginning, where the curve changes fast), the best,
		// the median, the mean and the worst of the values the executions had by then.
		void writeCurve(const Instance& in, bool byTime) const {
			double last = 0.0, firstStep = std::numeric_limits<double>::infinity();
			bool any = false;
			for(const PTRunResult& r : in.runs){
				if(!std::isfinite(r.value) || r.convergence.empty()) continue;
				any = true;
				last = std::max(last, byTime ? r.elapsedMs/1000.0 : (double)r.rounds);
				for(const PTImprovement& p : r.convergence){
					const double x = byTime ? p.timeMs/1000.0 : (double)p.ptl;
					last = std::max(last, x);
					if(x > 0.0) firstStep = std::min(firstStep, x);
				}
			}
			const std::string file = dir + "/" + safe(in.name) + (byTime ? "/convergence_time.csv" : "/convergence.csv");
			if(!any || !(last > 0.0)){ std::filesystem::remove(file); return; }
			std::vector<double> grid;
			grid.push_back(0.0);
			const double lo = byTime ? std::min(firstStep, last) : 1.0;
			for(int j = 0; j < curvePoints; ++j){
				double x = (curvePoints > 1) ? lo*std::pow(last/lo, (double)j/(double)(curvePoints - 1)) : last;
				if(!byTime) x = std::floor(x + 0.5);
				if(x > grid.back()) grid.push_back(x);
			}
			if(grid.back() < last) grid.push_back(last);
			const std::string s(1, sep);
			std::ofstream o(file.c_str(), std::ios::trunc);
			o << (byTime ? "seconds" : "round") << s << "runs" << s << "best" << s << "median" << s << "mean" << s << "worst\n";
			for(double x : grid){
				std::vector<double> v;
				for(const PTRunResult& r : in.runs){
					double y;
					if(std::isfinite(r.value) && valueAt(r, byTime, x, y)) v.push_back(y);
				}
				if(v.empty()) continue;
				std::sort(v.begin(), v.end());
				double mean = 0.0;
				for(double y : v) mean += y/(double)v.size();
				const double median = (v.size() % 2) ? v[v.size()/2] : 0.5*(v[v.size()/2 - 1] + v[v.size()/2]);
				o << num(x) << s << v.size() << s << num(v.front(), 17) << s << num(median) << s << num(mean) << s << num(v.back(), 17) << "\n";
			}
		}
		// Time-to-target of an instance with a reference: the times of the executions
		// that reached it, in ascending order, with the empirical probability
		// (i - 0.5)/R of having reached it by then (R = executions of the instance).
		void writeTimeToTarget(const Instance& in) const {
			const std::string file = dir + "/" + safe(in.name) + "/time_to_target.csv";
			std::vector<std::pair<double, int>> t;
			int total = 0;
			if(in.hasReference) for(size_t k = 0; k < in.runs.size(); ++k){
				if(!std::isfinite(in.runs[k].value)) continue;
				++total;
				const double x = secondsToReference(in, in.runs[k]);
				if(x >= 0.0) t.emplace_back(x, (int)k + 1);
			}
			if(t.empty()){ std::filesystem::remove(file); return; }
			std::sort(t.begin(), t.end());
			const std::string s(1, sep);
			std::ofstream o(file.c_str(), std::ios::trunc);
			o << "seconds" << s << "probability" << s << "run\n";
			for(size_t i = 0; i < t.size(); ++i) o << num(t[i].first) << s << num(((double)i + 0.5)/(double)total) << s << t[i].second << "\n";
		}

		// Scripts that draw the files above: figures.gp (gnuplot) and figures.tex
		// (LaTeX with pgfplots, one figure per page). The API does not run them:
		//     cd <outputDir> && gnuplot figures.gp        ->  <instance>/convergence.pdf, time_to_target.pdf
		//     cd <outputDir> && lualatex figures.tex      ->  figures.pdf
		void writeFigureScripts() const {
			std::vector<std::string> names, folders;
			std::vector<bool> target;
			for(const Instance& in : instances){
				const std::string f = safe(in.name);
				if(!std::filesystem::exists(dir + "/" + f + "/convergence.csv")) continue;
				names.push_back(in.name); folders.push_back(f);
				target.push_back(std::filesystem::exists(dir + "/" + f + "/time_to_target.csv"));
			}
			if(names.empty()) return;
			const std::string gpSep = (sep == '\t') ? "\\t" : std::string(1, sep);
			const std::string texSep = (sep == ',') ? "comma" : (sep == ';' ? "semicolon" : (sep == '\t' ? "tab" : "space"));
			{
				std::ofstream o((dir + "/figures.gp").c_str(), std::ios::trunc);
				o << "# PTAPI - figures of the experiment. Run inside this folder:  gnuplot figures.gp\n"
				  << "set datafile separator \"" << gpSep << "\"\nset terminal pdfcairo size 12cm,8cm font \",10\"\nset grid\nset key top right\n";
				for(size_t i = 0; i < names.size(); ++i){
					o << "\nset output \"" << folders[i] << "/convergence.pdf\"\nset title \"" << names[i] << "\" noenhanced\n"
					  << "set xlabel \"round\"\nset ylabel \"value\"\nset logscale x\nset yrange [*:*]\n"
					  << "f = \"" << folders[i] << "/convergence.csv\"\n"
					  << "plot f using 1:3:6 skip 1 with filledcurves fc rgb \"#c8d8ee\" title \"best to worst\", \\\n"
					  << "     f using 1:4 skip 1 with lines lw 2 lc rgb \"#1f4e9c\" title \"median\", \\\n"
					  << "     f using 1:5 skip 1 with lines dt 2 lc rgb \"#444444\" title \"mean\"\n";
					if(target[i])
						o << "\nset output \"" << folders[i] << "/time_to_target.pdf\"\nset xlabel \"seconds to the target\"\nset ylabel \"probability\"\nunset logscale x\nset yrange [0:1]\n"
						  << "plot \"" << folders[i] << "/time_to_target.csv\" using 1:2 skip 1 with linespoints pt 7 ps 0.5 lc rgb \"#1f4e9c\" notitle\n";
				}
				o << "\nunset output\n";
			}
			{
				std::ofstream o((dir + "/figures.tex").c_str(), std::ios::trunc);
				o << "% PTAPI - figures of the experiment, one per page. Run inside this folder:  lualatex figures.tex\n"
				  << "% To use one figure in a paper, copy its tikzpicture (packages: pgfplots, library fillbetween).\n"
				  << "\\documentclass[tikz,border=4pt]{standalone}\n\\usepackage{pgfplots}\n\\usepgfplotslibrary{fillbetween}\n\\pgfplotsset{compat=1.16}\n\\begin{document}\n";
				for(size_t i = 0; i < names.size(); ++i){
					const std::string title = texEscape(names[i]);
					o << "\\begin{tikzpicture}\n\\begin{axis}[width=12cm, height=8cm, xmode=log, xlabel={round}, ylabel={value}, title={" << title << "}, grid=major,\n"
					  << "    legend pos=north east, unbounded coords=discard, filter discard warning=false, scaled y ticks=false, yticklabel style={/pgf/number format/fixed}]\n"
					  << "  \\addplot[name path=best, draw=none, forget plot] table[col sep=" << texSep << ", x=round, y=best]{" << folders[i] << "/convergence.csv};\n"
					  << "  \\addplot[name path=worst, draw=none, forget plot] table[col sep=" << texSep << ", x=round, y=worst]{" << folders[i] << "/convergence.csv};\n"
					  << "  \\addplot[blue!18] fill between[of=best and worst]; \\addlegendentry{best to worst}\n"
					  << "  \\addplot[blue!70!black, thick] table[col sep=" << texSep << ", x=round, y=median]{" << folders[i] << "/convergence.csv}; \\addlegendentry{median}\n"
					  << "  \\addplot[black!70, dashed] table[col sep=" << texSep << ", x=round, y=mean]{" << folders[i] << "/convergence.csv}; \\addlegendentry{mean}\n"
					  << "\\end{axis}\n\\end{tikzpicture}\n";
					if(target[i])
						o << "\\begin{tikzpicture}\n\\begin{axis}[width=12cm, height=8cm, xlabel={seconds to the target}, ylabel={probability}, title={" << title << "}, grid=major, ymin=0, ymax=1]\n"
						  << "  \\addplot[blue!70!black, mark=*, mark size=1.2pt] table[col sep=" << texSep << ", x=seconds, y=probability]{" << folders[i] << "/time_to_target.csv};\n"
						  << "\\end{axis}\n\\end{tikzpicture}\n";
				}
				o << "\\end{document}\n";
			}
		}

		static PTInstanceSummary summarize(const Instance& in){
			PTInstanceSummary s;
			s.name = in.name;
			std::vector<double> v;
			for(const PTRunResult& r : in.runs){
				if(!std::isfinite(r.value)) continue;
				v.push_back(r.value);
				s.meanSeconds += r.elapsedMs/1000.0;
				s.meanSecondsToBest += r.timeToBestMs/1000.0;
				s.meanEvaluations += (double)r.evaluations;
			}
			s.runs = (int)v.size();
			s.hasReference = in.hasReference;
			s.reference = in.reference;
			if(v.empty()) return s;
			const double n = (double)v.size();
			std::sort(v.begin(), v.end());
			s.best = v.front();
			s.worst = v.back();
			s.median = (v.size() % 2) ? v[v.size()/2] : 0.5*(v[v.size()/2 - 1] + v[v.size()/2]);
			for(double x : v) s.mean += x;
			s.mean /= n;
			double ss = 0.0;
			for(double x : v) ss += (x - s.mean)*(x - s.mean);
			s.sd = (v.size() > 1) ? std::sqrt(ss/(n - 1.0)) : 0.0;   // sample standard deviation
			s.meanSeconds /= n; s.meanSecondsToBest /= n; s.meanEvaluations /= n;
			if(in.hasReference){
				const double den = (in.reference != 0.0) ? std::fabs(in.reference) : 1.0;
				s.gapBest = 100.0*(s.best - in.reference)/den;
				s.gapMean = 100.0*(s.mean - in.reference)/den;
				for(const PTRunResult& r : in.runs){
					const double t = secondsToReference(in, r);
					if(t >= 0.0){ ++s.hits; s.meanSecondsToReference += t; }
				}
				if(s.hits) s.meanSecondsToReference /= (double)s.hits;
			}
			return s;
		}

		static std::string texEscape(const std::string& s){
			std::string o;
			for(char c : s){
				if(c == '_' || c == '%' || c == '&' || c == '#' || c == '$') o += '\\';
				o += c;
			}
			return o;
		}

		void writeReports(){
			if(saveConvergence) for(const Instance& in : instances){
				if(!std::filesystem::exists(dir + "/" + safe(in.name))) continue;
				writeCurve(in, false);
				writeCurve(in, true);
				writeTimeToTarget(in);
			}
			if(saveConvergence) writeFigureScripts();
			summary.clear();
			bool anyRef = false;
			for(const Instance& in : instances){ summary.push_back(summarize(in)); anyRef = anyRef || in.hasReference; }
			const std::string s(1, sep);

			{	// every execution, one per line: the chosen columns, then the values of the user
				std::vector<std::string> extras;
				for(const Instance& in : instances) for(const PTRunResult& r : in.runs) for(const std::pair<std::string, std::string>& e : r.extra)
					if(std::find(extras.begin(), extras.end(), e.first) == extras.end() && std::find(columns.begin(), columns.end(), e.first) == columns.end()) extras.push_back(e.first);
				std::ofstream o((dir + "/runs.csv").c_str(), std::ios::trunc);
				bool first = true;
				for(const std::string& c : columns){ o << (first ? "" : s) << field(c); first = false; }
				for(const std::string& c : extras){ o << (first ? "" : s) << field(c); first = false; }
				o << "\n";
				auto extraValue = [](const PTRunResult& r, const std::string& name){
					for(const std::pair<std::string, std::string>& e : r.extra) if(e.first == name) return e.second;
					return std::string();
				};
				for(const Instance& in : instances) for(size_t k = 0; k < in.runs.size(); ++k){
					const PTRunResult& r = in.runs[k];
					if(!std::isfinite(r.value)) continue;
					first = true;
					for(const std::string& c : columns){
						std::string v;
						if(!builtin(c, in, k, r, v)) v = extraValue(r, c);
						o << (first ? "" : s) << field(v); first = false;
					}
					for(const std::string& c : extras){ o << (first ? "" : s) << field(extraValue(r, c)); first = false; }
					o << "\n";
				}
			}
			{	// one line per instance
				std::ofstream o((dir + "/summary.csv").c_str(), std::ios::trunc);
				o << "instance" << s << "runs" << s << "best" << s << "mean" << s << "median" << s << "sd" << s << "worst" << s << "mean_seconds" << s
				  << "mean_seconds_to_best" << s << "mean_evaluations";
				if(anyRef) o << s << "reference" << s << "gap_best_percent" << s << "gap_mean_percent" << s << "hits" << s << "mean_seconds_to_reference";
				o << "\n";
				for(const PTInstanceSummary& m : summary){
					o << m.name << s << m.runs << s << num(m.best, 17) << s << num(m.mean) << s << num(m.median) << s << num(m.sd) << s << num(m.worst, 17) << s
					  << num(m.meanSeconds) << s << num(m.meanSecondsToBest) << s << num(m.meanEvaluations);
					if(anyRef){
						if(m.hasReference) o << s << num(m.reference, 17) << s << num(m.gapBest) << s << num(m.gapMean) << s << m.hits << s << (m.hits ? num(m.meanSecondsToReference) : std::string());
						else o << s << s << s << s << s;
					}
					o << "\n";
				}
			}
			{	// the same table, aligned, to read on the terminal
				std::ofstream o((dir + "/summary.txt").c_str(), std::ios::trunc);
				o << summaryTable();
			}
			{	// and as a LaTeX table (booktabs)
				std::ofstream o((dir + "/summary.tex").c_str(), std::ios::trunc);
				o << "% \\usepackage{booktabs}\n\\begin{tabular}{lrrrrrr" << (anyRef ? "rr" : "") << "}\n\\toprule\n";
				o << "Instance & Runs & Best & Mean & Median & SD & Time (s)" << (anyRef ? " & Gap best (\\%) & Gap mean (\\%)" : "") << " \\\\\n\\midrule\n";
				for(const PTInstanceSummary& m : summary){
					o << texEscape(m.name) << " & " << m.runs << " & " << fixed(m.best, 2) << " & " << fixed(m.mean, 2) << " & " << fixed(m.median, 2) << " & "
					  << fixed(m.sd, 2) << " & " << fixed(m.meanSeconds, 2);
					if(anyRef) o << " & " << (m.hasReference ? fixed(m.gapBest, 2) : std::string("--")) << " & " << (m.hasReference ? fixed(m.gapMean, 2) : std::string("--"));
					o << " \\\\\n";
				}
				o << "\\bottomrule\n\\end{tabular}\n";
			}
		}

	public:
		// The files are written under 'outputDir' (created if needed).
		explicit PTExperiment(const std::string& outputDir): dir(outputDir) {}

		// 'path' is what the runner receives (usually the instance file); the name used
		// in the reports is the file name without extension, unless 'name' is given.
		void addInstance(const std::string& path, const std::string& name = ""){
			Instance in;
			in.path = path;
			in.name = name.empty() ? baseName(path) : name;
			for(const Instance& other : instances) if(other.name == in.name) throw std::invalid_argument("PTAPI: two instances named " + in.name);
			instances.push_back(in);
		}
		// Best known value or lower bound of the last instance added: the summaries get
		// the gap to it and the number of executions that reached it.
		void setReference(double value){
			if(instances.empty()) throw std::logic_error("PTAPI: setReference before addInstance");
			instances.back().hasReference = true;
			instances.back().reference = value;
		}
		// One instance per line: "path" or "path reference". Lines starting with # are ignored.
		void addInstancesFromFile(const std::string& listFile){
			std::ifstream f(listFile.c_str());
			if(!f) throw std::runtime_error("PTAPI: cannot read " + listFile);
			std::string line;
			while(std::getline(f, line)){
				std::istringstream in(line);
				std::string path;
				double ref;
				if(!(in >> path) || path[0] == '#') continue;
				addInstance(path);
				if(in >> ref) setReference(ref);
			}
		}

		void setRuns(int runs){ if(runs < 1) throw std::invalid_argument("PTAPI: the number of runs must be >= 1"); nRuns = runs; }
		// Execution k (from 0) of every instance uses the seed firstSeed + k.
		void setFirstSeed(unsigned long long s){ firstSeed = s; }
		// true (default): executions whose file already exists are not repeated, so an
		// interrupted experiment continues where it stopped. false: everything runs again.
		void setResume(bool on){ resume = on; }
		void setVerbose(bool on){ verbose = on; }
		void setCsvSeparator(char c){ sep = c; }

		// Columns of runs.csv, in this order. Names of csvColumnNames(), or the name of a
		// value the runner puts in PTRunResult::extra. The values of 'extra' that are
		// not listed are appended as the last columns. Default: instance, run, seed,
		// value, seconds, seconds_to_best, rounds, best_round, evaluations, threads,
		// stop, round_trips.
		void setCsvColumns(const std::vector<std::string>& names){
			if(names.empty()) throw std::invalid_argument("PTAPI: setCsvColumns needs at least one column");
			columns = names;
		}
		// Adds columns to the current ones.
		void addCsvColumns(const std::vector<std::string>& names){ columns.insert(columns.end(), names.begin(), names.end()); }
		// Every column the API can fill by itself.
		static std::vector<std::string> csvColumnNames(){
			return {"instance", "path", "run", "seed", "value", "seconds", "seconds_to_best", "rounds", "best_round", "evaluations",
			        "evaluations_per_second", "threads", "stop", "round_trips", "improvements", "replicas", "temp_first", "temp_last",
			        "temperatures", "temperatures_initial", "moves_proposed", "moves_accepted", "moves_rate", "seconds_to_reference", "swap_rate_mean", "swap_rate_min", "swap_rate_max", "barrier", "swap_rates", "move_rates",
			        "construction_seconds", "chains_ms", "exchanges_ms", "updates_ms", "waiting_ms", "reference", "gap_percent", "solution"};
		}
		// false: the individual files do not list the convergence curve, and the
		// convergence and time-to-target files of the instances are not written.
		void setSaveConvergence(bool on){ saveConvergence = on; }
		// Number of points of the convergence curves of the instances (default 100).
		void setCurvePoints(int n){ if(n < 2) throw std::invalid_argument("PTAPI: setCurvePoints needs at least 2 points"); curvePoints = n; }
		// Divides the experiment among 'count' processes or machines that write in the
		// same folder (or in folders joined afterwards): this one only runs the
		// executions of its part, index in 0..count-1. The reports of each part are
		// partial; report() in the folder with every file writes the complete ones.
		void setShard(int index, int count){
			if(count < 1 || index < 0 || index >= count) throw std::invalid_argument("PTAPI: setShard needs 0 <= index < count");
			shardIndex = index; shardCount = count;
		}

		// Runs everything and writes the reports. They are rewritten after each
		// instance, so partial results can be read while the experiment runs.
		void run(Runner runner){
			if(instances.empty()) throw std::logic_error("PTAPI: experiment without instances");
			std::filesystem::create_directories(dir);
			long long job = 0;
			for(Instance& in : instances){
				if(runner) std::filesystem::create_directories(dir + "/" + safe(in.name));
				in.runs.assign(nRuns, PTRunResult());
				for(int k = 0; k < nRuns; ++k, ++job){
					PTRunResult r;
					const bool had = (resume || !runner) && readRun(in, k, r);
					if(!had){
						if(!runner || (job % shardCount) != shardIndex) continue;      // not of this part: stays missing
						r = runner(in.path, firstSeed + (unsigned long long)k);
						writeRun(in, k, r);
					}
					in.runs[k] = r;
					if(verbose){
						std::printf("%-24s run %3d/%d  value %-14s %8.2f s  stop: %s%s\n", in.name.c_str(), k + 1, nRuns, num(r.value).c_str(),
							r.elapsedMs/1000.0, PT<solution>::stopReasonName(r.stopReason), had ? "  (already done)" : "");
						std::fflush(stdout);
					}
				}
				if(runner) writeReports();
			}
			if(!runner) writeReports();
			if(verbose) std::cout << "\n" << summaryTable() << "\nfiles in " << dir << ": runs.csv, summary.csv, summary.txt, summary.tex and one folder per instance\n";
		}
		// Runs nothing: reads the executions that are in the folder and writes the
		// reports (after the parts of setShard, or to rebuild them with other columns).
		void report(){ run(Runner()); }

		const std::vector<PTInstanceSummary>& getSummary() const { return summary; }
		// Executions of instance i, in the order of the seeds.
		const std::vector<PTRunResult>& getRuns(size_t i) const { return instances.at(i).runs; }
		const std::string& getOutputDir() const { return dir; }

		std::string summaryTable() const {
			bool anyRef = false;
			size_t w = 8;
			for(const PTInstanceSummary& m : summary){ anyRef = anyRef || m.hasReference; w = std::max(w, m.name.size()); }
			std::ostringstream o;
			o << std::left << std::setw((int)w) << "instance" << std::right << std::setw(6) << "runs" << std::setw(14) << "best" << std::setw(14) << "mean"
			  << std::setw(14) << "median" << std::setw(12) << "sd" << std::setw(14) << "worst" << std::setw(10) << "time(s)" << std::setw(12) << "to best(s)";
			if(anyRef) o << std::setw(11) << "gap best%" << std::setw(11) << "gap mean%" << std::setw(6) << "hits";
			o << "\n";
			for(const PTInstanceSummary& m : summary){
				o << std::left << std::setw((int)w) << m.name << std::right << std::setw(6) << m.runs << std::setw(14) << fixed(m.best, 2) << std::setw(14) << fixed(m.mean, 2)
				  << std::setw(14) << fixed(m.median, 2) << std::setw(12) << fixed(m.sd, 2) << std::setw(14) << fixed(m.worst, 2) << std::setw(10) << fixed(m.meanSeconds, 2)
				  << std::setw(12) << fixed(m.meanSecondsToBest, 2);
				if(anyRef){
					if(m.hasReference) o << std::setw(11) << fixed(m.gapBest, 2) << std::setw(11) << fixed(m.gapMean, 2) << std::setw(6) << m.hits;
					else o << std::setw(11) << "-" << std::setw(11) << "-" << std::setw(6) << "-";
				}
				o << "\n";
			}
			return o.str();
		}
};

#endif
