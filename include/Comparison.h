#ifndef __PTAPI_COMPARISON_H__
#define __PTAPI_COMPARISON_H__

// ---------------------------------------------------------------------------
// PTAPI v2 - comparison of experiments (optional header)
//
//   ptstat::         the statistical tests (usable on their own)
//   PTComparison     compares the experiments (PTExperiment) of several
//                    configurations over the same instances
//
// The first configuration added is the CONTROL: every other configuration is
// compared with it.
//   - In each instance, with the Mann-Whitney U test on the values of the
//     executions (are the two samples of the instance different?).
//   - Over all the instances, with the Wilcoxon signed-rank test on the means
//     of the instances (paired by instance); with more than one comparison the
//     p-values are also given corrected by the method of Holm.
//   - With three or more configurations, the mean rank of each one and the
//     Friedman test (is any configuration different from the others?).
// All the tests are two-sided and use ranks, so they do not suppose normal
// values. Minimization: smaller is better.
// ---------------------------------------------------------------------------

#include <map>
#include "Experiment.h"

namespace ptstat{

// Ranks from 1, ties with the mean of their ranks. 'tieTerm' receives the sum of
// t^3 - t over the groups of t equal values.
inline std::vector<double> ranks(const std::vector<double>& v, double* tieTerm = nullptr){
	const size_t n = v.size();
	std::vector<size_t> idx(n);
	for(size_t i = 0; i < n; ++i) idx[i] = i;
	std::sort(idx.begin(), idx.end(), [&v](size_t a, size_t b){ return v[a] < v[b]; });
	std::vector<double> r(n, 0.0);
	double ties = 0.0;
	for(size_t i = 0; i < n; ){
		size_t j = i;
		while(j + 1 < n && v[idx[j+1]] == v[idx[i]]) ++j;
		const double t = (double)(j - i + 1);
		for(size_t k = i; k <= j; ++k) r[idx[k]] = 0.5*((double)i + (double)j) + 1.0;
		ties += t*t*t - t;
		i = j + 1;
	}
	if(tieTerm) *tieTerm = ties;
	return r;
}

// P(|Z| >= z) of the standard normal distribution
inline double normalTwoSided(double z){ return std::erfc(std::fabs(z)/std::sqrt(2.0)); }

// P(X >= x) of the chi-square distribution with 'df' degrees of freedom
inline double chiSquareUpper(double x, double df){
	if(!(x > 0.0) || !(df > 0.0)) return 1.0;
	const double a = 0.5*df, y = 0.5*x, lg = std::lgamma(a);
	if(y < a + 1.0){                       // series of the lower incomplete gamma function
		double term = 1.0/a, sum = term;
		for(int n = 1; n < 10000; ++n){
			term *= y/(a + (double)n);
			sum += term;
			if(std::fabs(term) < std::fabs(sum)*1e-15) break;
		}
		return std::min(1.0, std::max(0.0, 1.0 - sum*std::exp(-y + a*std::log(y) - lg)));
	}
	const double tiny = 1e-300;            // continued fraction of the upper one (Lentz)
	double b = y + 1.0 - a, c = 1.0/tiny, d = 1.0/b, h = d;
	for(int n = 1; n < 10000; ++n){
		const double an = -(double)n*((double)n - a);
		b += 2.0;
		d = an*d + b; if(std::fabs(d) < tiny) d = tiny;
		c = b + an/c; if(std::fabs(c) < tiny) c = tiny;
		d = 1.0/d;
		const double del = d*c;
		h *= del;
		if(std::fabs(del - 1.0) < 1e-15) break;
	}
	return std::min(1.0, std::max(0.0, h*std::exp(-y + a*std::log(y) - lg)));
}

// Mann-Whitney U test (Wilcoxon rank-sum) of two independent samples: two-sided
// p-value of "the values of one sample tend to be smaller than those of the other".
// Exact (distribution of the rank sum over all the ways of dividing the values,
// ties included) up to 60 values in total; normal approximation with correction
// for ties and continuity above that. NaN if a sample is empty.
inline double mannWhitney(const std::vector<double>& a, const std::vector<double>& b){
	const size_t n1 = a.size(), n2 = b.size(), n = n1 + n2;
	if(n1 == 0 || n2 == 0) return std::numeric_limits<double>::quiet_NaN();
	std::vector<double> all(a);
	all.insert(all.end(), b.begin(), b.end());
	double ties = 0.0;
	const std::vector<double> r = ranks(all, &ties);
	double r1 = 0.0;
	for(size_t i = 0; i < n1; ++i) r1 += r[i];
	const double meanR1 = 0.5*(double)n1*((double)n + 1.0);
	if(n <= 60){
		// ways[c][s]: subsets of c values whose doubled ranks add up to s
		const int maxSum = (int)(n*(n + 1));
		std::vector<std::vector<double>> ways(n1 + 1, std::vector<double>(maxSum + 1, 0.0));
		ways[0][0] = 1.0;
		for(size_t i = 0; i < n; ++i){
			const int w = (int)std::lround(2.0*r[i]);
			for(size_t c = std::min(i + 1, n1); c >= 1; --c)
				for(int s = maxSum; s >= w; --s) ways[c][s] += ways[c-1][s-w];
		}
		const double dev = std::fabs(2.0*r1 - 2.0*meanR1) - 1e-9;
		double total = 0.0, tail = 0.0;
		for(int s = 0; s <= maxSum; ++s){
			total += ways[n1][s];
			if(std::fabs((double)s - 2.0*meanR1) >= dev) tail += ways[n1][s];
		}
		return std::min(1.0, tail/total);
	}
	const double N = (double)n;
	const double var = (double)n1*(double)n2/12.0*((N + 1.0) - ties/(N*(N - 1.0)));
	if(!(var > 0.0)) return 1.0;
	const double z = std::max(0.0, std::fabs(r1 - meanR1) - 0.5)/std::sqrt(var);
	return normalTwoSided(z);
}

// Wilcoxon signed-rank test of paired samples (a[i] with b[i]): two-sided p-value
// of "the differences a[i] - b[i] are centred on zero". Pairs with difference
// zero are discarded. Exact up to 50 pairs (ties included); normal approximation
// with correction for ties and continuity above that. 1 if every pair is equal.
inline double wilcoxonSignedRank(const std::vector<double>& a, const std::vector<double>& b){
	if(a.size() != b.size()) throw std::invalid_argument("PTAPI: wilcoxonSignedRank needs samples of the same size");
	std::vector<double> absd;
	std::vector<bool> positive;
	for(size_t i = 0; i < a.size(); ++i){
		const double d = a[i] - b[i];
		if(d == 0.0 || !std::isfinite(d)) continue;
		absd.push_back(std::fabs(d));
		positive.push_back(d > 0.0);
	}
	const size_t n = absd.size();
	if(n == 0) return 1.0;
	double ties = 0.0;
	const std::vector<double> r = ranks(absd, &ties);
	double wPlus = 0.0;
	for(size_t i = 0; i < n; ++i) if(positive[i]) wPlus += r[i];
	const double mean = 0.25*(double)n*((double)n + 1.0);
	if(n <= 50){
		const int maxSum = (int)(n*(n + 1));
		std::vector<double> ways(maxSum + 1, 0.0);
		ways[0] = 1.0;
		for(size_t i = 0; i < n; ++i){
			const int w = (int)std::lround(2.0*r[i]);
			for(int s = maxSum; s >= w; --s) ways[s] += ways[s-w];
		}
		const double dev = std::fabs(2.0*wPlus - 2.0*mean) - 1e-9;
		double total = 0.0, tail = 0.0;
		for(int s = 0; s <= maxSum; ++s){
			total += ways[s];
			if(std::fabs((double)s - 2.0*mean) >= dev) tail += ways[s];
		}
		return std::min(1.0, tail/total);
	}
	const double N = (double)n;
	const double var = N*(N + 1.0)*(2.0*N + 1.0)/24.0 - ties/48.0;
	if(!(var > 0.0)) return 1.0;
	const double z = std::max(0.0, std::fabs(wPlus - mean) - 0.5)/std::sqrt(var);
	return normalTwoSided(z);
}

// Friedman test: values[i][j] is the result of configuration j in instance i.
// The configurations are ranked inside each instance (1 = smallest) and the
// p-value is that of "every configuration has the same mean rank" (chi-square
// approximation, corrected for ties). 'meanRanks' receives the mean rank of each
// configuration. NaN with fewer than 2 instances or 2 configurations.
inline double friedman(const std::vector<std::vector<double>>& values, std::vector<double>* meanRanks = nullptr){
	const size_t n = values.size(), k = n ? values[0].size() : 0;
	std::vector<double> sum(k, 0.0);
	double ties = 0.0;
	for(const std::vector<double>& row : values){
		if(row.size() != k) throw std::invalid_argument("PTAPI: friedman needs the same number of configurations in every instance");
		double t = 0.0;
		const std::vector<double> r = ranks(row, &t);
		ties += t;
		for(size_t j = 0; j < k; ++j) sum[j] += r[j];
	}
	if(meanRanks){
		meanRanks->assign(k, 0.0);
		for(size_t j = 0; j < k; ++j) (*meanRanks)[j] = n ? sum[j]/(double)n : 0.0;
	}
	if(n < 2 || k < 2) return std::numeric_limits<double>::quiet_NaN();
	const double N = (double)n, K = (double)k;
	double ss = 0.0;
	for(size_t j = 0; j < k; ++j) ss += (sum[j]/N)*(sum[j]/N);
	double chi = 12.0*N/(K*(K + 1.0))*(ss - K*(K + 1.0)*(K + 1.0)/4.0);
	const double correction = 1.0 - ties/(N*K*(K*K - 1.0));
	if(!(correction > 0.0)) return 1.0;                 // every instance is a complete tie
	chi /= correction;
	return chiSquareUpper(chi, K - 1.0);
}

// Method of Holm for several comparisons: the adjusted p-values, in the given order.
inline std::vector<double> holm(const std::vector<double>& p){
	const size_t m = p.size();
	std::vector<size_t> idx(m);
	for(size_t i = 0; i < m; ++i) idx[i] = i;
	std::sort(idx.begin(), idx.end(), [&p](size_t a, size_t b){ return p[a] < p[b]; });
	std::vector<double> adj(m, 1.0);
	double running = 0.0;
	for(size_t i = 0; i < m; ++i){
		running = std::max(running, std::min(1.0, (double)(m - i)*p[idx[i]]));
		adj[idx[i]] = running;
	}
	return adj;
}

}   // namespace ptstat

// One configuration in one instance.
struct PTComparisonCell{
	int runs = 0;
	double best = 0.0, mean = 0.0, median = 0.0, sd = 0.0, meanSeconds = 0.0;
	double pValue = std::numeric_limits<double>::quiet_NaN();   // Mann-Whitney against the control
	int sign = 0;                 // +1 better than the control, -1 worse, 0 no significant difference (or the control)
};
// One configuration over all the instances.
struct PTComparisonTotal{
	std::string name;
	double meanRank = 0.0;        // 1 = best mean in every instance
	double meanGap = 0.0;         // mean over the instances of the gap of the mean value (%)
	int bestMean = 0;             // instances in which it has the best mean (ties count)
	int wins = 0, ties = 0, losses = 0;   // against the control, by the test of each instance
	double pValue = std::numeric_limits<double>::quiet_NaN();   // Wilcoxon signed-rank against the control
	double pHolm = std::numeric_limits<double>::quiet_NaN();
};

class PTComparison{
	private:
		struct Config{
			std::string name;
			std::map<std::string, std::vector<double>> values;   // instance -> values of the executions
			std::map<std::string, PTInstanceSummary> summary;
			std::vector<std::string> order;
		};
		std::vector<Config> configs;
		double alpha = 0.05;
		char sep = ',';
		// results of compute()
		std::vector<std::string> instances;
		std::vector<std::vector<PTComparisonCell>> cells;   // [instance][configuration]
		std::vector<PTComparisonTotal> totals;
		double friedmanP = std::numeric_limits<double>::quiet_NaN();
		bool gapToReference = false;

		static std::string fixed(double v, int decimals){ std::ostringstream o; o << std::fixed << std::setprecision(decimals) << v; return o.str(); }
		static std::string num(double v, int precision = 10){ std::ostringstream o; o << std::setprecision(precision) << v; return o.str(); }
		static std::string pText(double p){
			if(!std::isfinite(p)) return "-";
			if(p < 0.0001) return "<0.0001";
			return fixed(p, 4);
		}
		static const char* signText(int s){ return s > 0 ? "+" : (s < 0 ? "-" : "="); }
		static std::string texEscape(const std::string& s){
			std::string o;
			for(char c : s){ if(c == '_' || c == '%' || c == '&' || c == '#' || c == '$') o += '\\'; o += c; }
			return o;
		}

	public:
		// The executions of 'experiment' (after run() or report()) under the name
		// 'name'. The first configuration added is the control.
		void add(const std::string& name, const PTExperiment& experiment){
			for(const Config& c : configs) if(c.name == name) throw std::invalid_argument("PTAPI: two configurations named " + name);
			Config c;
			c.name = name;
			const std::vector<PTInstanceSummary>& sum = experiment.getSummary();
			for(size_t i = 0; i < sum.size(); ++i){
				std::vector<double> v;
				for(const PTRunResult& r : experiment.getRuns(i)) if(std::isfinite(r.value)) v.push_back(r.value);
				if(v.empty()) continue;
				c.values[sum[i].name] = v;
				c.summary[sum[i].name] = sum[i];
				c.order.push_back(sum[i].name);
			}
			if(c.order.empty()) throw std::invalid_argument("PTAPI: the experiment of " + name + " has no executions (call run() or report() first)");
			configs.push_back(c);
			instances.clear();
		}
		// Significance level of the marks +, - and = (default 0.05).
		void setAlpha(double a){ if(!(a > 0.0) || !(a < 1.0)) throw std::invalid_argument("PTAPI: setAlpha needs 0 < alpha < 1"); alpha = a; instances.clear(); }
		void setCsvSeparator(char c){ sep = c; }

		// Computes everything (called by the functions below when needed). Only the
		// instances that every configuration has are compared.
		void compute(){
			if(configs.size() < 2) throw std::logic_error("PTAPI: a comparison needs at least two configurations");
			instances.clear(); cells.clear(); totals.clear();
			for(const std::string& name : configs[0].order){
				bool everywhere = true;
				for(const Config& c : configs) everywhere = everywhere && c.values.count(name);
				if(everywhere) instances.push_back(name);
			}
			if(instances.empty()) throw std::logic_error("PTAPI: the configurations have no instance in common");
			const size_t k = configs.size(), n = instances.size();

			gapToReference = true;
			for(const std::string& name : instances) gapToReference = gapToReference && configs[0].summary.at(name).hasReference;

			totals.assign(k, PTComparisonTotal());
			for(size_t j = 0; j < k; ++j) totals[j].name = configs[j].name;
			std::vector<std::vector<double>> means(n, std::vector<double>(k, 0.0));
			for(size_t i = 0; i < n; ++i){
				std::vector<PTComparisonCell> row(k);
				double bestMean = std::numeric_limits<double>::infinity();
				for(size_t j = 0; j < k; ++j){
					const PTInstanceSummary& s = configs[j].summary.at(instances[i]);
					row[j].runs = s.runs; row[j].best = s.best; row[j].mean = s.mean; row[j].median = s.median; row[j].sd = s.sd; row[j].meanSeconds = s.meanSeconds;
					means[i][j] = s.mean;
					bestMean = std::min(bestMean, s.mean);
					if(j > 0){
						row[j].pValue = ptstat::mannWhitney(configs[j].values.at(instances[i]), configs[0].values.at(instances[i]));
						if(row[j].pValue < alpha && s.mean != row[0].mean) row[j].sign = (s.mean < row[0].mean) ? 1 : -1;
						if(row[j].sign > 0) ++totals[j].wins; else if(row[j].sign < 0) ++totals[j].losses; else ++totals[j].ties;
					}
				}
				const PTInstanceSummary& s0 = configs[0].summary.at(instances[i]);
				const double base = gapToReference ? s0.reference : bestMean;
				const double den = (base != 0.0) ? std::fabs(base) : 1.0;
				for(size_t j = 0; j < k; ++j){
					totals[j].meanGap += 100.0*(row[j].mean - base)/den/(double)n;
					if(row[j].mean == bestMean) ++totals[j].bestMean;
				}
				cells.push_back(row);
			}
			std::vector<double> meanRanks;
			friedmanP = ptstat::friedman(means, &meanRanks);
			if(k < 3) friedmanP = std::numeric_limits<double>::quiet_NaN();
			std::vector<double> p;
			for(size_t j = 0; j < k; ++j){
				totals[j].meanRank = meanRanks[j];
				if(j == 0) continue;
				std::vector<double> a, b;
				for(size_t i = 0; i < n; ++i){ a.push_back(means[i][j]); b.push_back(means[i][0]); }
				totals[j].pValue = ptstat::wilcoxonSignedRank(a, b);
				p.push_back(totals[j].pValue);
			}
			const std::vector<double> adj = ptstat::holm(p);
			for(size_t j = 1; j < k; ++j) totals[j].pHolm = adj[j-1];
		}

		const std::vector<std::string>& getInstances(){ if(instances.empty()) compute(); return instances; }
		// [instance][configuration], in the order they were added
		const std::vector<std::vector<PTComparisonCell>>& getCells(){ if(instances.empty()) compute(); return cells; }
		const std::vector<PTComparisonTotal>& getTotals(){ if(instances.empty()) compute(); return totals; }
		// p-value of the Friedman test (NaN with fewer than three configurations or two instances)
		double getFriedmanP(){ if(instances.empty()) compute(); return friedmanP; }

		// The comparison as text.
		std::string table(){
			if(instances.empty()) compute();
			const size_t k = configs.size();
			size_t w = 8, wc = 12;
			for(const std::string& s : instances) w = std::max(w, s.size());
			for(const Config& c : configs) wc = std::max(wc, c.name.size() + 2);
			std::ostringstream o;
			o << "Mean value of each configuration in each instance. Control: " << configs[0].name << ".\n"
			  << "p: Mann-Whitney U test against the control on the executions of the instance (two-sided).\n"
			  << "+ better than the control, - worse, = no significant difference (alpha = " << alpha << ").\n\n";
			o << std::left << std::setw((int)w) << "instance" << std::right;
			for(size_t j = 0; j < k; ++j){ o << std::setw((int)wc) << configs[j].name; if(j) o << std::setw(10) << "p" << "  "; }
			o << "\n";
			for(size_t i = 0; i < instances.size(); ++i){
				o << std::left << std::setw((int)w) << instances[i] << std::right;
				for(size_t j = 0; j < k; ++j){
					o << std::setw((int)wc) << fixed(cells[i][j].mean, 2);
					if(j) o << std::setw(10) << pText(cells[i][j].pValue) << " " << signText(cells[i][j].sign);
				}
				o << "\n";
			}
			o << "\nOver the " << instances.size() << " instances. Gap: mean value against " << (gapToReference ? "the reference of the instance" : "the best mean of the instance")
			  << " (%).\np: Wilcoxon signed-rank test against the control on the means of the instances (two-sided)";
			if(k > 2) o << "; Holm: the same, corrected for the " << (k - 1) << " comparisons";
			o << ".\n\n";
			o << std::left << std::setw((int)wc) << "configuration" << std::right << std::setw(11) << "mean rank" << std::setw(10) << "gap %" << std::setw(11) << "best mean"
			  << std::setw(7) << "+" << std::setw(5) << "=" << std::setw(5) << "-" << std::setw(10) << "p";
			if(k > 2) o << std::setw(10) << "Holm";
			o << "\n";
			for(size_t j = 0; j < k; ++j){
				const PTComparisonTotal& t = totals[j];
				o << std::left << std::setw((int)wc) << t.name << std::right << std::setw(11) << fixed(t.meanRank, 2) << std::setw(10) << fixed(t.meanGap, 2) << std::setw(11) << t.bestMean;
				if(j) o << std::setw(7) << t.wins << std::setw(5) << t.ties << std::setw(5) << t.losses << std::setw(10) << pText(t.pValue);
				else o << std::setw(7) << "" << std::setw(5) << "" << std::setw(5) << "" << std::setw(10) << "control";
				if(k > 2) o << std::setw(10) << (j ? pText(t.pHolm) : std::string(""));
				o << "\n";
			}
			if(std::isfinite(friedmanP)) o << "\nFriedman test (the configurations have the same mean rank): p = " << pText(friedmanP) << "\n";
			if(instances.size() < 6) o << "\nWith " << instances.size() << " instances the tests over the instances have little power: fewer than 6 pairs cannot give p < 0.05.\n";
			return o.str();
		}

		// Writes comparison.txt, comparison.csv (one line per instance and configuration),
		// comparison_total.csv (one line per configuration) and comparison.tex in 'dir'.
		void write(const std::string& dir){
			if(instances.empty()) compute();
			std::filesystem::create_directories(dir);
			const std::string s(1, sep);
			const size_t k = configs.size();
			{ std::ofstream o((dir + "/comparison.txt").c_str(), std::ios::trunc); o << table(); }
			{
				std::ofstream o((dir + "/comparison.csv").c_str(), std::ios::trunc);
				o << "instance" << s << "configuration" << s << "runs" << s << "best" << s << "mean" << s << "median" << s << "sd" << s << "mean_seconds" << s << "p_against_control" << s << "sign\n";
				for(size_t i = 0; i < instances.size(); ++i) for(size_t j = 0; j < k; ++j){
					const PTComparisonCell& c = cells[i][j];
					o << instances[i] << s << configs[j].name << s << c.runs << s << num(c.best, 17) << s << num(c.mean) << s << num(c.median) << s << num(c.sd) << s << num(c.meanSeconds) << s
					  << (j ? num(c.pValue) : std::string()) << s << (j ? signText(c.sign) : "") << "\n";
				}
			}
			{
				std::ofstream o((dir + "/comparison_total.csv").c_str(), std::ios::trunc);
				o << "configuration" << s << "instances" << s << "mean_rank" << s << "mean_gap_percent" << s << "best_mean" << s << "wins" << s << "ties" << s << "losses" << s << "p_wilcoxon" << s << "p_holm" << s << "p_friedman\n";
				for(size_t j = 0; j < k; ++j){
					const PTComparisonTotal& t = totals[j];
					o << t.name << s << instances.size() << s << num(t.meanRank) << s << num(t.meanGap) << s << t.bestMean << s;
					if(j) o << t.wins << s << t.ties << s << t.losses << s << num(t.pValue) << s << num(t.pHolm); else o << s << s << s << s;
					o << s << (std::isfinite(friedmanP) ? num(friedmanP) : std::string()) << "\n";
				}
			}
			{	// LaTeX: the means, the best of each instance in bold, the mark of the test
				std::ofstream o((dir + "/comparison.tex").c_str(), std::ios::trunc);
				o << "% \\usepackage{booktabs}\n% Mean of each configuration; bold: best mean of the instance.\n"
				  << "% +, -, =: better, worse or not different from " << texEscape(configs[0].name) << " (Mann-Whitney U, alpha = " << alpha << ").\n";
				o << "\\begin{tabular}{l" << std::string(k, 'r') << "}\n\\toprule\nInstance";
				for(const Config& c : configs) o << " & " << texEscape(c.name);
				o << " \\\\\n\\midrule\n";
				for(size_t i = 0; i < instances.size(); ++i){
					double bestMean = cells[i][0].mean;
					for(size_t j = 1; j < k; ++j) bestMean = std::min(bestMean, cells[i][j].mean);
					o << texEscape(instances[i]);
					for(size_t j = 0; j < k; ++j){
						const std::string v = fixed(cells[i][j].mean, 2);
						o << " & " << (cells[i][j].mean == bestMean ? "\\textbf{" + v + "}" : v);
						if(j) o << " $" << signText(cells[i][j].sign) << "$";
					}
					o << " \\\\\n";
				}
				o << "\\midrule\nMean rank";
				for(size_t j = 0; j < k; ++j) o << " & " << fixed(totals[j].meanRank, 2);
				o << " \\\\\nGap (\\%)";
				for(size_t j = 0; j < k; ++j) o << " & " << fixed(totals[j].meanGap, 2);
				o << " \\\\\n$+$ / $=$ / $-$";
				for(size_t j = 0; j < k; ++j){ if(j) o << " & " << totals[j].wins << " / " << totals[j].ties << " / " << totals[j].losses; else o << " & --"; }
				o << " \\\\\n$p$ (Wilcoxon)";
				for(size_t j = 0; j < k; ++j) { std::string t = j ? pText(totals[j].pValue) : std::string("--"); if(t[0] == '<') t = "$<$" + t.substr(1); o << " & " << t; }
				o << " \\\\\n\\bottomrule\n\\end{tabular}\n";
			}
		}
};

#endif
