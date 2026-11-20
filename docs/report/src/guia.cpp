// Programa-guia da PTAPI v2: particao de numeros.
// Dividir N numeros em dois grupos com somas tao proximas quanto possivel.
//   g++ guia.cpp -std=c++17 -O3 -pthread -o guia
#include <cstdio>
#include <cmath>
#include "../../../include/PT.h"

// Passo 1: a solucao. Herda de 'solution', que traz o campo evalSol (custo).
struct solPart: public solution{
	std::vector<char> lado;   // lado[i] = 0 ou 1: grupo do numero i
	double diff = 0;          // soma do grupo 0 menos soma do grupo 1
};

// Passo 2: o problema. Tudo o que depende do problema fica nesta classe.
class Particao: public Problem<solPart>{
	private:
		std::vector<double> num;
	public:
		Particao(int n, unsigned sementeDaInstancia){
			std::mt19937 g(sementeDaInstancia);
			std::uniform_int_distribution<int> u(1, 1000000);
			for(int i = 0; i < n; i++) num.push_back(u(g));
		}

		// solucao inicial, ja avaliada; ptapi::randInt usa o gerador da replica
		solPart construction(){
			solPart s;
			for(size_t i = 0; i < num.size(); i++) s.lado.push_back((char)ptapi::randInt(0, 1));
			s.evalSol = evaluate(s);
			s.diff = diferenca(s);
			return s;
		}

		// vizinho: um numero troca de grupo
		solPart neighbor(solPart s){
			int i = ptapi::randInt(0, (int)num.size()-1);
			s.lado[i] = 1 - s.lado[i];
			s.diff = diferenca(s);
			return s;
		}

		// custo: |soma do grupo 0 - soma do grupo 1| (minimizacao)
		double evaluate(const solPart& s){ return std::fabs(diferenca(s)); }

		double diferenca(const solPart& s) const {
			double d = 0;
			for(size_t i = 0; i < num.size(); i++) d += s.lado[i] ? -num[i] : num[i];
			return d;
		}

		// Opcional: um passo sem copias e com avaliacao incremental, O(1) em vez de O(N)
		void neighborEval(const solPart& atual, solPart& cand){
			int i = ptapi::randInt(0, (int)num.size()-1);
			cand.lado = atual.lado;                       // reaproveita a memoria de 'cand'
			cand.lado[i] = 1 - cand.lado[i];
			cand.diff = atual.diff + (atual.lado[i] ? 2.0 : -2.0)*num[i];
			cand.evalSol = std::fabs(cand.diff);
		}
};

int main(){
	Particao prob(60, 2026);

	// Passo 3: o PT.
	// PT(Tmin, Tmax, replicas, MCL, PTL, distribuicao, tipoDeAtualizacao, intervaloDeAtualizacao)
	PT<solPart> pt(10.0f, 1000000.0f, 16, 400, 2000, 4, 0, 20);
	pt.setSeed(42);                                       // execucao reprodutivel
	pt.setTimeLimit(5.0);                                 // no maximo 5 segundos

	solPart melhor = pt.start(0, &prob);                  // 0 = todos os nucleos

	printf("melhor diferenca: %.0f (rodada %d de %d, %.1f ms de %.1f ms)\n", melhor.evalSol,
		pt.getBestPTLSolIndex(), pt.getPTLEnd(), pt.getTimeToBestMs(), pt.getElapsedMs());
	printf("%lld avaliacoes em %d threads, semente %llu\n", pt.getEvaluations(), pt.getThreadsUsed(), pt.getSeed());
	printf("idas e voltas: %lld   barreira: %.2f\n", pt.getRoundTrips(), pt.getCommunicationBarrier());
	printf("escada final :"); for(double x : pt.getTemperatures()) printf(" %.0f", x); printf("\n");
	printf("taxa de troca:"); for(double x : pt.getSwapAcceptanceRates()) printf(" %.2f", x); printf("\n");
	return std::fabs(prob.evaluate(melhor) - melhor.evalSol) < 1e-6 ? 0 : 1;
}
