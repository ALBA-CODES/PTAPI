#include "SD.h"

//
thread_local std::vector<std::pair <int,int>> SD::twoSwap;

SD::SD(std::string filename, int neighborType, int stype){
	
	type = stype;
	nType = neighborType;
	fn = filename;
	std::string line; 
    std::ifstream ifs;
	//int i=0;
	int jo = 0;
	ifs.open(filename);		

	if ( ifs.is_open()){
		
		getline(ifs,line);
		std::stringstream data(line);
		getline (data ,line, ' ');		
		numberJobs = std::stoi(line);
		getline (data ,line, ' ');		
		numberTools = std::stoi(line);
		getline (data ,line, ' ');				
		capacityMagazine = std::stoi(line);		
		getline (data ,line, ' ');				
		capacityMin = std::stoi(line);		
		getline (data ,line, ' ');				
		capacityMax = std::stoi(line);		
		
		toolJob.assign(numberTools, std::vector<bool>(numberJobs,false));
		changeTax.assign(numberTools, std::vector<int>(numberTools,0));
		JobTools.resize(numberJobs);
		iniTax.resize(numberTools);
		
		for(int i = 0;i < numberTools; ++i){
			getline(ifs,line);
			iniTax[i] = std::stoi(line);
		}
		
		for(int i = 0;i < numberTools; ++i){
			getline(ifs,line);
			std::stringstream ss(line);
			for(int j = 0;j < numberTools; ++j){
				getline (ss ,line, ' ');
				changeTax[i][j] = std::stoi(line);
			}
		}
		
		while(getline(ifs,line)){
			std::stringstream ss(line);
			while ( getline (ss ,line, ' ')){
				int t = std::stoi(line);
				if(t > 0){
						toolJob[(t-1)][jo] = true;
						JobTools[jo].push_back(t-1);
						++sum;
					}
				jo++;
			}
			jo=0;
		}
		
	ifs.close();
	
	
	//~ std::cout<<"Jobs:"<< numberJobs<<"\n";
	//~ std::cout<<"Tools:"<< numberTools<<"\n";
	//~ std::cout<<"C total:"<< capacityMagazine<<"\n";
	//~ std::cout<<"C Min:"<< capacityMin<<"\n";
	//~ std::cout<<"C Max:"<< capacityMax<<"\n";

	//~ std::cout<<"Change tax ini \n";
	//~ for(auto& c:iniTax) std::cout<<c<<" ";
	//~ std::cout<<"\n";
	
	//~ std::cout<<"Change tax:\n";
	//~ for(auto& c:changeTax) {
		//~ for(auto& cc:c) std::cout<<cc<<" ";
		//~ std::cout<<"\n";
	//~ }
	
	//~ std::cout<<"\n";
	//~ std::cout<<"Jobs x Tools:\n";
	//~ for(auto& c:JobTools){
		//~ for(auto& cc:c) std::cout<<cc<<" ";
		//~ std::cout<<"\n";	
	//~ }
	//~ std::cout<<"\n";
	
	
	}else{
		std::cout << "Could not open file! \n";
	}	
}

SD::~SD(){	
}

solSD SD::construction(){
	solSD ss;
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; 
    	
	for(int i=0; i< numberJobs; i++){
		ss.sol.push_back(i);
	}
	
	std::shuffle(begin(ss.sol), end(ss.sol), mersenne_engine);
	
	ss = evaluate(ss);
	ss.Nup = false;
	ss.Ndown = false;
	 
	return ss;
}

solSD SD::neighbor(solSD sol){
	switch (nType)
	{
		case 1:{ 
			return neighbor2opt(sol); 
		}
		break;
		case 2:{ 
			return neighbor2Swap(sol); 
		}
		break;
		case 3:{ 
			return neighborRandIns(sol); 
		}
		break;
		default:
			return neighbor2opt(sol); 
		break;	
	}
}

// 2-opt neighbor moviment
solSD SD::neighbor2opt(solSD sol){
	
	solSD s;
	s.sol = sol.sol;
	int first = 0;
	int last = 0;
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};
	std::uniform_int_distribution<int> dist {0, (numberJobs-1)};
	
	do {		
		first = dist(mersenne_engine);
		last = dist(mersenne_engine);		
	} while(first == last);
	
	if (first > last)std::swap(first,last);
	
	while (first < last) {
		std::swap(s.sol[first], s.sol[last]);
		first++;
		last--;
	}
	
	s.Nup = sol.Nup;
	s.Ndown = sol.Ndown;

			
	return s;
}

solSD SD::neighborRandIns(solSD sol){
	
	solSD s;
	s.sol = sol.sol;
	int first = 0;
	int last = 0;
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};
	std::uniform_int_distribution<int> dist {0, (numberJobs-1)};
	
	
	do {		
		first = dist(mersenne_engine);
		last = dist(mersenne_engine);		
	} while(first == last);
	
	
	int value = s.sol[first];
	s.sol.erase(s.sol.begin()+first);
	
	if(last != (numberJobs-1)){
		s.sol.insert(s.sol.begin()+last, value);
	}else{
		s.sol.insert(s.sol.end(), value);	
	}

	
	s.Nup = sol.Nup;
	s.Ndown = sol.Ndown;
			
	return s;
}

solSD SD::neighbor2Swap(solSD sol){
	solSD s;
	s.sol = sol.sol;
	std::random_device rnd_device;
	std::mt19937 mersenne_engine {rnd_device()};

	
	if (twoSwap.size() == 0){
		for(int pos1 = 0; pos1 < numberJobs; pos1++){		
			for(int pos2 = (pos1+1); pos2 < numberJobs; pos2++){		
				twoSwap.push_back(std::make_pair (pos1,pos2));
			}
		}
	}

	std::shuffle(begin(twoSwap), end(twoSwap), mersenne_engine);	
	std::swap(s.sol[twoSwap[0].first], s.sol[twoSwap[0].second]);
		
	s.Nup = sol.Nup;
	s.Ndown = sol.Ndown;
			
	return s;
}


// KTNS proposto
solSD SD::evaluate(solSD sol){
	switch (type)
	{
		case 1:{ 
			return evaluate1(sol);
		}
		break;
		case 2:{ 
			return evaluate2(sol);
		}
		break;
		case 3:{ 
			return evaluate3(sol);
		}
		break;
		default:
			return evaluate2(sol);
		break;	
	}	
}


// KTNS+H Primeira parte é preenchida de forma aleatória
solSD SD::evaluate1(solSD sol){
	
	// Variaveis
	std::vector<bool> magazineL(numberTools, true);	
	double evalSol = 0.0; 
	
	// Percorre todas as ferramentas
	for(int jL= 0; jL < numberJobs; ++jL){
		
		std::vector<bool> magazineCL(numberTools,false);
		int left = jL;
		int cmL = 0;

		while((cmL < capacityMagazine) && (left < numberJobs)){
			
				for (auto it=JobTools[sol.sol[left]].begin(); ((it!=JobTools[sol.sol[left]].end()) && (cmL < capacityMagazine)); ++it){
					if((magazineL[*it]) && (!magazineCL[*it])){

						// Add ao magazine atual
						magazineCL[*it] = true;
						++cmL;
					
						// Custo de troca inicial
						if(jL == 0) evalSol += iniTax[*it]; 
					
					}else if((jL == left) && (!magazineCL[*it])){
						magazineCL[*it] = true;
						//newTools.push_back(*it);
						//		++evalSol;
						++cmL;
					}
				}
//			}
			++left;
		}
		
		for(int t=0; ((t < numberTools) && (cmL < capacityMagazine)); t++){
			if((magazineL[t]) && (!magazineCL[t])){
				magazineCL[t] = true;
				//mag.push_back(t);
				++cmL;			
			}
		}		
		
		// --- 2. PREPARAÇÃO PARA O ALGORITMO HÚNGARO ---
        
        std::vector<int> insertedTools;
        std::vector<int> removedTools;
        
        // Identifica quem entrou e quem saiu
        for (int t = 0; t < numberTools; ++t) {
            bool inPrev = magazineL[t];
            bool inCurr = magazineCL[t];

            if (inCurr && !inPrev) {
                insertedTools.push_back(t);
            } else if (!inCurr && inPrev) {
                removedTools.push_back(t);
            }
        }
		
		
// --- CÁLCULO DE CUSTO OTIMIZADO ---
        
        // Se houver trocas para fazer (entradas substituindo saídas)
        if (!insertedTools.empty() && !removedTools.empty()) {
            
            size_t n = removedTools.size();
            
            std::vector<std::vector<double>> costMatrix(n, std::vector<double>(n));

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    int toolOut = removedTools[i];
                    int toolIn = insertedTools[j];
                    
                    // Preenche com o custo da matriz changeTax
                    costMatrix[i][j] = changeTax[toolOut][toolIn];
                }
            }

            // Resolve o problema de designação
            evalSol += hungarianSolve(costMatrix);         
			
			if(jL == 0) std::cout<<"ERRO\n";

        } 	
        
		magazineL = magazineCL;
	}

	sol.evalSol = evalSol;
	
	return sol;

}


// KTNS+H com primeiro utilizando o H
solSD SD::evaluate2(solSD sol){
	
	// Variaveis
	std::vector<bool> magazineL(numberTools, true);	
	std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; 
	std::vector<int> firststage(capacityMagazine, -1);
	std::vector<int> secondstage(capacityMagazine, -1);

	double evalSol = 0.0; 
	
	// Percorre todas as ferramentas
	for(int jL= 0; jL < numberJobs; ++jL){
		
		std::vector<bool> magazineCL(numberTools,false);
		std::vector<int> newTools;		
		std::vector<int> toolsSort;

		int left = jL;
		int cmL = 0;
		int cmLID = 0;
		int smL = 0;

		while((cmL < capacityMagazine) && (left < numberJobs)){

			if ((jL == 0) && (left > 0) && ((cmL+JobTools[sol.sol[left]].size()) > capacityMagazine)){
				toolsSort = JobTools[sol.sol[left]];
				for(auto it = toolsSort.begin(); ((it != toolsSort.end()) && (smL < capacityMagazine));++it){
					if(!magazineCL[*it]){	
						++cmL;
					}
					secondstage[smL] = *it; 
					++smL;	
				}
				
				if(cmL < capacityMagazine){
					
					for(auto it = toolsSort.begin(); ((it != toolsSort.end()) && (smL < capacityMagazine));++it){
						if(!magazineCL[*it]){	
							magazineCL[*it] = true;
							firststage[cmLID] = *it;
							++cmLID;
						}
					}					
					std::fill(secondstage.begin(), secondstage.end(), -1); 
					smL = 0;
				}				
				
			}else{
			
				for (auto it=JobTools[sol.sol[left]].begin(); ((it!=JobTools[sol.sol[left]].end()) && (cmL < capacityMagazine)); ++it){
					if((magazineL[*it]) && (!magazineCL[*it])){

						// Custo de troca inicial
						if(jL == 0){
//							evalSol += iniTax[*it];
							firststage[cmLID] = *it;
							++cmLID;
						}  

						// Add ao magazine atual
						magazineCL[*it] = true;
						++cmL;					
					
					}else if((jL == left) && (!magazineCL[*it])){
						magazineCL[*it] = true;
						newTools.push_back(*it);
						//		++evalSol;
						++cmL;
					}
				}
			}
			++left;
		}
		
		for(int t=0; ((t < numberTools) && (cmL < capacityMagazine)); t++){
			if((magazineL[t]) && (!magazineCL[t])){
				magazineCL[t] = true;
				//mag.push_back(t);
				++cmL;			
			}
		}		
		
		// --- 2. PREPARAÇÃO PARA O ALGORITMO HÚNGARO ---
        
        std::vector<int> insertedTools;
        std::vector<int> removedTools;
        
        // Identifica quem entrou e quem saiu
        for (int t = 0; t < numberTools; ++t) {
            bool inPrev = magazineL[t];
            bool inCurr = magazineCL[t];

            if (inCurr && !inPrev) {
                insertedTools.push_back(t);
            } else if (!inCurr && inPrev) {
                removedTools.push_back(t);
            }
        }
		
		
// --- CÁLCULO DE CUSTO OTIMIZADO ---
        
        // Se houver trocas para fazer (entradas substituindo saídas)
        if (!insertedTools.empty() && !removedTools.empty()) {
            
            size_t n = removedTools.size();
            
            std::vector<std::vector<double>> costMatrix(n, std::vector<double>(n));

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    int toolOut = removedTools[i];
                    int toolIn = insertedTools[j];
                    
                    // Preenche com o custo da matriz changeTax
                    costMatrix[i][j] = changeTax[toolOut][toolIn];
                }
            }

            // Resolve o problema de designação
            evalSol += hungarianSolve(costMatrix); 
       if(debug) std::cout<<"["<<jL<<"] "<<"Valor intermediario: " <<evalSol<<"\n"; 
                   
        } else if(jL == 0){
			
		if(debug)std::cout<<"left: "<<(left-1)<<"\n";
			
		if(debug){	
			for(size_t v=0; v < sol.sol.size(); ++v){
				std::cout<<"tarefa: "<<sol.sol[v]<<"\n"<<"Ferramentas:[";			
		 
					for (auto it=JobTools[sol.sol[v]].begin(); it!=JobTools[sol.sol[v]].end(); ++it){
						std::cout<< *it <<" ";			
					}
				std::cout<< "]" <<"\n";			
			}
		}

			std::vector<std::vector<double>> costMatrix(capacityMagazine, std::vector<double>(capacityMagazine));
			 for (int i = 0; i < capacityMagazine; ++i) {
                for (int j = 0; j < capacityMagazine; ++j) {
                    int toolOut = firststage[i];
                    int toolIn = secondstage[j];
                    
                    // Preenche com o custo da matriz changeTax
                    if ((toolOut == -1) && (toolIn == -1)){
						costMatrix[i][j] = std::numeric_limits<int>::max();
					}else if((toolOut != -1) && (toolIn != -1)){
						costMatrix[i][j] = changeTax[toolOut][toolIn];
					}else if((toolOut == -1) && (toolIn != -1)){
						if(magazineCL[toolIn]){
							costMatrix[i][j] = std::numeric_limits<int>::max();
						}else{
							costMatrix[i][j] = iniTax[toolIn];
						}					
					}else if((toolOut != -1) && (toolIn == -1)){
						costMatrix[i][j] = 0;										
					}
                    
                }
            }
			int valor = 0;
			 // Resolve o problema de designação
  //DESCOMENTAR SE FOR UTILIZAR
            magazineCL = hungarian2(costMatrix,firststage,secondstage,&valor);
			evalSol += valor;
			sol.firstM = magazineCL;
			jL = (left-1);
			
			if(debug){
				std::cout<<"Ferramentas Magazine:\n";
				for(int vv = 0; vv<numberTools;++vv){
					if(magazineCL[vv])std::cout<<vv<<"\n";
				}
			
			}
			
		if(debug) std::cout<<"valor:"<<evalSol<<"\n";
		}	
        
		
		magazineL = magazineCL;
		//magL = mag;
		//if(jL == 0) sol.firstM = magazineCL;
	}
	
	//std::cout<<"valor Total:"<<evalSol<<"\n";

	sol.evalSol = evalSol;
	
	return sol;

}


// Somente H
solSD SD::evaluate3(solSD sol){
	
	// Variaveis
	std::vector<int> lastStage(capacityMagazine, -1);
	double evalSol = 0.0;
	int lsmL = 0;
	
	if(debug){	for(size_t v=0; v < sol.sol.size(); ++v){
				std::cout<<"tarefa: "<<sol.sol[v]<<"\n"<<"Ferramentas:[";
			
					for (auto it=JobTools[sol.sol[v]].begin(); it!=JobTools[sol.sol[v]].end(); ++it){
						std::cout<< *it <<" ";			
					}
			
					std::cout<< "]" <<"\n";			
	}}
	
	// Percorre todas as ferramentas
	for(int jL= 0; jL < numberJobs; ++jL){
		
		std::vector<bool> magazineCL(numberTools,false);
		std::vector<int> CurrStage(capacityMagazine, -1);
		
		
		int cmL = 0;
		
		
		if(jL == 0){
			int left = jL;
			int csmL = 0;
		
			while((cmL < capacityMagazine) && (left < numberJobs)){
				std::vector<int> auxTools;
				for (auto it=JobTools[sol.sol[left]].begin(); (it!=JobTools[sol.sol[left]].end()); ++it){
					if(!magazineCL[*it]){
						++cmL;
						auxTools.push_back(*it);
					}
				}
				if(cmL < capacityMagazine){
					for (auto it=auxTools.begin(); (it!=auxTools.end()); ++it){
						magazineCL[*it] = true;
						lastStage[csmL] = *it;
						++csmL;
					}
				}else{
					jL = left;
					for (size_t i = 0; i < JobTools[sol.sol[jL]].size(); ++i){
						CurrStage[i] = JobTools[sol.sol[jL]][i];
					}		
				}
			++left;
			}
		}else{
			for (size_t i = 0; i < JobTools[sol.sol[jL]].size(); ++i){
				CurrStage[i] = JobTools[sol.sol[jL]][i];
			}		
		}	
		
        		
// --- CÁLCULO DE CUSTO OTIMIZADO ---
        
        
        // Se houver trocas para fazer (entradas substituindo saídas)
		if(jL > 0){

			if(debug)std::cout<<"left: "<<jL<<"\n";

			std::vector<std::vector<double>> costMatrix(capacityMagazine, std::vector<double>(capacityMagazine));

			 for (int i = 0; i < capacityMagazine; ++i) {
                for (int j = 0; j < capacityMagazine; ++j){
                    int toolOut = lastStage[i];
                    int toolIn = CurrStage[j];
                    
                    // Preenche com o custo da matriz changeTax
                    if ((toolOut == -1) && (toolIn == -1)){
						costMatrix[i][j] = std::numeric_limits<int>::max();
					}else if((toolOut != -1) && (toolIn != -1)){
						costMatrix[i][j] = changeTax[toolOut][toolIn];
					}else if((toolOut == -1) && (toolIn != -1)){
						if(magazineCL[toolIn]){
							costMatrix[i][j] = std::numeric_limits<int>::max();
						}else{
							costMatrix[i][j] = iniTax[toolIn];
						}												
					}else if((toolOut != -1) && (toolIn == -1)){
						costMatrix[i][j] = 0;										
					}
                    
                }
            }
			int valor = 0;
			 // Resolve o problema de designação      
  // COMENTADO SE FOR UTILIZAR DESCOMENTAR
            CurrStage = hungarian3(costMatrix,lastStage,CurrStage,lsmL,&valor);
			evalSol += valor;
			
			if(debug){
				std::cout<<"Ferramentas Magazine: ";
				for(int vv = 0; vv < capacityMagazine;++vv){
					std::cout<<CurrStage[vv]<<" ";
				}
				std::cout<<"\n";
			}
			
		if(debug) std::cout<<"valor:"<<evalSol<<"\n";
		
		} //FIM HUNGARO
		
		lastStage = CurrStage;
		lsmL = jL;
	
	} //FIM TODOS OS ESTAGIOS
	
	sol.evalSol = evalSol;
	
	return sol;
}


// Função auxiliar: Algoritmo Húngaro O(n^3)
double SD::hungarianSolve(std::vector<std::vector<double>> costMatrix) {
    int n = costMatrix.size();
    if (n == 0) return 0.0;
    int m = costMatrix[0].size();
    
    // Vetores auxiliares (1-based indexing para facilitar a lógica interna)
    std::vector<double> u(n + 1, 0), v(m + 1, 0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    const double INF = std::numeric_limits<double>::max();

    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<bool> used(m + 1, false);
        
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            
            for (int j = 1; j <= m; ++j) {
                if (!used[j]) {
                    double cur = costMatrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    
    // O custo mínimo é o valor acumulado (negativo do potencial)
    return -v[0];
}


void SD::print(solSD sol, std::string filename, float tempMin, float tempMax, int tempL, int MCL, int PTL, int time, int indexPTL, int indexEnd){
    
    std::string base_filename = filename.substr(filename.find_last_of("/\\") + 1);
	std::string::size_type const p(base_filename.find_last_of('.'));
	std::string file_without_extension = base_filename.substr(0, p);
	std::ofstream timeF;
	timeF.open("PT_"+file_without_extension+".txt");
	timeF << "Instance filename: " << filename << "\n";
	timeF << "Min temp: " << tempMin << "\n";
	timeF << "Max temp: " << tempMax << "\n";
	timeF << "Best sol: " << sol.evalSol << "\n";
	timeF << "Best sol PTL: " << indexPTL << "\n";
	timeF << "Temp div: " << tempL << "\n";
	timeF << "MCL: " << MCL << "\n";
	timeF << "PTL: " << PTL << "\n";
	timeF << "Elapsed time was: " << time << "\n";
	timeF << "End PTL: " << indexEnd << "\n";
	timeF << "Validate:\n";
	for(int j:sol.sol)timeF <<j<<" ";
	timeF << "\n";
	timeF.close();		
}

// KTNS validate
void SD::validateAndPrint(solSD sol, std::ofstream* file){
	
	// Variaveis
	std::vector<bool> magazineL(numberTools, true);	
	std::vector<int> magL;
	double evalSol = 0.0; 
	
	// Percorre todas as ferramentas
	for(int jL= 0; jL < numberJobs; ++jL){
		
		std::vector<bool> magazineCL(numberTools,false);
		std::vector<int> newTools;		
		std::vector<int> mag;
		
		int trocas = 0;
		int left = jL;
		int cmL = 0;

		if(jL == 0){
		for(int i = 0; i < numberTools; ++i){
			if(sol.firstM[i]){
				evalSol += iniTax[i]; 
			//	std::cout<<"Fini: "<<i<<" tax: "<<iniTax[i]<<"\n";
				magazineCL[i] = true;
				++cmL;
			}
		}
		
		}else{
		while((cmL < capacityMagazine) && (left < numberJobs)){
			for (auto it=JobTools[sol.sol[left]].begin(); ((it!=JobTools[sol.sol[left]].end()) && (cmL < capacityMagazine)); ++it){
				if((magazineL[*it]) && (!magazineCL[*it])){
					magazineCL[*it] = true;
					mag.push_back(*it);
					++cmL;
					
					// Custo de troca inicial
					//if(jL == 0) {
					//	evalSol += iniTax[*it]; 
					//	std::cout<<"Fini: "<<*it<<" tax: "<<iniTax[*it]<<"\n";
					//}
				}else if((jL == left) && (!magazineCL[*it])){
					magazineCL[*it] = true;
					mag.push_back(*it);
					newTools.push_back(*it);
					++trocas;
					++cmL;
				}
			}
			++left;
		}
		}
		for(int t=0; ((t < numberTools) && (cmL < capacityMagazine)); t++){
			if((magazineL[t]) && (!magazineCL[t])){
				magazineCL[t] = true;
				mag.push_back(t);
				++cmL;			
			}
		}		
		// --- 2. PREPARAÇÃO PARA O ALGORITMO HÚNGARO ---
        
        std::vector<int> insertedTools;
        std::vector<int> removedTools;
        
        // Identifica quem entrou e quem saiu
        for (int t = 0; t < numberTools; ++t) {
            bool inPrev = magazineL[t];
            bool inCurr = magazineCL[t];

            if (inCurr && !inPrev) {
                insertedTools.push_back(t);
            } else if (!inCurr && inPrev) {
                removedTools.push_back(t);
            }
        }
		
		(*file) << trocas << " ";
		
// --- CÁLCULO DE CUSTO OTIMIZADO ---
        
        // Se houver trocas para fazer (entradas substituindo saídas)
        if (!insertedTools.empty() && !removedTools.empty()) {
            
            size_t n = removedTools.size();
            
            std::vector<std::vector<double>> costMatrix(n, std::vector<double>(n));

            for (size_t i = 0; i < n; ++i) {
                for (size_t j = 0; j < n; ++j) {
                    int toolOut = removedTools[i];
                    int toolIn = insertedTools[j];
                    
                    // Preenche com o custo da matriz changeTax
                    costMatrix[i][j] = changeTax[toolOut][toolIn];
                }
            }

            // Resolve o problema de designação
            evalSol += hungarianSolvePrint(costMatrix, removedTools, insertedTools, file);        
        } else{
	//	std::cout<<"Valor inicial: " <<evalSol<<" cml:"<<cmL<<"\n";	
		}	
		(*file) << "\n";
		magazineL = magazineCL;
		magL = mag;
	}
	//std::cout<<"Valor Final: " <<evalSol<<"\n";
}


double SD::hungarianSolvePrint(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID, std::ofstream* file) 
{
    int n = costMatrix.size();
    if (n == 0) return 0.0;
    int m = costMatrix[0].size();
    
    std::vector<double> u(n + 1, 0), v(m + 1, 0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    const double INF = std::numeric_limits<double>::max();

    // --- Lógica do Algoritmo Húngaro (Inalterada) ---
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<bool> used(m + 1, false);
        
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            
            for (int j = 1; j <= m; ++j) {
                if (!used[j]) {
                    double cur = costMatrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    
    // --- IMPRESSÃO DOS RESULTADOS DENTRO DA FUNÇÃO ---
  //  std::cout << ">>> Solucao Otima Hungaro (Custo Total: " << -v[0] << "):\n";
    
    // Percorre as colunas (j) para ver qual linha (p[j]) foi associada
    for (int j = 1; j <= m; ++j) {
        int linhaEscolhida = p[j]; // Índice 1-based
        
        // Converte para 0-based para acessar os vetores e a matriz
        int rowIdx = linhaEscolhida - 1; 
        int colIdx = j - 1;

        // Pega os IDs reais
        int toolOut = rowsID[rowIdx];
        int toolIn  = colsID[colIdx];
//        double custoIndividual = costMatrix[rowIdx][colIdx];
        
        (*file) << toolOut<<" "<< toolIn<<" ";

//        std::cout << "    [Sai Ferramenta " << toolOut << "]"
 //                 << " substituida por [Entra Ferramenta " << toolIn << "]"
 //                 << " | Custo: " << custoIndividual << "\n";
    }
//   (*file) << "\n";
    
   // std::cout << "------------------------------------------------\n";

    return -v[0];
}

std::vector<int> SD::hungarian3(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID,int lsml, int *value) 
{
    int n = costMatrix.size();
    int m = costMatrix[0].size();
	std::vector<int> magazineCL;
//	std::vector<bool> magazineCL(numberTools,false);
    
    if(debug)std::cout<<"n: "<<n<<" m:"<<m<<" rowsID:"<<rowsID.size()<<" colsID: "<<colsID.size()<<"\n";
    
    std::vector<double> u(n + 1, 0), v(m + 1, 0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    const double INF = std::numeric_limits<double>::max();

    // --- Lógica do Algoritmo Húngaro (Inalterada) ---
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<bool> used(m + 1, false);
        
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            
            for (int j = 1; j <= m; ++j) {
                if (!used[j]) {
                    double cur = costMatrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    
    // --- IMPRESSÃO DOS RESULTADOS DENTRO DA FUNÇÃO ---
  if(debug) std::cout << ">>> Solucao Otima Hungaro (Custo Total: " << -v[0] << "):\n";
    
    // Percorre as colunas (j) para ver qual linha (p[j]) foi associada
    int soma = 0;
    int somatroca = 0;

    for (int j = 1; j <= m; ++j) {
        int linhaEscolhida = p[j]; // Índice 1-based
        
        // Converte para 0-based para acessar os vetores e a matriz
        int rowIdx = linhaEscolhida - 1; 
        int colIdx = j - 1;

        // Pega os IDs reais
        int toolOut = rowsID[rowIdx];
        int toolIn  = colsID[colIdx];
        double custoIndividual = costMatrix[rowIdx][colIdx];        

   if(debug)    std::cout << "    [Sai Ferramenta " << toolOut << "]"
                  << " substituida por [Entra Ferramenta " << toolIn << "]"
                  << " | Custo: " << custoIndividual << "\n";
                  
       //CUSTO INICIAL DE CARREGAMENTO
	if(lsml == 0){	
		if(toolOut != -1){
			soma+= iniTax[toolOut];
		}else{
			soma+= iniTax[toolIn];
		}
	}
		if(toolIn != -1){
		   magazineCL.push_back(toolIn);
//		   magazineCL[toolIn] = true;
		}else{
			magazineCL.push_back(toolOut);
//		   magazineCL[toolOut] = true;
		}

		
		if((toolOut != -1) && (toolIn != -1)){
			somatroca+=changeTax[toolOut][toolIn];
		}	
    }
 if(debug)   std::cout << "Custo inicial: "<<soma<<" Custo troca: "<<somatroca<<"\n"; 
    *value = soma+somatroca;   
 if(debug)   std::cout << "------------------------------------------------\n";

   // return -v[0];
	  return magazineCL;
}


std::vector<bool> SD::hungarian2(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID, int *value) 
{
    int n = costMatrix.size();
    int m = costMatrix[0].size();
//	std::vector<int> magazineCL;
	std::vector<bool> magazineCL(numberTools,false);
    
    if(debug)std::cout<<"n: "<<n<<" m:"<<m<<" rowsID:"<<rowsID.size()<<" colsID: "<<colsID.size()<<"\n";
    
    std::vector<double> u(n + 1, 0), v(m + 1, 0);
    std::vector<int> p(m + 1, 0), way(m + 1, 0);
    const double INF = std::numeric_limits<double>::max();

    // --- Lógica do Algoritmo Húngaro (Inalterada) ---
    for (int i = 1; i <= n; ++i) {
        p[0] = i;
        int j0 = 0;
        std::vector<double> minv(m + 1, INF);
        std::vector<bool> used(m + 1, false);
        
        do {
            used[j0] = true;
            int i0 = p[j0];
            int j1 = 0;
            double delta = INF;
            
            for (int j = 1; j <= m; ++j) {
                if (!used[j]) {
                    double cur = costMatrix[i0 - 1][j - 1] - u[i0] - v[j];
                    if (cur < minv[j]) {
                        minv[j] = cur;
                        way[j] = j0;
                    }
                    if (minv[j] < delta) {
                        delta = minv[j];
                        j1 = j;
                    }
                }
            }
            
            for (int j = 0; j <= m; ++j) {
                if (used[j]) {
                    u[p[j]] += delta;
                    v[j] -= delta;
                } else {
                    minv[j] -= delta;
                }
            }
            j0 = j1;
        } while (p[j0] != 0);
        
        do {
            int j1 = way[j0];
            p[j0] = p[j1];
            j0 = j1;
        } while (j0);
    }
    
    // --- IMPRESSÃO DOS RESULTADOS DENTRO DA FUNÇÃO ---
  if(debug) std::cout << ">>> Solucao Otima Hungaro (Custo Total: " << -v[0] << "):\n";
    
    // Percorre as colunas (j) para ver qual linha (p[j]) foi associada
    int soma = 0;
    int somatroca = 0;

    for (int j = 1; j <= m; ++j) {
        int linhaEscolhida = p[j]; // Índice 1-based
        
        // Converte para 0-based para acessar os vetores e a matriz
        int rowIdx = linhaEscolhida - 1; 
        int colIdx = j - 1;

        // Pega os IDs reais
        int toolOut = rowsID[rowIdx];
        int toolIn  = colsID[colIdx];
        double custoIndividual = costMatrix[rowIdx][colIdx];        

   if(debug)    std::cout << "    [Sai Ferramenta " << toolOut << "]"
                  << " substituida por [Entra Ferramenta " << toolIn << "]"
                  << " | Custo: " << custoIndividual << "\n";
                  
       //CUSTO INICIAL DE CARREGAMENTO
		if(toolOut != -1){
			soma+= iniTax[toolOut];
		}else{
			soma+= iniTax[toolIn];
		}

		if(toolIn != -1){
		   magazineCL[toolIn] = true;
		}else{
		   magazineCL[toolOut] = true;
		}

		
		if((toolOut != -1) && (toolIn != -1)){
			somatroca+=changeTax[toolOut][toolIn];
		}	
    }
 if(debug)   std::cout << "Custo inicial: "<<soma<<" Custo troca: "<<somatroca<<"\n"; 
    *value = soma+somatroca;   
 if(debug)   std::cout << "------------------------------------------------\n";

   // return -v[0];
	  return magazineCL;
}

void SD::setDebug(bool value){
	debug = value;
}
