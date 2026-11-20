#include "IPMTC.h"

//

IPMTC::IPMTC(std::string filename){
	
	std::string line; 
    std::ifstream ifs;
	int i=0;
	int j=0;
	ifs.open(filename);		

	if ( ifs.is_open()){
		
		getline(ifs,line);
		std::stringstream s(line);
		getline (s,line, ' ');
		numberMachine = std::stoi(line);
		std::cout<<"Machine number:"<< numberMachine<<"\n";
		getline (s,line, ' ');		
		numberJobs = std::stoi(line);
		std::cout<<"Jobs:"<< numberJobs<<"\n";
		getline(s,line, ' ');		
		numberTools = std::stoi(line);
		std::cout<<"Tools:"<< numberTools<<"\n";
		getline(s,line, ' ');		
		capacityMagazine = std::stoi(line);
		std::cout<<"Magazine:"<< capacityMagazine<<"\n";
		
		getline(ifs,line);
		changeTax = std::stoi(line);
		getline(ifs,line);		
		std::stringstream ss(line);

		toolJob.assign(numberTools, std::vector<bool>(numberJobs,false));
		JobTools.resize(numberJobs);
		
		while ( getline(ss,line, ' ')){
			totalJobTax += std::stoi(line);
			JobTax.push_back(std::stoi(line));
		}
					
		while(getline(ifs,line)){
			
			std::stringstream sss(line);
			while ( getline (sss ,line, ' ')){
				toolJob[i][j] = std::stoi(line);
				if(toolJob[i][j]){ 
					JobTools[j].push_back(i);
					++sum;
					}
			j++;
			}
		j=0;
		i++;		
		}
	ifs.close();
	
	
	// Print instance informations 
	std::cout<<"Change tax:"<< changeTax<<"\n";
	std::cout<<"Job tax:";
	for(auto& c:JobTax) std::cout<<c<<" ";
	std::cout<<"\n";
	std::cout<<"Jobs x Tools:\n";
	
	for(const auto &t:toolJob){
		for(const auto &jj:t) std::cout<<jj<<" ";
		std::cout<<"\n";
	}
	
	
	}else{
		std::cout << "Could not open! \n";
	}		
}

IPMTC::~IPMTC(){	
}

solIPMTC IPMTC::construction(){
	
	solIPMTC ss;
	// v2: random engine of the replica (fast and reproducible with PT::setSeed)
	ptapi::Rng& mersenne_engine = ptapi::rng(); 
	
	for(int i=0; i< numberJobs; i++){
		ss.sol.push_back(i);
	}
	
	std::shuffle(begin(ss.sol), end(ss.sol), mersenne_engine);
	
	int taxa = numberJobs/numberMachine;
	ss.cutStart.push_back(0);
	
	for(int np = 1; np < numberMachine; ++np){
		int p = np * taxa;
		ss.cutEnd.push_back(p);
		ss.cutStart.push_back(p);
	}
	
	ss.cutEnd.push_back(numberJobs);
	ss.evalSol = evaluate(ss);

	return ss;
}

// 2-opt neighbor moviment
solIPMTC IPMTC::neighbor(solIPMTC sol){
	solIPMTC& s = sol;   // v2: the argument is already a copy; it is changed and returned

	int first = 0;
	int last = 0;
	// v2: random engine of the replica (fast and reproducible with PT::setSeed)
	ptapi::Rng& mersenne_engine = ptapi::rng();
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
	
	return sol;
}

double IPMTC::evaluate(const solIPMTC& sol){
	
	int pipes_count = 0;
	int last_full = 0;
	std::vector<int> last_seen(numberTools);
	std::vector<std::vector<int>> M;
	

	//Completa o last_seen
	for(int i = 0; i < numberTools; ++i){
		if(toolJob[i][sol.sol[0]]) last_seen[i] = 0;
		else last_seen[i] = -1;
	}
	M.push_back(JobTools[sol.sol[0]]);
	
	
	for(int e = 1; e < numberJobs; ++e){
		
		M.push_back(JobTools[sol.sol[e]]);
				
		for (auto t = JobTools[sol.sol[e]].begin(); t != JobTools[sol.sol[e]].end(); ++t){
			
			if(last_full <= last_seen[*t]){
				
				++pipes_count;
								
				for(int i = (last_seen[*t]+1); i < e; ++i){
					
					M[i].push_back(*t);					
					if((int)M[i].size() == capacityMagazine) last_full = i;
				}
			}
			last_seen[*t] = e; 	
		}
		
		if((int)M[e].size() == capacityMagazine) last_full = e;

	}
 
	float result = (sum - capacityMagazine - pipes_count);
	if  (result < 0.0) result = 0.0;

	return result; 
}

solIPMTC IPMTC::setCut(solIPMTC sol){
	
	solIPMTC best = sol;
	
	float lowerBound = (float)(totalJobTax + (best.evalSol * changeTax))/(float)numberMachine;
	int job = 0;
	for(int np = 1; np < numberMachine; ++np){
	
		int makespamJ = 0;
		int makespamT = 0;
		std::vector<int> s;
	
		do{
			makespamJ += JobTax[best.sol[job]];
			s.push_back(best.sol[job]);
			makespamT = changeTax * GPCAIPMTC(s);
			++job;
		}while((makespamJ + makespamT) < lowerBound);
	
		best.cutStart[np] = job;
		best.cutEnd[np -1] = job;
	}


	best.evalSol = INT_MAX; 
	solIPMTC s = best;
	
	int max = 0;
	int cm = -1;
	
	while(cm != (numberMachine-1)){
		max = 0;
		cm = -1;
		
		for(int i = 0; i < numberMachine; ++i){
			std::vector<int> sv = std::vector<int>(s.sol.begin() + s.cutStart[i], s.sol.begin() + s.cutEnd[i]);		
			int tax = GPCAIPMTC(sv);				
			if(tax > max){
				max = tax;
				cm = i;
			}
		}
		
		if(best.evalSol > max){
			s.evalSol = max;
			best = s;
		}
		
		if((cm+1) < numberMachine){
			--s.cutEnd[cm]; 
			--s.cutStart[cm+1];
		}
		
	}	
		s = best;
		
		while(cm != 0){
		max = 0;
		cm = -1;
		
		for(int i = 0; i < numberMachine; ++i){
			std::vector<int> sv = std::vector<int>(s.sol.begin() + s.cutStart[i], s.sol.begin() + s.cutEnd[i]);		
			int tax = GPCAIPMTC(sv);				
			if(tax > max){
				max = tax;
				cm = i;
			}
		}
	
		if(best.evalSol > max){
			s.evalSol = max;
			best = s;
		}
		
		if((cm-1) >= 0){
			++s.cutEnd[cm-1]; 
			++s.cutStart[cm];
		}
		
	}

	return best;
}

int IPMTC::GPCAIPMTC(std::vector<int> sol){
	
	int pipes_count = 0;
	int custo = 0;
	int somaUm = 0;
	int last_full = 0;
	int jobSize = sol.size();
	
	std::vector<int> last_seen(numberTools);
	std::vector<std::vector<int>> M;
	

	//Completa o last_seen
	for(int i = 0; i < numberTools; ++i){
		if(toolJob[i][sol[0]]) last_seen[i] = 0;
		else last_seen[i] = -1;
	}
	M.push_back(JobTools[sol[0]]);
	custo += JobTax[sol[0]];
	somaUm += JobTools[sol[0]].size();
	
	
	for(int e = 1; e < jobSize; ++e){
		custo += JobTax[sol[e]];
		somaUm += JobTools[sol[e]].size();

		M.push_back(JobTools[sol[e]]);
				
		for (auto t = JobTools[sol[e]].begin(); t != JobTools[sol[e]].end(); ++t){
			
			if(last_full <= last_seen[*t]){
				
				++pipes_count;
								
				for(int i = (last_seen[*t]+1); i < e; ++i){
					
					M[i].push_back(*t);					
					if((int)M[i].size() == capacityMagazine) last_full = i;
				}
			}
			last_seen[*t] = e; 	
		}
		
		if((int)M[e].size() == capacityMagazine) last_full = e;

	}
	
 int trocas = (somaUm - capacityMagazine - pipes_count);
 if	(trocas < 0) trocas = 0;
 
 return (custo+(changeTax * trocas)); 
}

