#include "RCPMS.h"

//

RCPMS::RCPMS(std::string filename){

	std::string line; 
    std::ifstream ifs;
	ifs.open(filename);		

	if ( ifs.is_open()){
		
		getline(ifs,line);
		std::stringstream s(line);
		getline (s,line, ' ');
		numberJobs = std::stoi(line);
		getline (s,line, ' ');		
		numberMachine = std::stoi(line);		
		getline(s,line, ' ');		
		numberTools = std::stoi(line);
		
		getline(ifs,line);
		changeTax = std::stoi(line);
		
		getline(ifs,line);		
		std::stringstream ss(line);		
		while ( getline(ss,line, ' ')){
			JobTools.push_back(std::stoi(line));
		}
		
		getline(ifs,line);
		std::stringstream sss(line);		
		while ( getline(sss,line, ' ')){
			JobTax.push_back(std::stoi(line));
		}
					
	ifs.close();
	
	
	std::cout<<"Machines:"<< numberMachine<<"\n";
	std::cout<<"Jobs:"<< numberJobs<<"\n";
	std::cout<<"Tools:"<< numberTools<<"\n";
	std::cout<<"Change tax:"<< changeTax<<"\n";
	std::cout<<"Job tax:";
	for(auto& c:JobTax) std::cout<<c<<" ";
	std::cout<<"\n";
	std::cout<<"Jobs x Tools:";
	for(auto& c:JobTools) std::cout<<c<<" ";
	std::cout<<"\n";	
	
	}else{
		std::cout << "Could not open! \n";
	}	
}

RCPMS::~RCPMS(){	
}

solRCPMS RCPMS::construction(){
	
	solRCPMS ss;
	std::vector<int> v(numberTools);
	// v2: random engine of the replica (fast and reproducible with PT::setSeed)
	ptapi::Rng& mersenne_engine = ptapi::rng();
	
	// create vector
	std::iota(v.begin(), v.end(), 0);	
	std::shuffle(v.begin(), v.end(), mersenne_engine);
	
	int posFirst = 0;
	int posEnd = 0;
	
	for (int &t: v) {
		for(int i= 0; i < numberJobs; i++){
			if(JobTools[i] == t){
				ss.sol.push_back(i);
				++posEnd;
			}
		}
		ss.blockPos.push_back( std::make_pair(posFirst,(posEnd-1)));
		posFirst = posEnd;
	}

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
solRCPMS RCPMS::neighbor(solRCPMS sol){
	
	solRCPMS s;	
	s.sol = sol.sol;
	s.cutEnd = sol.cutEnd;
	s.cutStart = sol.cutStart;

	
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
	

	solRCPMS best = s;
	
	int cm = -1;
	
	while(cm != (numberMachine-1)){		

		cm = criticalMachine(&s);				
	
		if(best.evalSol > s.evalSol){
			best = s;
		}
		
		if((cm+1) < numberMachine){
			--s.cutEnd[cm]; 
			--s.cutStart[cm+1];
		}
	}

	s = best;
	
	while(cm != 0){
		
		cm = criticalMachine(&s);				
	
		if(best.evalSol > s.evalSol){
			best = s;
		}
		
		if((cm-1) >= 0){
			++s.cutEnd[cm-1]; 
			++s.cutStart[cm];
		}		
	}	
			
	return best;	
}

double RCPMS::evaluate(const solRCPMS& sol){
	
	int max = 0;
	int m = 0;
	int cm = -1;	
	
	std::queue<int> exec;

	std::vector<int> mSpam(numberMachine,0); //makespan da maquina
	std::vector<int> magazine(numberMachine,0);	//magazine atual da maquina
	std::vector<int> p(numberMachine,0);	
	std::vector<int> fSpamStart(numberTools,0);
	std::vector<int> fSpamEnd(numberTools,0);
	std::vector<int> fMachine(numberTools,-1);
	
	
	for(int i = 0; i < numberMachine;++i){
		magazine[i] = JobTools[sol.sol[sol.cutStart[i]]];
		exec.push(i);
	}

	
	while(!exec.empty()){
		m = exec.front();
		exec.pop();
							
		// verifica se ainda tem tarefas a ser processada
		if((sol.cutStart[m]+p[m]) < sol.cutEnd[m]){
				
			// verifica se o makespam onde a ferramenta foi utilizada é maior que o makespam da maquina
			if((fSpamStart[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] <= mSpam[m]) && (fSpamEnd[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] > mSpam[m])){
				// o makespam da maquina recebe o makespam da ferramenta
				mSpam[m] = fSpamEnd[JobTools[sol.sol[sol.cutStart[m]+p[m]]]];
			//	std::cout<<"Adionar gap \n";

			}else{				
				fSpamStart[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] = mSpam[m];
			}
			
			// Verifica se há troca
			if((magazine[m] != JobTools[sol.sol[sol.cutStart[m]+p[m]]]) || ((fMachine[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] != m) && (fMachine[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] != -1))){
				mSpam[m] += changeTax;
			}
				
			mSpam[m] += JobTax[sol.sol[sol.cutStart[m]+p[m]]];
				
			magazine[m] = JobTools[sol.sol[sol.cutStart[m]+p[m]]];
			fSpamEnd[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] = mSpam[m];
			fMachine[JobTools[sol.sol[sol.cutStart[m]+p[m]]]] = m;
				
		}
		
		++p[m];
	
		if(exec.empty()){
			int min = std::numeric_limits<int>::max();
			int machine = -1;
		
			for(int i = 0; i < numberMachine;++i){
				if((min >= mSpam[i]) && ((sol.cutStart[i]+p[i]) < sol.cutEnd[i])){
					machine = i;
					min = mSpam[i];
				}
				if(max < mSpam[i]){ 
					max = mSpam[i];
					cm = i;
				} 			
			}
			
			if(machine != -1) exec.push(machine);
		}
	}
	
		
	return max;
}

int RCPMS::criticalMachine(solRCPMS* sol){
	
	int max = 0;
	int m = 0;
	int cm = -1;	
	
	std::queue<int> exec;

	std::vector<int> mSpam(numberMachine,0); //makespan da maquina
	std::vector<int> magazine(numberMachine,0);	//magazine atual da maquina
	std::vector<int> p(numberMachine,0);	
	std::vector<int> fSpamStart(numberTools,0);
	std::vector<int> fSpamEnd(numberTools,0);
	std::vector<int> fMachine(numberTools,-1);
	
	
	for(int i = 0; i < numberMachine;++i){
		magazine[i] = JobTools[(*sol).sol[(*sol).cutStart[i]]];
		exec.push(i);
	}

	
	while(!exec.empty()){
		m = exec.front();
		exec.pop();
							
		// verifica se ainda tem tarefas a ser processada
		if(((*sol).cutStart[m]+p[m]) < (*sol).cutEnd[m]){
				
			// verifica se o makespam onde a ferramenta foi utilizada é maior que o makespam da maquina
			if((fSpamStart[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] <= mSpam[m]) && (fSpamEnd[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] > mSpam[m])){
				// o makespam da maquina recebe o makespam da ferramenta
				mSpam[m] = fSpamEnd[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]];
			//	std::cout<<"Adionar gap \n";

			}else{				
				fSpamStart[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] = mSpam[m];
			}
			
			// Verifica se há troca
			if((magazine[m] != JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]) || ((fMachine[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] != m) && (fMachine[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] != -1))){
				mSpam[m] += changeTax;
			}
				
			mSpam[m] += JobTax[(*sol).sol[(*sol).cutStart[m]+p[m]]];
				
			magazine[m] = JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]];
			fSpamEnd[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] = mSpam[m];
			fMachine[JobTools[(*sol).sol[(*sol).cutStart[m]+p[m]]]] = m;
				
		}
		
		++p[m];
	
		if(exec.empty()){
			int min = std::numeric_limits<int>::max();
			int machine = -1;
		
			for(int i = 0; i < numberMachine;++i){
				if((min >= mSpam[i]) && (((*sol).cutStart[i]+p[i]) < (*sol).cutEnd[i])){
					machine = i;
					min = mSpam[i];
				}
				if(max < mSpam[i]){ 
					max = mSpam[i];
					cm = i;
				} 			
			}
			
			if(machine != -1) exec.push(machine);
		}
	}
	
	
	(*sol).evalSol = max;
	
	return cm;
}


