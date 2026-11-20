#include "SSP.h"

//

SSP::SSP(std::string filename){
	
	fn = filename;
	std::string line; 
    std::ifstream ifs;
	int i=0;
	int j=0;
	ifs.open(filename);		

	if ( ifs.is_open()){
		
		getline(ifs,line);
		numberJobs = std::stoi(line);
		
		getline(ifs,line);
		numberTools = std::stoi(line);
		
		getline(ifs,line);
		capacityMagazine = std::stoi(line);		
		toolJob.assign(numberTools, std::vector<bool>(numberJobs,false));
		JobTools.resize(numberJobs);
			
		while(getline(ifs,line)){
			
			std::stringstream ss(line);
			while ( getline(ss ,line, ' ')){
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
	
	}else{
		std::cout << "Could not open file! \n";
	}	
}

SSP::~SSP(){	
}

solSSP SSP::construction(){
	solSSP ss;
    // v2: random engine of the replica (fast and reproducible with PT::setSeed)
    ptapi::Rng& mersenne_engine = ptapi::rng(); 
    	
	for(int i=0; i< numberJobs; i++){
		ss.sol.push_back(i);
	}
	
	std::shuffle(begin(ss.sol), end(ss.sol), mersenne_engine);
	
	ss.evalSol = evaluate(ss);
	 
	return ss;
}

// 2-opt neighbor moviment
solSSP SSP::neighbor(solSSP sol){
	solSSP& s = sol;   // v2: the argument is already a copy; it is changed and returned
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

double SSP::evaluate(const solSSP& sol){
	
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

