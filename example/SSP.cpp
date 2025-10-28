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
    std::random_device rnd_device;
    std::mt19937 mersenne_engine {rnd_device()}; 
    	
	for(int i=0; i< numberJobs; i++){
		ss.sol.push_back(i);
	}
	
	std::shuffle(begin(ss.sol), end(ss.sol), mersenne_engine);
	
	ss.evalSol = evaluate(ss);
	ss.Nup = false;
	ss.Ndown = false;
	 
	return ss;
}

// 2-opt neighbor moviment
solSSP SSP::neighbor(solSSP sol){
	solSSP s;
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

double SSP::evaluate(solSSP sol){
	
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

