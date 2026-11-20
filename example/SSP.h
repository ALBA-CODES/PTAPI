
#ifndef __SSP_H__
#define __SSP_H__

#include <string.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <queue>
#include <utility>
#include <forward_list>
#include <list>
#include <set>
#include <map>
#include <limits.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <thread>
#include <atomic>
#include "../include/Problem.h"


struct solSSP: public solution{ 
  std::vector<int> sol; 
};


class SSP: public Problem<solSSP>{
	private:

		int numberJobs;
		int numberTools;
		int capacityMagazine;
		int sum = 0;
		std::string fn;
		std::vector<std::vector<bool>>  toolJob;
		std::vector<std::vector<int>>  JobTools;

		
	public:
		SSP(std::string filename);
		~SSP();
		solSSP construction();
		solSSP neighbor(solSSP sol);
		double evaluate(const solSSP& sol); // v2: no copy of the solution

};

#endif 
