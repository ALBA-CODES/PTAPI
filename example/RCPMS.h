
#ifndef __RCPMS_H__
#define __RCPMS_H__

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


struct solRCPMS: public solution{ 
  std::vector<int> sol;
  std::vector<int> cutEnd; 
  std::vector<int> cutStart;
  std::vector<std::pair<int,int>> blockPos;  
};


class RCPMS: public Problem<solRCPMS>{
	private:

		int numberJobs = 0;
		int numberMachine = 0;
		int numberTools = 0;
		int changeTax = 0;
		int sum = 0;
		std::string fn;
		std::vector<int>  JobTools;
		std::vector<int>  JobTax;
		
	public:
		RCPMS(std::string filename);
		~RCPMS();
		solRCPMS construction();
		solRCPMS neighbor(solRCPMS sol);
		double evaluate(const solRCPMS& sol); // v2: no copy of the solution
		int criticalMachine(solRCPMS* sol);
};

#endif 
