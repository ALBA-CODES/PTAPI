
#ifndef __IPMTC_H__
#define __IPMTC_H__

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


struct solIPMTC: public solution{ 
  std::vector<int> sol; 
  std::vector<int> cutEnd; 
  std::vector<int> cutStart; 
};


class IPMTC: public Problem<solIPMTC>{
	private:

		int numberJobs = 0;
		int numberTools = 0;
		int capacityMagazine = 0;
		int numberMachine = 0;
		int changeTax = 0;
		int totalJobTax = 0;
		int sum = 0;
		std::string fn;
		std::vector<std::vector<bool>>  toolJob;
		std::vector<std::vector<int>>  JobTools;
		std::vector<int>  JobTax;
		
		
	public:
		IPMTC(std::string filename);
		~IPMTC();
		solIPMTC construction();
		solIPMTC neighbor(solIPMTC sol);
		double evaluate(const solIPMTC& sol); // v2: no copy of the solution
		solIPMTC setCut(solIPMTC sol);
		int GPCAIPMTC(std::vector<int> sol);

};

#endif 
