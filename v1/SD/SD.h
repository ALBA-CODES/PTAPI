
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
#include <omp.h>
#include <limits.h>
#include <random>
#include <cmath>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <thread>
#include <atomic>
#include "../include/Problem.h"


struct solSD: public solution{ 
  std::vector<int> sol; 
  std::vector<bool> firstM; 

};


class SD: public Problem<solSD>{
	private:
		int nType;
		int numberJobs;
		int numberTools;
		int capacityMagazine;
		int capacityMin;
		int capacityMax;
		int type = 2;
		int sum = 0;
		bool debug = false;
		std::string fn;
		std::vector<std::vector<bool>>  toolJob;
		std::vector<std::vector<int>>  JobTools;
		std::vector<std::vector<int>>  changeTax;
		std::vector<int>  iniTax;
		thread_local static std::vector<std::pair <int,int>> twoSwap;
		
		void validateAndPrint(solSD sol, std::ofstream* file);
		double hungarianSolve(std::vector<std::vector<double>> costMatrix);
		std::vector<bool> hungarian2(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID, int *value);
		std::vector<int> hungarian3(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID, int lsml, int *value);
		double hungarianSolvePrint(const std::vector<std::vector<double>>& costMatrix,const std::vector<int>& rowsID,const std::vector<int>& colsID, std::ofstream* file);
	public:
		SD(std::string filename, int neighborType, int stype);
		~SD();
		solSD construction();
		solSD neighbor(solSD sol);
		solSD neighbor2opt(solSD sol);
		solSD neighbor2Swap(solSD sol);
		solSD neighborRandIns(solSD sol);
		solSD evaluate(solSD sol);
		solSD evaluate1(solSD sol);
		solSD evaluate2(solSD sol);
		solSD evaluate3(solSD sol);
		void setDebug(bool value);
		void print(solSD sol, std::string filename, float tempMin, float tempMax, int tempL, int MCL, int PTL, int time, int indexPTL, int indexEnd);
		
};

#endif 
