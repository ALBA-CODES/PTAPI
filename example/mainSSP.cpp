/*
 * main.cpp
 * 
 * Copyright 2021 André Luis <André Luis@DESKTOP-HDL2CBS> 
 *
 * Compilation Instructions 
 * g++ mainSSP.cpp ../include/*.h ./SSP.cpp -std=c++20 -Wshadow -Wall -o mainSSP -O3 -Wno-unused-result -lpthread
 *
 * Example Execution
 * ./mainSSP ../instances/SSP/F1001.txt --TEMP_INIT 0.2 --TEMP_FIM 1 --N_REPLICAS 16 --MCL 400 --PTL 1000 --TEMP_DIST 1 --TYPE_UPDATE 2 --TEMP_UPDATE 35000
 *
 */
#include <cstdlib>
#include "../include/ExecTime.h"
#include "./SSP.h"
#include "../include/PT.h"


using namespace std;

int main(int argc, char* argv[])
{
	//varibles
	float tempIni = 0.01;
	float tempfim = 2.0;
	int tempN = 10;
	int MCL = 10;
	int PTL = 10;	
	int tempUp = 50;
	int tempD = 1;
	int uType = 0;
	int thN = thread::hardware_concurrency();	
	vector<std::string> arguments(argv + 1, argv + argc);	
	
	// Instance file name
	std::string fn = arguments[0];
	
	// Read arguments
	for(unsigned int i=1; i<arguments.size(); i+=2)
	{		            
        if(arguments[i]== "--TEMP_FIM")
            tempfim = std::stof(arguments[i+1]);
        else if(arguments[i]== "--TEMP_INIT")
            tempIni = std::stof(arguments[i+1]);
        else if(arguments[i]== "--N_REPLICAS")
            tempN = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--MCL")
            MCL  = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--PTL")
            PTL = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TEMP_DIST")
            tempD = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TYPE_UPDATE")
            uType = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--TEMP_UPDATE")
            tempUp = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--THREAD_USED")
            thN = std::stoi(arguments[i+1]);
    }
	
	// Create SSP object
	SSP* prob = new SSP(fn);
	
	// Create and start PT 

	PT<solSSP> algo(tempIni,tempfim,tempN,MCL,PTL,tempD,uType,tempUp);
	// algo.setPTLStopEnd(0.1); Parar o PT se não encontrar uma solução melhor após (0.1 * PTL)  
	// algo.setTimeEnd(50); Parar por tempo em segundos
	ExecTime et;
	solSSP sol = algo.start(thN, prob);
	int index = algo.getBestPTLSolIndex();
	et.printTime();
	std::cout<<"Best sol value:"<<sol.evalSol<<" PTL index:"<<index<<"\n";

	return 0;
}

