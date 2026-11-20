/*
 * main.cpp
 * 
 * Copyright 2021 André Luis <André Luis@DESKTOP-HDL2CBS>
 *
 * Compilation Instructions (PTAPI v2, C++17 or newer)
 * g++ mainIPMTC.cpp IPMTC.cpp -std=c++17 -O3 -pthread -o mainIPMTC
 *
 * Example Execution
 * ./mainIPMTC ../instances/IPMTC/instanceLarge1_m\=3_n\=50_l\=30_c\=10_s\=0.PMTC --TEMP_INIT 0.2 --TEMP_FIM 1 --N_REPLICAS 16 --MCL 400 --PTL 1000 --TEMP_DIST 1 --TYPE_UPDATE 2 --TEMP_UPDATE 35000
 *
 * Options added in v2: --SEED n (reproducible run), --TIME_LIMIT s (seconds, may be fractional),
 * --CALIBRATE 1|2|3 (temperatures chosen by the API), --SWAP_SCHEME 0|1 (0 sequential, 1 even-odd)
 *
 */
#include <cstdlib>
#include "../include/ExecTime.h"
#include "./IPMTC.h"
#include "../include/PT.h"


using namespace std;

int main(int argc, char* argv[])
{
	//varibles
	float tempIni = 0.01;
	float tempfim = 2.0;
	int tempN = 10;
	int MCL = 0;
	int PTL = 2;	
	int tempUp = 50;
	int tempD = 1;
	int uType = 0;
	int thN = thread::hardware_concurrency();
	long long seed = -1;
	int calibrate = 0;
	double timeLimit = 0.0;
	int swapScheme = 0;	
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
        else if(arguments[i]== "--SEED") // v2: reproducible run
            seed = std::stoll(arguments[i+1]);
        else if(arguments[i]== "--TIME_LIMIT") // v2: seconds
            timeLimit = std::stod(arguments[i+1]);
        else if(arguments[i]== "--SWAP_SCHEME") // v2: 0 = sequential, 1 = even-odd
            swapScheme = std::stoi(arguments[i+1]);
        else if(arguments[i]== "--CALIBRATE") // v2: temperatures chosen by the API: 1 = first and last, 2 = ladder for 23% of exchanges, 3 = both
            calibrate = std::stoi(arguments[i+1]);
    }
	
	// Create SSP object
	IPMTC* prob = new IPMTC(fn);
	
	// Create and start PT 
	PT<solIPMTC> algo(tempIni,tempfim,tempN,MCL,PTL,tempD,uType,tempUp);
	if(seed >= 0) algo.setSeed((unsigned long long)seed);
	if(timeLimit > 0.0) algo.setTimeLimit(timeLimit);
	algo.setSwapScheme(swapScheme);
	if(calibrate & 1){
		std::pair<double,double> t = algo.calibrateTemperatures(prob);
		std::cout<<"Calibrated temperatures: "<<t.first<<" "<<t.second<<"\n";
	}
	if(calibrate & 2){
		std::vector<double> l = algo.calibrateLadder(prob, thN);
		std::cout<<"Calibrated ladder: "<<l.front()<<" "<<l.back()<<"\n";
	}
	ExecTime et;
	solIPMTC sol = algo.start(thN, prob);
	// v2: statistics of the run
	std::cout<<"PTL executed:"<<algo.getPTLEnd()<<" evaluations:"<<algo.getEvaluations()<<" time to best (ms):"<<algo.getTimeToBestMs()
	         <<" threads:"<<algo.getThreadsUsed()<<" seed:"<<algo.getSeed()<<" round trips:"<<algo.getRoundTrips()<<"\n";
	// Divides between machines
	sol = prob->setCut(sol); 
	et.printTime();
	std::cout<<sol.evalSol<<"\n";

	return 0;
}

