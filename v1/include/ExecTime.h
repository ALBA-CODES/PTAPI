#ifndef __EXECTIME_H__
#define __EXECTIME_H__

#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

class ExecTime{
	private:
		std::chrono::high_resolution_clock::time_point start;

	public:
		ExecTime();
		void printTime();
		int getTimeMs();
		~ExecTime();
};


ExecTime::ExecTime()
	: start(std::chrono::high_resolution_clock::now()) 
	{}

void ExecTime::printTime(){	
	auto duration = std::chrono::high_resolution_clock::now() - start;
	std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(duration).count() << " ms\n";	
} 

int ExecTime::getTimeMs(){
	auto duration = std::chrono::high_resolution_clock::now() - start;
	return std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();	
}

ExecTime::~ExecTime(){
}

#endif
