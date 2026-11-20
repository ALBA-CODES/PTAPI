#ifndef __EXECTIME_H__
#define __EXECTIME_H__

#include <chrono>
#include <iostream>
#include <string>
#include <fstream>

class ExecTime{
	private:
		std::chrono::steady_clock::time_point start;

	public:
		ExecTime(): start(std::chrono::steady_clock::now()) {}
		~ExecTime(){}

		void reset(){ start = std::chrono::steady_clock::now(); }

		void printTime(){
			std::cout << getTimeMs() << " ms\n";
		}

		int getTimeMs(){
			auto duration = std::chrono::steady_clock::now() - start;
			return (int)std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();
		}
};

#endif
