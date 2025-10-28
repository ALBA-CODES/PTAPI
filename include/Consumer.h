#ifndef __CONSUMER_H__
#define __CONSUMER_H__

#include <memory>
#include <thread>  
#include <vector> 
#include <tuple>      
#include <iostream>
#include <sstream>
#include "Pool.h"

using namespace std;

template<typename S>
class Consumer{
	private:
		unsigned threadCount;
		vector<thread> threads;
		vector<Pool> buff;
		atomic<unsigned> index{0};
		atomic<int> bestV{INT_MAX};
		atomic<int> indexPT{0};
		atomic<int> maxDif{0};
		bool finish{true};
		mutex mtxSol;
		atomic_int end{0};
		int maxEnd{0};
		S bestSol;
	public:		
		Consumer(unsigned nTh);
		void execAsync(Node* item);
		void run(unsigned i);
		void finished();
		bool theEnd();
		bool theEnd(S sol_);
		void setMaxEnd(int maxEnd_);
		void setBestSol(S b, unsigned i);
		void setMaxEnd();
		S getBestSol();
		int getEnd();
		atomic<int>* getIndexPT();
		atomic<int>* getMaxDif();
		int getStopC();
};



template<typename S>
Consumer<S>::Consumer(unsigned nTh){
	
	threadCount = nTh;
	
	for (unsigned n = 0; n != threadCount; ++n) {
		buff.emplace_back();
	}
	
	for (unsigned n = 0; n != threadCount; ++n) {
		threads.emplace_back([&, n]{ run(n); });
	}
}

template<typename S>
void Consumer<S>::run(unsigned i){
	
	Node* c;
	while(finish){
		for(unsigned n = 0; n != threadCount; ++n){
			if((c = buff[(i + n) % threadCount].pop())){
				c->run();
				break;
			}
		}
	}

}

template<typename S>
void Consumer<S>::finished(){
// End the thread	
	while(!(maxEnd <= end)) finish = true;	
	finish = false;
	for(auto& t:threads) t.join();
}

template<typename S>
bool Consumer<S>::theEnd(){
	end++;
	return (maxEnd <= end);	
}

template<typename S>
bool Consumer<S>::theEnd(S sol_){
 
	{
		unique_lock<mutex> lock{mtxSol};
		end++;
		if (sol_.evalSol < bestSol.evalSol) bestSol=sol_;
	}
return (maxEnd <= end);	
}

template<typename S>
void Consumer<S>::setBestSol(S b, unsigned i){
	if(b.evalSol < bestV){
		bestV = b.evalSol;
		indexPT = i;
	}
}

template<typename S>
void Consumer<S>::setMaxEnd(int maxEnd_){
	maxEnd = maxEnd_;	
}

template<typename S>
void Consumer<S>::setMaxEnd(){
	++maxEnd;	
}

template<typename S>
S Consumer<S>::getBestSol(){
 return bestSol;	
}

template<typename S>
int Consumer<S>::getEnd(){
 return end;	
}

template<typename S>
atomic<int>* Consumer<S>::getIndexPT(){
 return &indexPT;
}

template<typename S>
atomic<int>* Consumer<S>::getMaxDif(){
 return &maxDif;
}

template<typename S>
void Consumer<S>::execAsync(Node* item){
	auto i = index++;
	buff[i % threadCount].push(item);
}

#endif
