#ifndef __PT_H__
#define __PT_H__


#include <vector>
#include <deque> 
#include <thread>         
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <tuple>
#include <math.h>
#include <algorithm>
#include <bits/stdc++.h>
#include "Problem.h"
#include "Consumer.h"
#include "NodeMCMC.h"
#include "NodeSwap.h"

template<typename S>
class PT{
	private:
		// PT varibles
		float tempMin_; // init temp
		float tempMax_; // end temp
		float perPTLEnd = 0.0; // PT end PTL best sol
		int tempL_; // temp size
		int MKL_; // MC length
		int PTL_; // number of swaps
		int ut; // update type
		int tempUpdate; // number of swaps to update de temp
		int indexBestPTL; // Best PTL solution
		int timePTEnd = 0; // Time to finish the PT execution
		std::deque<double> allTemps;
		S best;
	public:
		PT(float tempMin, float tempMax, int tempL, int MKL, int PTL, int tempD, int upType, int tempUp); // constructor

		~PT();
		S start(int thN, Problem<S>* prob);
		S getBestSol();
		int getBestPTLSolIndex();
		void setTimeEnd(int timeEnd);
		void setPTLStopEnd(float per);
		std::deque<double> tempPG(float tempMin, float tempMax, int tempL); 
		std::deque<double> tempExp(float tempMin, float tempMax, int tempL); 
		std::deque<double> tempLinear(float tempMin, float tempMax, int tempL); 
		std::deque<double> tempILinear(float tempMin, float tempMax, int tempL); 
};


// inicialização do PT. valores da temp. tamanho da MCMC e o Número de propostas de troca entre temp
template<typename S>
PT<S>::PT(float tempMin, float tempMax, int tempL, int MKL, int PTL, int tempD, int upType, int tempUp)
	: tempMin_(tempMin)
	, tempMax_(tempMax)
	, tempL_(tempL)
	, MKL_(MKL)
	, PTL_(PTL)
	, ut(upType)
	, tempUpdate(tempUp)
{
	switch (tempD){
		case 1:{
			allTemps = tempLinear(tempMin,tempMax,tempL);
		}
		break;
		case 2:{
			allTemps = tempILinear(tempMin,tempMax,tempL);
		}
		break;
		case 3:{
			allTemps = tempExp(tempMin,tempMax,tempL);
		}
		break;
		case 4:{
			allTemps = tempPG(tempMin,tempMax,tempL);
		}
		break;
		default:
			cout << "Invalid temperature initial distribution! \n";
		break;	
	}

} 

template<typename S>
PT<S>::~PT(){

}

template<typename S>
S PT<S>::start(int thN, Problem<S>* prob){	
	Consumer<S>* consumer = new Consumer<S>(thN);
	atomic<int> PTLEnd = PTL_;


	//variaveis
	Node* nMCMC;
	Node* nSwap;
	Node* nUpTemp;	
	Node* nMCMCAux;
	Node* nSwapAux;
	Node* nUpTempAux;
		    
    // Creates the first MCMC node
	nMCMC = new NodeMCMC<S>(MKL_,&PTLEnd,allTemps.front(), prob, consumer, timePTEnd, perPTLEnd);
	((NodeMCMC<S>*)nMCMC)->setFirstTemp(); // check First temp 
	consumer->setMaxEnd();
	allTemps.pop_front();
	nMCMCAux = nMCMC; 
		
    // Creates the second MCMC node
	nMCMC = new NodeMCMC<S>(MKL_,&PTLEnd,allTemps.front(), prob, consumer, timePTEnd, perPTLEnd);
	consumer->setMaxEnd();
	allTemps.pop_front();
	
	// Creates the first exchange node
	if(!ut) nSwap = new NodeSwap<S>(nMCMC,nMCMCAux, consumer, &PTLEnd);
	else{ 
		nUpTemp = new NodeTempUp<S>(nMCMC,nMCMCAux, NULL, consumer,ut, tempL_);
		nSwap = new NodeSwap<S>(nMCMC,nMCMCAux,nUpTemp, consumer, &PTLEnd, tempUpdate);	
	}
	
	consumer->setMaxEnd();


	// edges arriving in swap node
	nSwap->addEdge(nMCMCAux,nSwap);
	nSwap->addEdge(nMCMC,nSwap);
	
	// edges leaving the swap node
	nSwap->addEdge(nSwap, nMCMCAux);
	nSwap->addEdge(nSwap, nMCMC);
	
	// updates the edges
	nMCMCAux->addEdge(nMCMCAux,nSwap);
	nMCMC->addEdge(nMCMC,nSwap);
	nMCMCAux->addEdge(nSwap,nMCMCAux);
	nMCMC->addEdge(nSwap,nMCMC);

	
	// execution begins	
	consumer->execAsync(nMCMCAux);

	// Retrieve the node on the right
	nMCMCAux = nMCMC;
	nSwapAux = nSwap;
	nUpTempAux = nUpTemp;
	
	// Create the remaining nodes	
	while(!allTemps.empty()){

		nMCMC = new NodeMCMC<S>(MKL_,&PTLEnd, allTemps.front(),prob, consumer, timePTEnd, perPTLEnd);
		consumer->setMaxEnd();
		allTemps.pop_front();

		// Creates the exchange node
		if(!ut) nSwap = new NodeSwap<S>(nMCMC,nMCMCAux,consumer,&PTLEnd);
		else {
			nUpTemp = new NodeTempUp<S>(nMCMC,nMCMCAux, nUpTempAux, consumer,ut, tempL_);
			((NodeTempUp<S>*)nUpTempAux)->setUpTempUp(nUpTemp);
			nSwap = new NodeSwap<S>(nMCMC,nMCMCAux,nUpTemp,consumer, &PTLEnd, tempUpdate);
		}
		
		consumer->setMaxEnd();
		
		// Updates incoming edges
		nSwap->addEdge(nMCMCAux,nSwap);
		nSwap->addEdge(nMCMC,nSwap);
		nSwap->addEdge(nSwapAux,nSwap);
		nMCMCAux->addEdge(nSwap, nMCMCAux);
		nMCMC->addEdge(nSwap,nMCMC);		

		// Updates the edges come out 
		nMCMCAux->addEdge(nMCMCAux,nSwap);
		nMCMC->addEdge(nMCMC,nSwap);			
		nSwapAux->addEdge(nSwapAux,nSwap);
		nSwap->addEdge(nSwap, nMCMCAux);
		nSwap->addEdge(nSwap, nMCMC);

		// start exec
		consumer->execAsync(nMCMCAux);
		
		// update varibles
		nMCMCAux = nMCMC;
		nSwapAux = nSwap;
		nUpTempAux = nUpTemp;
	}

	//check last temp
	((NodeMCMC<S>*)nMCMCAux)->setLastTemp();
	if(ut)((NodeTempUp<S>*)nUpTempAux)->setUpTempUp(NULL);
	//start exec
	consumer->execAsync(nMCMCAux);


consumer->finished();
best = consumer->getBestSol();
std::atomic<int>*  indexPT = consumer->getIndexPT(); 
indexBestPTL = (*indexPT).load(); 

return best;

}

template<typename S>
std::deque<double> PT<S>::tempPG(float tempMin, float tempMax, int tempL){
	// variables
	std::deque<double> temps(tempL);
	double b = tempMax/tempMin;
	double ex = 1.0/((double)(tempL - 1));	
	double r =  std::pow(b,ex); 
	// Calc temps
	for(unsigned int ind = 0; ind < temps.size();++ind){
		temps[ind] = tempMin * std::pow(r,ind);	
	}

	return temps;
}

template<typename S>
std::deque<double> PT<S>::tempExp(float tempMin, float tempMax, int tempL){
	// variables
	std::deque<double>temps(tempL);
	double tempMaxExp = log(tempMax);
	double tempMinExp = log(tempMin);

	// Calc temps
	for(unsigned int ind = 0; ind < temps.size();++ind){
		double currentTemp = tempMinExp + ((tempMaxExp - tempMinExp) * ((double)ind/(double)(tempL - 1)));
		temps[ind] = exp(currentTemp);
	}
	return temps;
} 

template<typename S>
std::deque<double> PT<S>::tempLinear(float tempMin, float tempMax, int tempL){
	// variables
	std::deque<double> temps(tempL);
	
	// Calc temps
	for(unsigned int ind = 0; ind < temps.size();++ind){
		double currentTemp = tempMin + ((tempMax - tempMin) * ((double)ind/(double)(tempL - 1)));
		temps[ind] = currentTemp;
	}
	
	return temps;
}

template<typename S>
std::deque<double> PT<S>::tempILinear(float tempMin, float tempMax, int tempL){
	// variables
	std::deque<double> temps(tempL);
	double bm = 1.0/tempMax;
	double b1 = 1.0/tempMin;
	
	// Calc temps
	for(unsigned int ind = 0; ind < temps.size();++ind){
		// Calc temps
		double bi = bm + ((b1-bm) * ((double)ind/(double)(tempL - 1)));
		temps[ind] = 1.0 / bi;
	}
	return temps;
}

template<typename S>
S PT<S>::getBestSol(){
 return best;
}

template<typename S>
int PT<S>::getBestPTLSolIndex(){
 return indexBestPTL;
}

template<typename S>
void PT<S>::setTimeEnd(int timeEnd){
 timePTEnd = timeEnd;
}

template<typename S>
void PT<S>::setPTLStopEnd(float per){
 perPTLEnd = per;
}


#endif
