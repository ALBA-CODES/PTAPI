#ifndef __NODESWAP_H__
#define __NODESWAP_H__


#include <limits.h>
#include "NodeTempUp.h"
#include "ReplicaEst.h"

template<typename S>
class NodeSwap: public Node{
	private:
		Node* left; 
		Node* right;
		Node* tempUp;
		Consumer<S>* pool;
		int accept = 0;
		int upTime = INT_MAX;
		ReplicaEst* replicaEstL; 
		ReplicaEst* replicaEstR; 
		std::default_random_engine gen;
		std::uniform_real_distribution<double> dis;
	public:
		NodeSwap(Node* left_, Node* right_,Consumer<S>* pool_, atomic<int>* PTL_);
		NodeSwap(Node* left_, Node* right_, Node* tempUp_ ,Consumer<S>* pool_, atomic<int>* PTL_, int UPL_);
		~NodeSwap();
		void swap();
		void run();
		void start();
		bool ready();
		bool notify();
		bool observer(Node* from);
};

template<typename S>
NodeSwap<S>::NodeSwap(Node* left_, Node* right_, Consumer<S> *pool_, atomic<int>* PTL_)
:left(left_)
,right(right_)
,pool(pool_)
{
	execMax = PTL_;
	indexPT = pool_->getIndexPT();
	replicaEstL = new ReplicaEst();
	replicaEstR = new ReplicaEst();
}

template<typename S>
NodeSwap<S>::NodeSwap(Node* left_, Node* right_,Node* tempUp_ , Consumer<S>* pool_, atomic<int>* PTL_, int UPL_)
:left(left_)
,right(right_)
,tempUp(tempUp_)
,pool(pool_)
,upTime(UPL_)
{
	execMax = PTL_;
	indexPT = pool_->getIndexPT();
	replicaEstL = new ReplicaEst();
	replicaEstR = new ReplicaEst();
}

template<typename S>
NodeSwap<S>::~NodeSwap(){
}

template<typename S>
void NodeSwap<S>::run(){
	// Checks whether the execution has reached the end
	if(theEnd()){
		 pool->theEnd();
		endN = true;
	}

	double temp1 = ((NodeMCMC<S>*)left)->getTemp();
	double temp2 = ((NodeMCMC<S>*)right)->getTemp();
	double energia1 = ((NodeMCMC<S>*)left)->getEnergia();
	double energia2 = ((NodeMCMC<S>*)right)->getEnergia();
	replicaEstL->Push(energia1);
	replicaEstR->Push(energia2);
	
	// Calculation of the exponent of exp
	double deltaBeta = ((1.0/temp1) - (1.0/temp2)) * (energia1-energia2);
	//When delta is greater than and equal to 0, it always changes
	if(deltaBeta >= 0){									
		swap();							
		++accept;
	}else{
		//Calc exchange probability
		double probab;
		probab = exp(deltaBeta);
		//Troca aceita
		if(dis(gen) <= probab){
			swap();	
			++accept;
		} 
	} // End if/else
			
	// Reset to receve new request
	reset();
	
	if((execAtual%upTime) == 0){
		((NodeTempUp<S>*)tempUp)->setacceptRate((double)accept/(double)execAtual);
		((NodeTempUp<S>*)tempUp)->setEnergy(replicaEstL->Mean(), replicaEstR->Mean());
		tempUp->run();
		accept = 0;
		replicaEstL->Clear();
		replicaEstR->Clear();		
	}
		 
	// Notify nodes forward
	notify();
}

template<typename S>
bool NodeSwap<S>::ready(){
	
	for(vector<std::pair<Node*,bool>>::iterator it = edgeto.begin(); it != edgeto.end(); it++){
		if(!it->second) return false;
	}
	
	return true;
}

template<typename S>
bool NodeSwap<S>::notify(){
	
	for(vector<std::pair<Node*,bool>>::iterator it = edgeFrom.begin(); it != edgeFrom.end(); it++){
		Node* to = it->first;
		if(!to->observer(this))return false;	
	}
	
	return true;
}

template<typename S>
bool NodeSwap<S>::observer(Node* from){
	bool status = false;
	bool ready = true;
	
	mtxNode.lock();
	for(vector<std::pair<Node*,bool>>::iterator it = edgeto.begin(); it != edgeto.end(); it++){
		
		if(it->first == from){
			it->second = true;
			status = true;
		}
		
		if(!it->second) ready = false;
	}
	mtxNode.unlock();

	if((ready) && (!finish() || !endN)){ 
			pool->execAsync(this);

	}
	
	return status;
}

template<typename S>
void NodeSwap<S>::swap(){
	S aux = ((NodeMCMC<S>*)left)->getSol();
	((NodeMCMC<S>*)left)->setSol(((NodeMCMC<S>*)right)->getSol());
	((NodeMCMC<S>*)right)->setSol(aux); 
	
	// try to set the labels up and down
	((NodeMCMC<S>*)left)->trySetLabels();
	((NodeMCMC<S>*)right)->trySetLabels();
	
	// update flow
	((NodeMCMC<S>*)left)->updateFlow();
	((NodeMCMC<S>*)right)->updateFlow();
	
}

template<typename S>
void NodeSwap<S>::start(){
	pool->execAsync(this);
}

#endif
