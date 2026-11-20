#ifndef __NODETEMPUP_H__
#define __NODETEMPUP_H__

#include <random>
#include <math.h>
#include "Pool.h"
#include "NodeMCMC.h"

template<typename S>
class NodeTempUp: public Node{
	private:
		Node* left; 
		Node* right;
		Node* downTempUp;
		Node* upTempUp;
		Consumer<S>* pool;
		int m = 0;
		int ut = 0;
		double c = 0;
		double fl = 0;
		double fr = 0;
		double deltaE = 0;
		double energyL = 0;
		double energyR = 0;
		double acceptRate = 0;
		std::vector<Node *> allTemps;
	public:
		void run();
		void start();
		bool ready();
		bool notify();
		~NodeTempUp();
		double getC();
		void tempUpFB();
		void tempUp23();
		double getDeltaE();
		void tempUpEnergy();
		double getDownTemp();
		double getProbDens();
		double getLeftTemp();
		bool observer(Node* from);
		void setUpTempUp(Node* utu);
		void setacceptRate(double aR);
		void setEnergy(double l, double r);
		NodeTempUp(Node* left_, Node* right_,Node* downTempUp_, Consumer<S>* pool_, int ut_, int tL);
};

template<typename S>
NodeTempUp<S>::NodeTempUp(Node* left_, Node* right_,Node* downTempUp_, Consumer<S> *pool_, int ut_, int tL)
:left(left_)
,right(right_)
,downTempUp(downTempUp_)
,pool(pool_)
,m(tL)
,ut(ut_)
{}

template<typename S>
NodeTempUp<S>::~NodeTempUp(){
}

template<typename S>
void NodeTempUp<S>::run(){
    switch (ut){
		case 1:{
			tempUp23();
		}
		break;
		case 2:{
			tempUpEnergy();
		}
		break;
		case 3:{
			tempUpFB();
		}
		break;
		default:
			cout << "Invalid temp update type! \n";
		break;
	}
}

template<typename S>
void NodeTempUp<S>::tempUp23(){

    if(acceptRate < 0.2) {

		double k = (energyR - energyL)/((NodeMCMC<S>*)right)->getTemp();
		double newT = (energyL - energyR) / (log(0.23) - k);
		if(newT > 0.0)((NodeMCMC<S>*)left)->setTemp(newT);
	}
}

template<typename S>
void NodeTempUp<S>::tempUpEnergy(){
	deltaE = energyL - energyR;
	
	if(downTempUp){
		double k1 = ((NodeTempUp<S>*)downTempUp)->getDeltaE();
		double tMenos = ((NodeTempUp<S>*)downTempUp)->getDownTemp();
		double eq = ((k1/tMenos) + (deltaE/((NodeMCMC<S>*)left)->getTemp()));
		double newT = ((deltaE+k1)/eq);
		if(newT > 0.0)((NodeMCMC<S>*)right)->setTemp(newT);
	}
}

template<typename S>
void NodeTempUp<S>::tempUpFB(){
	
	if(!downTempUp){
		double density = this->getProbDens();
		double dfdt = ((fl - fr)/(((NodeMCMC<S>*)left)->getTemp() - ((NodeMCMC<S>*)right)->getTemp()));//(fr - fl)
		if((dfdt != 0) && (density != 0) && (!isnan(density))){
		c = 1.0/density;
		double eq = pow((1.0/(double)(m*c)),2)/abs(dfdt);
		double newT = eq + ((NodeMCMC<S>*)right)->getTemp();
		if((newT > ((NodeMCMC<S>*)right)->getTemp()) && ((newT < ((NodeTempUp<S>*)upTempUp)->getLeftTemp())))((NodeMCMC<S>*)left)->setTemp(newT);
		}
	}else if(upTempUp){
		c = ((NodeTempUp<S>*)downTempUp)->getC();
		double dfdt = ((fl - fr)/(((NodeMCMC<S>*)left)->getTemp() - ((NodeMCMC<S>*)right)->getTemp()));
		if((dfdt != 0) && (c != 0)){
		double eq = abs(pow((1.0/(double)(m*c)),2)/dfdt);
		double newT = eq + ((NodeMCMC<S>*)right)->getTemp();
		if((newT > ((NodeMCMC<S>*)right)->getTemp()) && ((newT < ((NodeTempUp<S>*)upTempUp)->getLeftTemp())))((NodeMCMC<S>*)left)->setTemp(newT);	
		}
	}
	
}

template<typename S>
bool NodeTempUp<S>::ready(){
	
	for(vector<std::pair<Node*,bool>>::iterator it = edgeto.begin(); it != edgeto.end(); it++){
		if(!it->second) return false;
	}
	
	return true;
}

template<typename S>
bool NodeTempUp<S>::notify(){
	
	for(vector<std::pair<Node*,bool>>::iterator it = edgeFrom.begin(); it != edgeFrom.end(); it++){
		Node* to = it->first;
		if(!to->observer(this))return false;	
	}
	
	return true;
}

template<typename S>
bool NodeTempUp<S>::observer(Node* from){
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

	if((ready) && (!finish())){ 
			pool->execAsync(this);

	}
	
	return status;
}

template<typename S>
void NodeTempUp<S>::start(){
	pool->execAsync(this);
}

template<typename S>
void NodeTempUp<S>::setacceptRate(double aR){
	acceptRate = aR;
}

template<typename S>
void NodeTempUp<S>::setEnergy(double l, double r){
 energyL = l;
 energyR = r;
}

template<typename S>
double NodeTempUp<S>::getDeltaE(){
 return deltaE;
}

template<typename S>
double NodeTempUp<S>::getDownTemp(){
	return ((NodeMCMC<S>*)right)->getTemp(); 
}

template<typename S>
void NodeTempUp<S>::setUpTempUp(Node* utu){
	upTempUp = utu;
}

template<typename S>
double NodeTempUp<S>::getProbDens(){
	double den;
	double atual;
	
	// Get flow
	fr = ((NodeMCMC<S>*)right)->getFlow();
	fl = ((NodeMCMC<S>*)left)->getFlow();
	//(fl > fr)
	if((fr == 0) || (fl==0)) return 0;

	((NodeMCMC<S>*)right)->resetFlow();
	double dfdt = abs(((fl - fr)/(((NodeMCMC<S>*)left)->getTemp() - ((NodeMCMC<S>*)right)->getTemp())));

	den = ((NodeMCMC<S>*)left)->getTemp() - ((NodeMCMC<S>*)right)->getTemp(); 
	//abs(fr - fl)
	atual = sqrt(dfdt/den)*den;
	if(upTempUp)atual += ((NodeTempUp<S>*)upTempUp)->getProbDens();
	return atual;
}

template<typename S>
double NodeTempUp<S>::getC(){
	return c;
}

template<typename S>
double NodeTempUp<S>::getLeftTemp(){
	return ((NodeMCMC<S>*)left)->getTemp();
}




#endif
