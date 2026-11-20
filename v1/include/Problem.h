#ifndef __PROBLEM_H__
#define __PROBLEM_H__

#include <vector>
#include <climits>


struct solution{ 
  double evalSol = INT_MAX; 
  bool Nup,Ndown;
};

template<typename S>
class Problem{
	private:
		solution bestSol;
	public:
		Problem();
		~Problem();
		void setBestSol(S sol);
		virtual S construction()=0;
		virtual S neighbor(S sol)=0;
		virtual S evaluate(S sol)=0;
};

template<typename S>
Problem<S>::Problem(){
}

template<typename S>
Problem<S>::~Problem(){
}

template<typename S>
void Problem<S>::setBestSol(S sol){
}


#endif
