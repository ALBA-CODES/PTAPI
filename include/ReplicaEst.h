#ifndef __REPLICAEST_H__
#define __REPLICAEST_H__

// Running mean of the energy of a replica between two temperature updates.
// (v1 also accumulated the 2nd..4th moments on every exchange, never used.)
class ReplicaEst{
	private:
		long long n = 0;
		double sum = 0.0;
	public:
		void Clear(){ n = 0; sum = 0.0; }
		void Push(double x){ ++n; sum += x; }
		long long NumDataValues() const { return n; }
		double Mean() const { return n ? sum/(double)n : 0.0; }
};

#endif
