#ifndef __NODE_H__
#define __NODE_H__

#include <vector>
#include <atomic>

// Dataflow node. A node becomes ready when all the nodes it depends on have
// finished their current execution. v2 replaces the mutex + flag vector of v1
// with a single atomic dependency counter.
class Node{
	private:
		std::vector<Node*> succ;       // nodes that depend on this one
		std::atomic<int> pending{0};   // dependencies not yet satisfied
		int nDeps = 0;
	public:
		Node(){}
		virtual ~Node(){}
		Node(const Node&) = delete;
		Node& operator=(const Node&) = delete;

		// Executes the node. Returns false when the node has finished for good
		// (nothing is notified and the graph drains naturally).
		virtual bool run() = 0;

		// Edge from -> to: 'to' can only run after 'from'.
		static void addEdge(Node* from, Node* to){
			from->succ.push_back(to);
			to->pending.store(++to->nDeps, std::memory_order_relaxed);
		}

		const std::vector<Node*>& successors() const { return succ; }
		void rearm(){ pending.store(nDeps, std::memory_order_release); }
		bool dependencyDone(){ return pending.fetch_sub(1, std::memory_order_acq_rel) == 1; }
};

#endif
