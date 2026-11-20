#ifndef __POOL_H__
#define __POOL_H__

#include <thread>		// std::this_thread::sleep_for
#include <chrono>		// std::chrono::seconds
#include <queue>		//std::queue
#include <deque>
#include <mutex>
#include <atomic>
#include <iostream> 	
#include <optional>
#include "Node.h"

using namespace std;

class Pool{
	private:
		queue<Node*> data;
		std::mutex mtx;
		atomic_int end = 0;
		int maxEnd = 0;
	public:
		Pool();
		Pool(Pool&& obj);
		Pool& operator=(const Pool& obj);
		~Pool();
		size_t size();
		void push(Node* item);
		Node* pop();
};

Pool::Pool(){}
Pool::Pool(Pool&& obj){}
Pool& Pool::operator=(const Pool& obj){
 return *this;
}
Pool::~Pool(){}
		
size_t Pool::size(){
	const std::lock_guard<std::mutex> lock(mtx);
	size_t result = data.size();
	return result;
}


void Pool::push(Node* item){
	const std::lock_guard<std::mutex> lock(mtx);
	data.push(item);
}

Node* Pool::pop()
{
	const std::lock_guard<std::mutex> lock(mtx);

		if (data.size() <= 0)
		{
			return NULL;
		}
		Node* item = std::move(data.front()); data.pop();
		return item;
}


#endif
