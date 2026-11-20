#ifndef __POOL_H__
#define __POOL_H__

#include <deque>
#include <mutex>
#include <thread>
#include <vector>
#include <atomic>
#include <condition_variable>
#include "Node.h"

#if defined(__x86_64__) || defined(__i386__) || defined(_M_X64)
#include <immintrin.h>
#define PTAPI_CPU_RELAX() _mm_pause()
#else
#define PTAPI_CPU_RELAX() std::this_thread::yield()
#endif

// Number of polls of the queue an idle worker does before it goes to sleep.
#ifndef PTAPI_SPIN
#define PTAPI_SPIN 5000
#endif

// Thread pool that executes the dataflow graph.
//  - idle workers spin for a short time and then SLEEP (v1: endless busy wait);
//  - a worker that makes other nodes ready runs one of them itself
//    (no queue round trip) and publishes the others;
//  - termination is detected by counting the nodes in flight.
class Pool{
	private:
		std::deque<Node*> data;
		std::mutex mtx;
		std::condition_variable cvWork;
		std::condition_variable cvDone;
		std::vector<std::thread> threads;
		std::atomic<int> queued{0};
		std::atomic<int> inFlight{1};   // 1 = token held by the owner until seal()
		int sleepers = 0;
		bool quit = false;
		bool done = false;

		Node* tryPop(){
			if(queued.load(std::memory_order_relaxed) == 0) return nullptr;
			std::lock_guard<std::mutex> lock(mtx);
			if(data.empty()) return nullptr;
			Node* item = data.front();
			data.pop_front();
			queued.fetch_sub(1, std::memory_order_relaxed);
			return item;
		}

		Node* pop(){
			for(int spin = 0; spin < PTAPI_SPIN; ++spin){
				if(Node* item = tryPop()) return item;
				PTAPI_CPU_RELAX();
			}
			std::unique_lock<std::mutex> lock(mtx);
			++sleepers;
			cvWork.wait(lock, [this]{ return quit || !data.empty(); });
			--sleepers;
			if(data.empty()) return nullptr;
			Node* item = data.front();
			data.pop_front();
			queued.fetch_sub(1, std::memory_order_relaxed);
			return item;
		}

		void release(){
			if(inFlight.fetch_sub(1, std::memory_order_acq_rel) == 1){
				std::lock_guard<std::mutex> lock(mtx);
				done = true;
				cvDone.notify_all();
			}
		}

		void work(){
			while(Node* c = pop()){
				while(c){
					Node* next = nullptr;
					if(c->run()){
						c->rearm();
						for(Node* to : c->successors()){
							if(to->dependencyDone()){
								if(!next){ next = to; inFlight.fetch_add(1, std::memory_order_relaxed); }
								else push(to);
							}
						}
					}
					release();
					c = next;
				}
			}
		}

	public:
		explicit Pool(unsigned nTh){
			if(nTh == 0) nTh = 1;
			for(unsigned n = 0; n != nTh; ++n) threads.emplace_back([this]{ work(); });
		}

		~Pool(){
			{
				std::lock_guard<std::mutex> lock(mtx);
				quit = true;
			}
			cvWork.notify_all();
			for(auto& t : threads) t.join();
		}

		unsigned size() const { return (unsigned)threads.size(); }

		void push(Node* item){
			inFlight.fetch_add(1, std::memory_order_relaxed);
			bool wake;
			{
				std::lock_guard<std::mutex> lock(mtx);
				data.push_back(item);
				queued.fetch_add(1, std::memory_order_relaxed);
				wake = sleepers > 0;
			}
			if(wake) cvWork.notify_one();
		}

		// The owner has pushed every initial node.
		void seal(){ release(); }

		// Blocks until the graph has drained. If 'deadline' is reached first,
		// onTimeout() is called once (it must make the nodes stop).
		template<typename Clock, typename F>
		void wait(bool hasDeadline, std::chrono::time_point<Clock> deadline, F onTimeout){
			std::unique_lock<std::mutex> lock(mtx);
			if(hasDeadline && !cvDone.wait_until(lock, deadline, [this]{ return done; })){
				lock.unlock();
				onTimeout();
				lock.lock();
			}
			cvDone.wait(lock, [this]{ return done; });
		}

		// Blocks until the graph has drained or 'until' is reached. Returns true if drained.
		template<typename Clock>
		bool waitUntil(std::chrono::time_point<Clock> until){
			std::unique_lock<std::mutex> lock(mtx);
			return cvDone.wait_until(lock, until, [this]{ return done; });
		}
		void waitDone(){
			std::unique_lock<std::mutex> lock(mtx);
			cvDone.wait(lock, [this]{ return done; });
		}
};

#endif
