/* 
Copyright (c) 2012 Jakob Progsch, Václav Zeman

This software is provided 'as-is', without any express or implied
warranty. In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

   1. The origin of this software must not be misrepresented; you must not
   claim that you wrote the original software. If you use this software
   in a product, an acknowledgment in the product documentation would be
   appreciated but is not required.

   2. Altered source versions must be plainly marked as such, and must not be
   misrepresented as being the original software.

   3. This notice may not be removed or altered from any source
   distribution.

Notice of Alteration
William Silversmith
May 2019, December 2023, January 2025

- The license file was moved from a seperate file to the top of this one.
- Created public "join" member function from destructor code.
- Created public "start" member function from constructor code.
- Used std::invoke_result_t to update to modern C++
- Added thread indicator as worker function parameter and
  removed futures since I'm not using them and it was difficult
  to combine that logic with a thread indicator.
*/

#ifndef THREAD_POOL_H
#define THREAD_POOL_H

#include <vector>
#include <queue>
#include <memory>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>

class ThreadPool {
public:
    ThreadPool(size_t);
    void enqueue(std::function<void(size_t)>); 
    void start(size_t);
    void join();
    ~ThreadPool();
private:
    // need to keep track of threads so we can join them
    std::vector< std::thread > workers;
    // the task queue
    std::queue< std::function<void(size_t)> > tasks;
    
    // synchronization
    std::mutex queue_mutex;
    std::condition_variable condition;
    bool stop;
};
 
// the constructor just launches some amount of workers
inline ThreadPool::ThreadPool(size_t threads)
    :   stop(false)
{
    start(threads);
}

void ThreadPool::start(size_t threads) {
    stop = false;
    for(size_t i = 0;i<threads;++i)
        workers.emplace_back(
            [this,i]
            {
                for(;;)
                {
                    std::function<void(size_t)> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock,
                            [this]{ return this->stop || !this->tasks.empty(); });
                        if(this->stop && this->tasks.empty())
                            return;
                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }
                    
                    task(i);
                }
            }
        );
}

// add new work item to the pool
void ThreadPool::enqueue(std::function<void(size_t)> f) {
    auto task = std::make_shared< std::packaged_task<void(size_t)> >(f);
    
    {
        std::unique_lock<std::mutex> lock(queue_mutex);

        // don't allow enqueueing after stopping the pool
        if(stop)
            throw std::runtime_error("enqueue on stopped ThreadPool");

        tasks.emplace([task](size_t t){ (*task)(t); });
    }
    condition.notify_one();
}

inline void ThreadPool::join () {
    {
        std::unique_lock<std::mutex> lock(queue_mutex);
        stop = true;
    }
    condition.notify_all();
    for(std::thread &worker: workers)
        worker.join();

    workers.clear();
}

// the destructor joins all threads
inline ThreadPool::~ThreadPool() {
    join();
}



#endif