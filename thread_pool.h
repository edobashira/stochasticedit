#include <mutex>
#include <thread>
#include <condition_variable>

/* A semaphore, from http://p9as.blogspot.com/2012/06/c11-semaphores.html */

class Semaphore {
  int value_;
  int wakeups_;
  std::mutex mutex_;
  std::condition_variable cond_;
  
 public:
  Semaphore(int value) : value_(value), wakeups_(0) {}
  
  void wait();
  void signal();

};

void Semaphore::signal() {
  std::lock_guard<std::mutex> lock(mutex_);
  value_++;
  if (value_ <= 0) {
      wakeups_++;
      cond_.notify_one();
    }
}

void Semaphore::wait() {
  std::unique_lock<std::mutex> lock(mutex_);
  value_--;
  if (value_ < 0) {
      cond_.wait(lock, [this] { return this->wakeups_ > 0; });
      wakeups_--;
    }
}

/* A simple "thread pool", which limits the amount of running threads to a fixed number */

class ThreadPool {
  unsigned n_threads;
  Semaphore semaphore;

 public:
  ThreadPool(unsigned n_threads) : semaphore(n_threads), n_threads(n_threads) {}

  template<typename F>
  void enqueue(F f) {
    semaphore.wait();
    std::thread([this, f] {
        f();
        this->semaphore.signal();
    }).detach();
  }

  void join() {
    for(unsigned k = 0; k < n_threads; k++)
        semaphore.wait();
  }
};
