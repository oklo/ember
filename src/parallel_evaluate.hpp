#pragma once
#include <algorithm>
#include <array>
#include <condition_variable>
#include <exception>
#include <functional>
#include <mutex>
#include <stdexcept>
#include <thread>
#include <vector>

namespace ember::detail {
// Each calling thread owns its workers. They retain exact physics caches
// between evaluations; no model or callback survives completion of a batch.
class EvaluationWorkers {
 public:
  EvaluationWorkers() { workers_.reserve(63); }
  EvaluationWorkers(const EvaluationWorkers&) = delete;
  EvaluationWorkers& operator=(const EvaluationWorkers&) = delete;
  ~EvaluationWorkers() {
    { std::lock_guard lock(mutex_); stopping_=true; }
    ready_.notify_all();
    for(auto& thread:workers_)thread.join();
  }

  void run(std::size_t count,const std::function<void(std::size_t)>& work) {
    std::unique_lock lock(mutex_);
    while(workers_.size()+1<count) {
      const auto id=workers_.size()+1,seen=generation_;
      workers_.emplace_back([this,id,seen] { worker(id,seen); });
    }
    work_=work;
    active_=count;
    remaining_=count-1;
    errors_.fill({});
    ++generation_;
    lock.unlock();
    ready_.notify_all();
    try {work(0);}catch(...) {errors_[0]=std::current_exception();}
    lock.lock();
    done_.wait(lock,[&] {return remaining_==0;});
    work_={};
    for(const auto& error:errors_)if(error)std::rethrow_exception(error);
  }

 private:
  void worker(std::size_t id,std::size_t seen) {
    std::unique_lock lock(mutex_);
    for(;;) {
      ready_.wait(lock,[&] {return stopping_ || generation_!=seen;});
      if(stopping_)return;
      seen=generation_;
      if(id>=active_)continue;
      lock.unlock();
      try {work_(id);}catch(...) {errors_[id]=std::current_exception();}
      lock.lock();
      if(--remaining_==0)done_.notify_one();
    }
  }
  std::mutex mutex_;
  std::condition_variable ready_,done_;
  std::vector<std::thread> workers_;
  std::function<void(std::size_t)> work_;
  std::array<std::exception_ptr,64> errors_{};
  std::size_t generation_{},active_{},remaining_{};
  bool stopping_{};
};

inline EvaluationWorkers& evaluation_workers() {
  thread_local EvaluationWorkers workers;
  return workers;
}
inline thread_local bool inside_parallel_evaluation=false;
struct EvaluationScope {
  EvaluationScope() {inside_parallel_evaluation=true;}
  ~EvaluationScope() {inside_parallel_evaluation=false;}
};

template<class Function>
void independent_evaluations(std::size_t count,std::size_t threads,const Function& evaluate) {
  if(threads==0 || threads>64)
    throw std::invalid_argument("evaluation thread count must be between 1 and 64");
  const auto workers=std::min(threads,(count+31)/32);
  // A callback can perform another independent evaluation. Keep that inner
  // call serial rather than waiting recursively on the same worker group.
  if(workers<=1 || inside_parallel_evaluation) {
    for(std::size_t i=0;i<count;++i)evaluate(i);
    return;
  }
  std::vector<std::exception_ptr> failures(count);
  evaluation_workers().run(workers,[&](std::size_t worker) {
    EvaluationScope scope;
    // Fixed interleaved chunks balance inner/outer stellar layers and keep
    // an unchanged zone on the same thread for subsequent exact-cache hits.
    for(auto begin=8*worker;begin<count;begin+=8*workers)
      for(auto i=begin;i<std::min(begin+8,count);++i)try {evaluate(i);}
      catch(...) {failures[i]=std::current_exception();}
  });
  for(const auto& error:failures)if(error)std::rethrow_exception(error);
}
} // namespace ember::detail
