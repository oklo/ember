#include "../src/parallel_evaluate.hpp"
#include <atomic>
#include <cstdio>
#include <future>
#include <set>
#include <string>

int main() {
  using ember::detail::independent_evaluations;
  int failures=0,checks=0;
  auto check=[&](bool ok,const char* label) {
    ++checks;if(!ok){++failures;std::printf("FAIL %s\n",label);}
  };
  constexpr std::size_t count=160;
  std::array<std::thread::id,count> first{},second{};
  std::array<unsigned,count> visits{};
  auto evaluate=[&](auto& ids) {
    independent_evaluations(count,4,[&](std::size_t i) {
      ids[i]=std::this_thread::get_id();++visits[i];
    });
  };
  evaluate(first);evaluate(second);
  check(first==second,"workers persist with stable index ownership");
  check(std::set(first.begin(),first.end()).size()==4,"four active threads");
  check(std::all_of(visits.begin(),visits.end(),[](auto n){return n==2;}),"each index evaluated once per batch");
  // Decreasing the requested count leaves unused workers asleep; increasing
  // it again must reuse them without counting an idle worker in the barrier.
  independent_evaluations(count,2,[&](std::size_t i){second[i]=std::this_thread::get_id();});
  check(std::set(second.begin(),second.end()).size()==2,"requested two threads honored");
  evaluate(second);check(first==second,"workers retained across thread-count changes");

  std::array<bool,count> nested_ok{};
  independent_evaluations(count,4,[&](std::size_t i) {
    const auto caller=std::this_thread::get_id();bool ok=true;
    independent_evaluations(128,4,[&](std::size_t) {ok=ok && std::this_thread::get_id()==caller;});
    nested_ok[i]=ok;
  });
  check(std::all_of(nested_ok.begin(),nested_ok.end(),[](bool v){return v;}),"nested evaluation stays on its calling thread");

  std::atomic<unsigned> completed{};std::string error;
  try {
    independent_evaluations(count,4,[&](std::size_t i) {
      ++completed;if(i==3 || i==90)throw std::runtime_error(std::to_string(i));
    });
  }catch(const std::exception& e){error=e.what();}
  check(error=="3","first failing index independent of worker order");
  check(completed==count,"all workers joined before error is returned");
  evaluate(second);check(first==second,"worker group usable after callback failure");
  for(auto invalid:{0u,65u}) {
    bool rejected=false;
    try {independent_evaluations(count,invalid,[](std::size_t){});}
    catch(const std::invalid_argument&){rejected=true;}
    check(rejected,"invalid thread count rejected");
  }
  unsigned empty_calls=0;independent_evaluations(0,4,[&](std::size_t){++empty_calls;});
  check(empty_calls==0,"empty batch returns without callbacks");
  auto concurrent=[] {
    std::array<unsigned,count> values{};
    for(unsigned repeat=0;repeat<20;++repeat)
      independent_evaluations(count,4,[&](std::size_t i){values[i]+=static_cast<unsigned>(i)+1;});
    for(std::size_t i=0;i<count;++i)if(values[i]!=20*(i+1))return false;
    return true;
  };
  auto a=std::async(std::launch::async,concurrent);
  auto b=std::async(std::launch::async,concurrent);
  check(a.get() && b.get(),"concurrent callers own separate reusable worker groups");
  std::printf("parallel evaluation: %d checks, %d failures\n",checks,failures);
  return failures?1:0;
}
