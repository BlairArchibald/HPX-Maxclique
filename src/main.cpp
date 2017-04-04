#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <map>
#include <chrono>
#include <memory>
#include <typeinfo>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/lcos/broadcast.hpp>
#include <hpx/include/iostreams.hpp>
#include <hpx/include/serialization.hpp>
#include <hpx/include/thread_executors.hpp>

#include "incumbent_component.hpp"
#include "workqueue_component.hpp"

#include "DimacsParser.hpp"
#include "BitGraph.hpp"
#include "BitSet.hpp"

// 64 bit words
// 8 words covers 400 vertices
// 13 words covers 800 vertices
// Later we can specialise this at compile time
#define NWORDS 13

// Forward action decls
namespace graph {
  hpx::promise<int>* foundPromise;
  template<unsigned n_words_>
  auto maxcliqueTask(std::uint64_t spawndepth,
                     hpx::naming::id_type incumbent,
                     std::vector<unsigned> c,
                     BitSet<n_words_> p,
                     hpx::naming::id_type promise) -> void;

  template<unsigned n_words_>
  auto findcliqueTask(std::uint64_t spawnDepth,
                      const BitGraph<n_words_> graph,
                      const hpx::naming::id_type incumbent,
                      const hpx::naming::id_type found,
                      std::vector<unsigned> c,
                      BitSet<n_words_> p,
                      const hpx::naming::id_type promise) -> void;

  std::atomic<int> globalBound(0);
  BitGraph<NWORDS> globalGraph;
  auto updateBound(int newBound) -> void;
  auto broadcastBound(int newBound) -> void;

  template<unsigned n_words_>
  auto setGlobalGraph(const BitGraph<n_words_> graph) -> void;
}
HPX_PLAIN_ACTION(graph::maxcliqueTask<NWORDS>,  maxcliqueTaskAction)
HPX_PLAIN_ACTION(graph::setGlobalGraph<NWORDS>, setGlobalGraphAction)
HPX_PLAIN_ACTION(graph::findcliqueTask<NWORDS>, findcliqueTaskAction)
HPX_PLAIN_ACTION(graph::updateBound, updateBoundAction)

namespace scheduler {
  std::atomic<bool> running(true);
  hpx::lcos::local::counting_semaphore tasks_required_sem;
  hpx::naming::id_type local_workqueue;
  auto cancelScheduler() -> void;
  auto scheduler(std::vector<hpx::naming::id_type> workqueues) -> void;
}
HPX_PLAIN_ACTION(scheduler::cancelScheduler, cancelSchedulerAction)
HPX_PLAIN_ACTION(scheduler::scheduler, schedulerAction)

namespace graph {
  // Order a graphFromFile and return an ordered graph alongside a map to invert
  // the vertex numbering at the end.
  template<unsigned n_words_>
  auto orderGraphFromFile(const dimacs::GraphFromFile & g, std::map<int,int> & inv) -> BitGraph<n_words_> {
    std::vector<int> order(g.first);
    std::iota(order.begin(), order.end(), 0);

    // Order by degree, tie break on number
    std::vector<int> degrees;
    std::transform(order.begin(), order.end(), std::back_inserter(degrees),
                    [&] (int v) { return g.second.find(v)->second.size(); });

    std::sort(order.begin(), order.end(),
            [&] (int a, int b) { return ! (degrees[a] < degrees[b] || (degrees[a] == degrees[b] && a > b)); });


    // Construct a new graph with this new ordering
    BitGraph<n_words_> graph;
    graph.resize(g.first);

    for (unsigned i = 0 ; i < g.first ; ++i)
      for (unsigned j = 0 ; j < g.first ; ++j)
        if (g.second.find(order[i])->second.count(order[j]))
          graph.add_edge(i, j);

    // Create inv map (maybe just return order?)
    for (int i = 0; i < order.size(); i++) {
      inv[i] = order[i];
    }

    return graph;
    }

  template<unsigned n_words_>
  auto colour_class_order(const BitGraph<n_words_> & graph,
                          const BitSet<n_words_> & p,
                          std::array<unsigned, n_words_ * bits_per_word> & p_order,
                          std::array<unsigned, n_words_ * bits_per_word> & p_bounds) -> void {
    BitSet<n_words_> p_left = p; // not coloured yet
    unsigned colour = 0;         // current colour
    unsigned i = 0;              // position in p_bounds

    // while we've things left to colour
    while (! p_left.empty()) {
      // next colour
      ++colour;
      // things that can still be given this colour
      BitSet<n_words_> q = p_left;

      // while we can still give something this colour
      while (! q.empty()) {
        // first thing we can colour
        int v = q.first_set_bit();
        p_left.unset(v);
        q.unset(v);

        // can't give anything adjacent to this the same colour
        graph.intersect_with_row_complement(v, q);

        // record in result
        p_bounds[i] = colour;
        p_order[i] = v;
        ++i;
      }
    }
  }

  template <unsigned n_words_>
  auto maximiseExpand(const BitGraph<n_words_> & graph,
                      const hpx::naming::id_type incumbent,
                      std::vector<unsigned> & c,
                      BitSet<n_words_> & p) -> void {
    // initial colouring
    std::array<unsigned, n_words_ * bits_per_word> p_order;
    std::array<unsigned, n_words_ * bits_per_word> p_bounds;
    colour_class_order(graph, p, p_order, p_bounds);

    // for each v in p... (v comes later)
    for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
      auto bnd = globalBound.load();
      if (c.size() + p_bounds[n] < bnd)
        return;

      auto v = p_order[n];

      // consider taking v
      c.push_back(v);

      // filter p to contain vertices adjacent to v
      BitSet<n_words_> new_p = p;
      graph.intersect_with_row(v, new_p);

      if (c.size() > bnd) {
        std::set<int> members;
        for (auto & v : c) {
          members.insert(v);
        }
        hpx::lcos::broadcast_apply<updateBoundAction>(hpx::find_all_localities(), c.size());
        hpx::async<globalBound::incumbent::updateBound_action>(incumbent, c.size(), members).get();
      }

      maximiseExpand(graph, incumbent, c, new_p);

      // now consider not taking v
      c.pop_back();
      p.unset(v);
    }
  }

  template <unsigned n_words_>
  auto findExpand(const BitGraph<n_words_> & graph,
              const hpx::naming::id_type incumbent,
              const hpx::naming::id_type found,
              std::vector<unsigned> & c,
              BitSet<n_words_> & p) -> void {
    // initial colouring
    std::array<unsigned, n_words_ * bits_per_word> p_order;
    std::array<unsigned, n_words_ * bits_per_word> p_bounds;
    colour_class_order(graph, p, p_order, p_bounds);

    // for each v in p... (v comes later)
    for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
      auto bnd = globalBound.load();
      if (c.size() + p_bounds[n] < bnd)
        return;

      auto v = p_order[n];

      // consider taking v
      c.push_back(v);

      // filter p to contain vertices adjacent to v
      BitSet<n_words_> new_p = p;
      graph.intersect_with_row(v, new_p);

      if (c.size() == bnd) {
          std::set<int> members;
          for (auto & v : c) {
            members.insert(v);
          }
          hpx::async<globalBound::incumbent::updateBound_action>(incumbent, bnd, members).get();

          try {
            hpx::async<hpx::lcos::base_lco_with_value<int>::set_value_action>(found, 1).get();
          } catch (hpx::exception const& e){
            // Failed on double write, this is fine since we will terminate very soon anyway.
          }
      }

      findExpand(graph, incumbent, found, c, new_p);

      // now consider not taking v
      c.pop_back();
      p.unset(v);
    }
  }

  // Max
  template <unsigned n_words_>
  auto maximiseExpandSpawn(const BitGraph<n_words_> & graph,
                           const hpx::naming::id_type incumbent,
                           std::uint64_t spawnDepth,
                           std::vector<unsigned> & c,
                           BitSet<n_words_> & p) -> void {
    if (spawnDepth == 0) {
      maximiseExpand(graph, incumbent, c, p);
    } else {
      // initial colouring
      std::array<unsigned, n_words_ * bits_per_word> p_order;
      std::array<unsigned, n_words_ * bits_per_word> p_bounds;
      colour_class_order(graph, p, p_order, p_bounds);

      std::vector<hpx::future<int>> futures;

      auto bnd = globalBound.load();
      for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
        if (c.size() + p_bounds[n] < bnd)
          continue; // Don't spawn tasks that can't find the bound

        auto v = p_order[n];

        // consider taking v
        c.push_back(v);

        // filter p to contain vertices adjacent to v
        BitSet<n_words_> new_p = p;
        graph.intersect_with_row(v, new_p);

        if (c.size() > bnd) {
          std::set<int> members;
          for (auto & v : c) {
            members.insert(v);
          }
          hpx::lcos::broadcast_apply<updateBoundAction>(hpx::find_all_localities(), c.size());
          hpx::async<globalBound::incumbent::updateBound_action>(incumbent, c.size(), members).get();
        }

        auto promise = std::shared_ptr<hpx::promise<int>>(new hpx::promise<int>());
        auto f = promise->get_future();
        auto promise_id = promise->get_id();
        futures.push_back(std::move(f));

        hpx::util::function<void(hpx::naming::id_type)> task = hpx::util::bind(maxcliqueTaskAction(), _1, spawnDepth - 1, incumbent, c, new_p, promise_id);
        hpx::apply<workstealing::workqueue::addWork_action>(scheduler::local_workqueue, task);

        // now consider not taking v
        c.pop_back();
        p.unset(v);
      }
      // This task is going to sleep until all children are done, tell the scheduler to get more work
      scheduler::tasks_required_sem.signal();
      hpx::wait_all(futures);
    }
  }

  // Find
  template <unsigned n_words_>
  auto findExpandSpawn(const BitGraph<n_words_> & graph,
                   const hpx::naming::id_type incumbent,
                   const hpx::naming::id_type found,
                   std::uint64_t spawnDepth,
                   std::vector<unsigned> & c,
                   BitSet<n_words_> & p) -> void {
    if (spawnDepth == 0) {
      findExpand(graph, incumbent, found, c, p);
    } else {
      std::array<unsigned, n_words_ * bits_per_word> p_order;
      std::array<unsigned, n_words_ * bits_per_word> p_bounds;
      colour_class_order(graph, p, p_order, p_bounds);

      std::vector<hpx::future<int>> futures;

      auto bnd = globalBound.load();
      for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
        if (c.size() + p_bounds[n] < bnd)
          continue; // Don't spawn tasks that can't find the bound

        auto v = p_order[n];

        // consider taking v
        c.push_back(v);

        // filter p to contain vertices adjacent to v
        BitSet<n_words_> new_p = p;
        graph.intersect_with_row(v, new_p);

        if (c.size() == bnd) {
          std::set<int> members;
          for (auto & v : c) {
            members.insert(v);
          }
          hpx::async<globalBound::incumbent::updateBound_action>(incumbent, bnd, members).get();
          hpx::async<hpx::lcos::base_lco_with_value<int>::set_value_action>(found, 1).get();
        }

        auto promise = std::shared_ptr<hpx::promise<int>>(new hpx::promise<int>());
        auto f = promise->get_future();
        auto promise_id = promise->get_id();
        futures.push_back(std::move(f));

        hpx::util::function<void(hpx::naming::id_type)> task = hpx::util::bind(findcliqueTaskAction(), _1, spawnDepth - 1, graph, incumbent, found, c, new_p, promise_id);
        hpx::apply<workstealing::workqueue::addWork_action>(scheduler::local_workqueue, task);

        // now consider not taking v
        c.pop_back();
        p.unset(v);
      }
      // This task is going to sleep until all children are done, tell the scheduler to get more work
      scheduler::tasks_required_sem.signal();
      hpx::wait_all(futures);
    }
  }

  template<unsigned n_words_>
  void maxcliqueTask(std::uint64_t spawndepth, hpx::naming::id_type incumbent, std::vector<unsigned> c, BitSet<n_words_> p, hpx::naming::id_type promise) {
    maximiseExpandSpawn(globalGraph, incumbent, spawndepth, c, p);
    hpx::apply<hpx::lcos::base_lco_with_value<int>::set_value_action>(promise, 1);
    scheduler::tasks_required_sem.signal();
  }

  template<unsigned n_words_>
  void findcliqueTask(std::uint64_t spawnDepth,
                         const BitGraph<n_words_> graph,
                         const hpx::naming::id_type incumbent,
                         const hpx::naming::id_type found,
                         std::vector<unsigned> c,
                         BitSet<n_words_> p,
                         const hpx::naming::id_type promise) {
    findExpandSpawn(graph, incumbent, found, spawnDepth, c, p);
    hpx::apply<hpx::lcos::base_lco_with_value<int>::set_value_action>(promise, 1);
    scheduler::tasks_required_sem.signal();
  }

  auto updateBound(int newBound) -> void {
    while(true) {
      auto curBnd = globalBound.load();
      if (newBound < curBnd) {
        break;
      }

      if (globalBound.compare_exchange_weak(curBnd, newBound)) {
        break;
      }
    }
  }

  template<unsigned n_words_>
  auto setGlobalGraph(const BitGraph<n_words_> graph) -> void {
    globalGraph = graph;
  }

  template<unsigned n_words_>
  auto runMaxClique(const BitGraph<n_words_> & graph, std::uint64_t spawnDepth, const hpx::naming::id_type incumbent) -> void {
    BitSet<n_words_> p;
    p.resize(graph.size());
    p.set_all();

    std::vector<unsigned> c;
    maximiseExpandSpawn(graph, incumbent, spawnDepth, c, p);
  }

  template<unsigned n_words_>
  auto doSearch(const BitGraph<n_words_> & graph, std::uint64_t spawnDepth, hpx::naming::id_type incumbent, hpx::naming::id_type found) -> void {
    std::vector<unsigned> c;
    BitSet<n_words_> p; // Need to initialise this
    p.resize(graph.size());
    p.set_all();

    findExpandSpawn(graph, incumbent, found, spawnDepth, c, p);

    try {
      hpx::async<hpx::lcos::base_lco_with_value<int>::set_value_action>(found, 1).get();
    } catch (hpx::exception const& e){
      // Failed on double write, this is fine since we will terminate very soon anyway.
    }
  }

  template<unsigned n_words_>
  auto runFindClique(const BitGraph<n_words_> & graph, std::uint64_t spawnDepth, hpx::naming::id_type incumbent, std::uint64_t searchsize) -> void {
    BitSet<n_words_> p;
    p.resize(graph.size());
    p.set_all();

    auto done = hpx::lcos::broadcast<updateBoundAction>(hpx::find_all_localities(), searchsize);
    hpx::wait_all(done);

    foundPromise = new hpx::promise<int>();
    auto foundF  = foundPromise->get_future();
    auto foundId = foundPromise->get_id();

    hpx::apply(hpx::util::bind(&doSearch<n_words_>, graph, spawnDepth, incumbent, foundId));

    foundF.get();
  }
}

namespace scheduler {
  auto cancelScheduler() -> void {
    auto cur = running.load();
    while (!running.compare_exchange_weak(cur, false)) {
      ;
    }
  }

  auto scheduler(std::vector<hpx::naming::id_type> workqueues) -> void {
    auto here = hpx::find_here();
    hpx::naming::id_type last_remote = here; // We use *here* as a sentinel value
    auto distributed = hpx::find_all_localities().size() > 1;

    // Figure out which workqueue is local to this scheduler
    for (auto it = workqueues.begin(); it != workqueues.end(); ++it) {
      if (hpx::get_colocation_id(*it).get() == here) {
        local_workqueue = *it;
        workqueues.erase(it);
        break;
      }
    }

    auto threads = hpx::get_os_thread_count() == 1 ? 1 : hpx::get_os_thread_count() - 1;
    hpx::threads::executors::current_executor scheduler;

    // Pre-init the sem
    tasks_required_sem.signal(threads);

    // Debugging
    std::cout << "Running with: " << threads << " scheduler threads" << std::endl;

    while (running) {
      tasks_required_sem.wait();

      // Try local queue first then distributed
      hpx::util::function<void(hpx::naming::id_type)> task;
      task = hpx::async<workstealing::workqueue::steal_action>(local_workqueue).get();
      if (distributed && !task) {
        if (last_remote != here) {
          task = hpx::async<workstealing::workqueue::steal_action>(last_remote).get();
        }

        if (!task) {
          // Possibly a biased random function
          auto victim = workqueues.begin();
          std::advance(victim, std::rand() % workqueues.size());
          task = hpx::async<workstealing::workqueue::steal_action>(*victim).get();
          if (task) {
            last_remote = *victim;
          } else {
            last_remote = here;
          }
        }
      }

      if (task) {
        scheduler.add(hpx::util::bind(task, here));
      } else {
        hpx::this_thread::suspend(10);
        tasks_required_sem.signal();
      }
    }
  }
}

int hpx_main(boost::program_options::variables_map & opts) {
  auto inputFile = opts["input-file"].as<std::string>();
  if (inputFile.empty()) {
    hpx::finalize();
    return EXIT_FAILURE;
  }

  auto gFile = dimacs::read_dimacs(inputFile);

  // Order the graph (keep a hold of the map)
  std::map<int, int> invMap;
  auto graph = graph::orderGraphFromFile<NWORDS>(gFile, invMap);
  hpx::lcos::broadcast<setGlobalGraphAction>(hpx::find_all_localities(), graph).get();

  auto incumbent = hpx::new_<globalBound::incumbent>(hpx::find_here()).get();
  // Launch one workqueue per node
  auto localities = hpx::find_all_localities();
  std::vector<hpx::naming::id_type> workqueues;
  for (auto const& loc : localities) {
    workqueues.push_back(hpx::new_<workstealing::workqueue>(loc).get());
  }
  hpx::lcos::broadcast_apply<schedulerAction>(hpx::find_all_localities(), workqueues);
  hpx::this_thread::suspend(200); // Fix for now - give time for local_workqueue to be set

  auto spawnDepth = opts["spawn-depth"].as<std::uint64_t>();
  auto findClique = opts.count("find-clique");

  auto start_time = std::chrono::steady_clock::now();

  if (findClique) {
    auto expectedSize = opts["size-required"].as<std::uint64_t>();
    graph::runFindClique(graph, spawnDepth, incumbent, expectedSize);
  } else {
    graph::runMaxClique(graph, spawnDepth, incumbent);

  }
  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time);

  // Output result
  hpx::cout << "Size: " << hpx::async<globalBound::incumbent::getBound_action>(incumbent).get() << hpx::endl;

  auto members = hpx::async<globalBound::incumbent::getMembers_action>(incumbent).get();
  std::vector<unsigned> clique;
  for (auto const& m : members) {
    clique.push_back(invMap[m]);
  }
  std::sort(clique.begin(), clique.end());

  hpx::cout << "Members: " << hpx::endl;
  for (auto const & c : clique) {
    hpx::cout << c << " ";
  }
  hpx::cout << hpx::endl;

  hpx::cout << "cpu = " << overall_time.count() << std::endl;

  hpx::lcos::broadcast<cancelSchedulerAction>(hpx::find_all_localities()).get();

  if (findClique) {
    //We are done, just kill everything (not a particularly nice way to end).
    hpx::terminate();
    return 0;
  } else {
    // Terminate when everything is finished
    return hpx::finalize();
  }
}

int main (int argc, char* argv[]) {
  boost::program_options::options_description
    desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

  desc_commandline.add_options()
    ( "spawn-depth,d",
      boost::program_options::value<std::uint64_t>()->default_value(0),
      "Depth in the tree to spawn at"
    )
    ( "size-required,s",
      boost::program_options::value<std::uint64_t>()->default_value(0),
      "Size of clear we are searching for"
    )
    ( "find-clique",
      "Find a known clique rather than running the maximisation version"
    )
    ( "input-file,f",
      boost::program_options::value<std::string>(),
      "DIMACS formatted input graph"
    );

  return hpx::init(desc_commandline, argc, argv);
}


