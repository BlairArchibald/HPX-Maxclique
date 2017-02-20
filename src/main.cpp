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
#include <hpx/include/iostreams.hpp>
#include <hpx/include/serialization.hpp>

#include "incumbent_component.hpp"
#include "workqueue_component.hpp"

#include "DimacsParser.hpp"
#include "BitGraph.hpp"
#include "BitSet.hpp"

// 64 bit words
// 8 words covers 400 vertices
// Later we can specialise this at compile time
#define NWORDS 8

// Forward action decls
namespace graph {
  template<unsigned n_words_>
  auto maxcliqueTask(const BitGraph<n_words_> graph, hpx::naming::id_type incumbent, std::vector<unsigned> c, BitSet<n_words_> p, hpx::naming::id_type promise) -> void;

  std::atomic<int> globalBound(0);
  auto updateBound(int newBound) -> void;
  auto broadcastBound(int newBound) -> void;

}
HPX_PLAIN_ACTION(graph::maxcliqueTask<NWORDS>, maxcliqueTask400Action)
HPX_PLAIN_ACTION(graph::updateBound, updateBoundAction)
HPX_PLAIN_ACTION(graph::broadcastBound, broadcastBoundAction)

// For distributed promises
HPX_REGISTER_ACTION(hpx::lcos::base_lco_with_value<int>::set_value_action, set_value_action_int);

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
  auto expand(const BitGraph<n_words_> & graph,
              const hpx::naming::id_type & incumbent,
              std::vector<unsigned> & c,
              BitSet<n_words_> & p) -> void {
    // initial colouring
    std::array<unsigned, n_words_ * bits_per_word> p_order;
    std::array<unsigned, n_words_ * bits_per_word> p_bounds;
    colour_class_order(graph, p, p_order, p_bounds);

    // for each v in p... (v comes later)
    for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
      auto bnd = globalBound.load();
      if (c.size() + p_bounds[n] <= bnd)
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
          hpx::async<globalBound::incumbent::updateBound_action>(incumbent, c.size(), members).get();
          hpx::async<broadcastBoundAction>(hpx::find_here(), c.size()).get();
      }

      expand(graph, incumbent, c, new_p);

      // now consider not taking v
      c.pop_back();
      p.unset(v);
    }
  }

  template<unsigned n_words_>
  auto spawnTasks(const BitGraph<n_words_> & graph,
                  std::uint64_t spawnDepth,
                  std::vector<unsigned> c,
                  BitSet<n_words_> p,
                  hpx::naming::id_type workqueue,
                  hpx::naming::id_type incumbent,
                  std::vector<hpx::future<int>> & futures) -> void {
    std::array<unsigned, n_words_ * bits_per_word> p_order;
    std::array<unsigned, n_words_ * bits_per_word> p_bounds;
    colour_class_order(graph, p, p_order, p_bounds);

    for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
      auto v = p_order[n];
      c.push_back(v);

      // filter p to contain vertices adjacent to v
      BitSet<n_words_> new_p = p;
      graph.intersect_with_row(v, new_p);

      if (spawnDepth == 0) {
        auto promise = std::shared_ptr<hpx::promise<int>>(new hpx::promise<int>());
        auto f = promise->get_future();
        auto promise_id = promise->get_id();

        futures.push_back(std::move(f));

        hpx::util::function<void(hpx::naming::id_type)> task = hpx::util::bind(maxcliqueTask400Action(), _1, graph, incumbent, c, new_p, promise_id);
        hpx::apply<workstealing::workqueue::addWork_action>(workqueue, task);
      } else {
        spawnTasks(graph, spawnDepth - 1, c, new_p, workqueue, incumbent, futures);
      }

      c.pop_back();
      p.unset(v);
    }
  }

  template<unsigned n_words_>
  void maxcliqueTask(const BitGraph<n_words_> graph, hpx::naming::id_type incumbent, std::vector<unsigned> c, BitSet<n_words_> p, hpx::naming::id_type promise) {
    expand(graph, incumbent, c, p);
    hpx::apply<hpx::lcos::base_lco_with_value<int>::set_value_action>(promise, 1);
    return;
  }

  // Atomically update the global bound on all localities
  auto broadcastBound(int newBound) -> void {
    auto localities = hpx::find_all_localities();
    for (auto const & node : localities) {
      hpx::apply<updateBoundAction>(node, newBound);
    }
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
  auto runMaxClique(const BitGraph<n_words_> & graph, std::uint64_t spawnDepth, hpx::naming::id_type incumbent, hpx::naming::id_type workqueue) -> void {
    BitSet<n_words_> p;
    p.resize(graph.size());
    p.set_all();

    std::vector<hpx::future<int>> futures;
    std::vector<unsigned> c;
    spawnTasks(graph, spawnDepth, c, p, workqueue, incumbent, futures);
    hpx::wait_all(futures);
  }
}


void scheduler(hpx::naming::id_type workqueue) {
  auto threads = hpx::get_os_thread_count() == 1 ? 1 : hpx::get_os_thread_count() - 1;
  hpx::threads::executors::local_queue_executor scheduler(threads);

  // Debugging
  std::cout << "Running with: " << threads << " scheduler threads" << std::endl;

  bool running = true;
  while (running) {
    auto pending = scheduler.num_pending_closures();
    if (pending < threads) {
      auto task = hpx::async<workstealing::workqueue::steal_action>(workqueue).get();
      if (task) {
        auto t = hpx::util::bind(task, hpx::find_here());
        scheduler.add(t);
      }
    } else {
      hpx::this_thread::suspend();
    }
  }
}
HPX_PLAIN_ACTION(scheduler, schedulerAction)

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

  // Run Maxclique
  auto incumbent = hpx::new_<globalBound::incumbent>(hpx::find_here()).get();
  auto workqueue = hpx::new_<workstealing::workqueue>(hpx::find_here()).get();

  // Start a scheduler on each node
  auto localities = hpx::find_all_localities();
  for (auto const & node : localities) {
    // Can't use async here since the scheduler never returns
    hpx::apply<schedulerAction>(node, workqueue);
  }

  auto spawnDepth = opts["spawn-depth"].as<std::uint64_t>();
  auto start_time = std::chrono::steady_clock::now();
  graph::runMaxClique(graph, spawnDepth, incumbent, workqueue);
  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - start_time);

  // Output result
  hpx::cout << "Size: " << hpx::async<globalBound::incumbent::getBound_action>(incumbent).get() << hpx::endl;

  auto members = hpx::async<globalBound::incumbent::getMembers_action>(incumbent).get();
  hpx::cout << "Members: " << hpx::endl;
  for (auto const& m : members) {
    hpx::cout << invMap[m] + 1 << " ";
  }
  hpx::cout << hpx::endl << hpx::flush;

  hpx::cout << "cpu = " << overall_time.count() << std::endl;

  // TODO: A nicer termination
  //hpx::finalize();
  hpx::this_thread::suspend(2000);
  hpx::terminate();
}

int main (int argc, char* argv[]) {
  boost::program_options::options_description
    desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

  desc_commandline.add_options()
    ( "spawn-depth,d",
      boost::program_options::value<std::uint64_t>()->default_value(1),
      "Depth in the tree to spawn at"
    )
    ( "input-file,f",
      boost::program_options::value<std::string>(),
      "DIMACS formatted input graph"
    );

  return hpx::init(desc_commandline, argc, argv);
}


