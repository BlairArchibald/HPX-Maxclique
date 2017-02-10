#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <map>
#include <chrono>

#include <hpx/hpx.hpp>
#include <hpx/hpx_init.hpp>
#include <hpx/include/lcos.hpp>
#include <hpx/include/iostreams.hpp>

#include "incumbent_component.hpp"

#include "DimacsParser.hpp"
#include "BitGraph.hpp"
#include "BitSet.hpp"

// 64 bit words
// 8 words covers 400 vertices
// Later we can specialise this at compile time
#define NWORDS 8

namespace graph {

  static unsigned long long colourings = 0;

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

    ++colourings;
    // initial colouring
    std::array<unsigned, n_words_ * bits_per_word> p_order;
    std::array<unsigned, n_words_ * bits_per_word> p_bounds;
    colour_class_order(graph, p, p_order, p_bounds);

    // for each v in p... (v comes later)
    for (int n = p.popcount() - 1 ; n >= 0 ; --n) {
      // bound, timeout or early exit?

      auto bnd = hpx::async<globalBound::incumbent::getBound_action>(incumbent).get();
      if (c.size() + p_bounds[n] <= bnd)
        return;

      auto v = p_order[n];

      // consider taking v
      c.push_back(v);

      // filter p to contain vertices adjacent to v
      BitSet<n_words_> new_p = p;
      graph.intersect_with_row(v, new_p);

      if (new_p.empty()) {
        if (c.size() > bnd) {
          std::set<int> members;
          for (auto & v : c) {
            members.insert(v);
          }
          // Fire and forget updates
          hpx::apply<globalBound::incumbent::updateBound_action>(incumbent, c.size(), members);
        }
      }
      else
        expand(graph, incumbent, c, new_p);

      // now consider not taking v
      c.pop_back();
      p.unset(v);
    }
  }

  template<unsigned n_words_>
  auto runMaxClique(const BitGraph<n_words_> & graph, hpx::naming::id_type incumbent) -> void {
    std::vector<unsigned> c;
    c.reserve(graph.size());

    BitSet<n_words_> p;
    p.resize(graph.size());
    p.set_all();

    // Updates happen to incumbent as a side effect
    expand(graph, incumbent, c, p);
  }

}


int hpx_main(int argc, char* argv[]) {
  if (2 != argc) {
    std::cout << "Usage: " << argv[0] << " file" << std::endl;
    hpx::finalize();
    return EXIT_FAILURE;
  }

  auto gFile = dimacs::read_dimacs(std::string(argv[1]));

  // Order the graph (keep a hold of the map)
  std::map<int, int> invMap;
  auto graph = graph::orderGraphFromFile<NWORDS>(gFile, invMap);

  // Run Maxclique
  auto incumbent = hpx::new_<globalBound::incumbent>(hpx::find_here()).get();

  auto start_time = std::chrono::steady_clock::now();

  graph::runMaxClique(graph, incumbent);

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
  hpx::cout << "colourings = " << graph::colourings << std::endl;

  hpx::finalize();
}

int main (int argc, char* argv[]) {
  return hpx::init(argc, argv);
}


