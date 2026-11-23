//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <gtest/gtest.h>
#include <random>
#include <algorithm>
#include <chrono>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <vector>
#include "fields2cover/types.h"
#include "fields2cover/objectives/rp_obj/direct_dist_path_obj.h"
#include "fields2cover/route_planning/route_planner_base.h"
#include "fields2cover.h"

TEST(fields2cover_rp_route_plan_base, simple_example) {
  f2c::Random rand(4);
  F2CCells cells {
    F2CCell(F2CLinearRing({
          F2CPoint(0,0), F2CPoint(2,0),F2CPoint(2,2),F2CPoint(0,2), F2CPoint(0,0)
    }))
  };
  cells.addRing(0, F2CLinearRing({
          F2CPoint(.4,.4), F2CPoint(.4,.6),F2CPoint(.6,.6),F2CPoint(.6,.4), F2CPoint(.4,.4)
  }));
  cells.addRing(0, F2CLinearRing({
          F2CPoint(1.2,1.2), F2CPoint(1.2,1.6),F2CPoint(1.6,1.6),F2CPoint(1.6,1.2), F2CPoint(1.2,1.2)
  }));

  cells *= 3e1;

  f2c::hg::ConstHL const_hl;
  F2CCells no_hl = const_hl.generateHeadlandArea(cells, 1, 3);
  auto hl_swaths = const_hl.generateHeadlandSwaths(cells, 1, 3, false);

  f2c::decomp::BoustrophedonDecomp decomp;
  decomp.setSplitAngle(M_PI/2.0);
  auto decomp_cells = decomp.decompose(no_hl);

  f2c::sg::BruteForce bf;
  F2CSwathsByCells swaths = bf.generateSwaths(M_PI/2.0, 5, decomp_cells);

  f2c::rp::RoutePlannerBase route_planner;
  F2CRoute route = route_planner.genRoute(hl_swaths[1], swaths);

  EXPECT_FALSE(route.isEmpty());
  EXPECT_GT(route.sizeVectorSwaths(), 1);
  EXPECT_EQ(route.sizeVectorSwaths(), route.sizeConnections());

  /*
  f2c::Visualizer::figure();
  f2c::Visualizer::plot(cells);
  //f2c::Visualizer::plot(no_hl);
  //f2c::Visualizer::plot(hl_swaths[1]);
  //f2c::Visualizer::plot(decomp_cells);
  //for (auto&& c : route.getConnections())
  f2c::Visualizer::plot(route);
  f2c::Visualizer::show();
  */
}

TEST(fields2cover_rp_route_plan_base, redirect_flag) {
  f2c::Random rand(4);
  F2CCells cells {
    F2CCell(F2CLinearRing({
          F2CPoint(0,0), F2CPoint(2,0),F2CPoint(2,2),F2CPoint(0,2), F2CPoint(0,0)
    }))
  };
  cells.addRing(0, F2CLinearRing({
          F2CPoint(.4,.4), F2CPoint(.4,.6),F2CPoint(.6,.6),F2CPoint(.6,.4), F2CPoint(.4,.4)
  }));
  cells.addRing(0, F2CLinearRing({
          F2CPoint(1.2,1.2), F2CPoint(1.2,1.6),F2CPoint(1.6,1.6),F2CPoint(1.6,1.2), F2CPoint(1.2,1.2)
  }));

  cells *= 3e1;

  f2c::hg::ConstHL const_hl;
  F2CCells no_hl = const_hl.generateHeadlandArea(cells, 1, 3);
  auto hl_swaths = const_hl.generateHeadlandSwaths(cells, 1, 3, false);

  f2c::decomp::BoustrophedonDecomp decomp;
  decomp.setSplitAngle(M_PI/2.0);
  auto decomp_cells = decomp.decompose(no_hl);

  f2c::sg::BruteForce bf;
  F2CSwathsByCells swaths = bf.generateSwaths(M_PI/2.0, 5, decomp_cells);

  f2c::rp::RoutePlannerBase route_planner;
  F2CRoute route = route_planner.genRoute(hl_swaths[1], swaths, false, 1e-4, false);
  
  F2CSwaths old_swaths = swaths.flatten();
  F2CSwaths new_swaths;
  for (size_t sbc = 0; sbc < route.sizeVectorSwaths(); ++sbc) {
    new_swaths.append(route.getSwaths(sbc));
  }

  EXPECT_EQ(new_swaths.size(), old_swaths.size());
  
  for (size_t s = 0; s < new_swaths.size(); ++s) {
    F2CSwath old_swath = old_swaths.at(s);
    F2CSwath new_swath = new_swaths.at(s);
    EXPECT_TRUE(new_swath.hasSameDir(old_swath));
  }

  EXPECT_FALSE(route.isEmpty());
  EXPECT_GT(route.sizeVectorSwaths(), 1);
  EXPECT_EQ(route.sizeVectorSwaths(), route.sizeConnections());

  /*
  f2c::Visualizer::figure();
  f2c::Visualizer::plot(cells);
  //f2c::Visualizer::plot(no_hl);
  //f2c::Visualizer::plot(hl_swaths[1]);
  //f2c::Visualizer::plot(decomp_cells);
  //for (auto&& c : route.getConnections())
  f2c::Visualizer::plot(route);
  f2c::Visualizer::show();
  */
}

// Cost of the closed tour depot -> route... -> depot.
int64_t tourCost(const F2CGraph2D& g, const std::vector<long long int>& route,
    size_t depot) {
  int64_t cost = 0;
  size_t prev = depot;
  for (auto id : route) {
    cost += g.getCostFromEdge(prev, static_cast<size_t>(id));
    prev = static_cast<size_t>(id);
  }
  return cost + g.getCostFromEdge(prev, depot);
}

// Graph2D numbers nodes by the order their points are first added, and
// numNodes() counts points. Register node i as point (i, 0) up front (a
// zero-cost self-loop, which the solver never uses) so index i is node i.
F2CGraph2D graphWithNodes(size_t n) {
  F2CGraph2D g;
  for (size_t i = 0; i < n; ++i) {
    g.addDirectedEdge(F2CPoint(i, 0), F2CPoint(i, 0), 0);
  }
  return g;
}

void addArc(F2CGraph2D& g, size_t from, size_t to, int64_t cost) {
  g.addDirectedEdge(F2CPoint(from, 0), F2CPoint(to, 0), cost);
}

// Complete directed graph with reproducible pseudo-random arc costs.
F2CGraph2D completeGraph(size_t n, unsigned int seed) {
  std::mt19937 gen(seed);
  std::uniform_int_distribution<int64_t> dist(10, 1000);
  F2CGraph2D g = graphWithNodes(n);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      if (i != j) {
        addArc(g, i, j, dist(gen));
      }
    }
  }
  return g;
}

TEST(fields2cover_rp_route_plan_base, best_route_with_params_is_optimal) {
  const size_t n = 7;
  const size_t depot = n - 1;
  F2CGraph2D g = completeGraph(n, 7);

  std::vector<long long int> nodes(n - 1);
  std::iota(nodes.begin(), nodes.end(), 0);
  int64_t optimum = std::numeric_limits<int64_t>::max();
  do {
    optimum = std::min(optimum, tourCost(g, nodes, depot));
  } while (std::next_permutation(nodes.begin(), nodes.end()));

  f2c::rp::RoutePlannerBase route_planner;
  auto route = route_planner.computeBestRouteWithParams(
      g, "time_limit { seconds: 5 }", true, 200);

  ASSERT_EQ(route.size(), n - 1);
  std::vector<long long int> sorted = route;
  std::sort(sorted.begin(), sorted.end());
  std::vector<long long int> expected(n - 1);
  std::iota(expected.begin(), expected.end(), 0);
  EXPECT_EQ(sorted, expected);
  EXPECT_EQ(tourCost(g, route, depot), optimum);
}

TEST(fields2cover_rp_route_plan_base, best_route_with_params_removes_arcs) {
  // Only one Hamiltonian cycle exists through the real arcs:
  // 4 -> 0 -> 1 -> 2 -> 3 -> 4. Cheap shortcuts do not exist as arcs, so a
  // route using them would mean a missing arc was taken.
  F2CGraph2D g = graphWithNodes(5);
  for (auto [from, to, cost] : std::vector<std::tuple<size_t, size_t, int64_t>>{
      {4, 0, 10}, {0, 1, 10}, {1, 2, 10}, {2, 3, 10}, {3, 4, 10},
      {0, 2, 1000}, {1, 3, 1000}, {2, 0, 1000}}) {
    addArc(g, from, to, cost);
  }

  f2c::rp::RoutePlannerBase route_planner;
  for (bool remove : {true, false}) {
    auto route = route_planner.computeBestRouteWithParams(
        g, "time_limit { seconds: 5 }", remove, 100);
    EXPECT_EQ(route, (std::vector<long long int>{0, 1, 2, 3}))
      << "remove_missing_arcs=" << remove;
  }
}

TEST(fields2cover_rp_route_plan_base,
    best_route_with_params_falls_back_when_pruning_is_infeasible) {
  // Node 2 has no outgoing arc, so no tour exists through real arcs only.
  // With pruning the model is infeasible; the function must solve again
  // with missing arcs kept at INF cost and still return every node.
  F2CGraph2D g = graphWithNodes(4);
  for (auto [from, to, cost] : std::vector<std::tuple<size_t, size_t, int64_t>>{
      {3, 0, 10}, {0, 1, 10}, {1, 2, 10}, {1, 3, 10}}) {
    addArc(g, from, to, cost);
  }

  f2c::rp::RoutePlannerBase route_planner;
  auto route = route_planner.computeBestRouteWithParams(
      g, "time_limit { seconds: 5 }", true, 100);
  std::sort(route.begin(), route.end());
  EXPECT_EQ(route, (std::vector<long long int>{0, 1, 2}));
}

TEST(fields2cover_rp_route_plan_base,
    best_route_with_params_stall_limit_is_deterministic_and_stops_early) {
  F2CGraph2D g = completeGraph(40, 11);
  f2c::rp::RoutePlannerBase route_planner;

  auto t0 = std::chrono::steady_clock::now();
  auto first = route_planner.computeBestRouteWithParams(
      g, "time_limit { seconds: 60 }", true, 500);
  const double first_s = std::chrono::duration<double>(
      std::chrono::steady_clock::now() - t0).count();
  auto second = route_planner.computeBestRouteWithParams(
      g, "time_limit { seconds: 60 }", true, 500);

  EXPECT_EQ(first.size(), 39u);
  EXPECT_EQ(first, second);
  // The stall limit, not the 60 s cap, must end the search.
  EXPECT_LT(first_s, 30.0);
}

TEST(fields2cover_rp_route_plan_base, best_route_with_params_rejects_bad_text) {
  F2CGraph2D g = completeGraph(4, 3);
  f2c::rp::RoutePlannerBase route_planner;
  EXPECT_THROW(
      route_planner.computeBestRouteWithParams(g, "not a field: 1"),
      std::invalid_argument);
}
