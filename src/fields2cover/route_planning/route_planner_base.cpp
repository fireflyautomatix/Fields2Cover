//=============================================================================
//    Copyright (C) 2021-2024 Wageningen University - All Rights Reserved
//                     Author: Gonzalo Mier
//                        BSD-3 License
//=============================================================================

#include <ortools/constraint_solver/routing.h>
#include <ortools/constraint_solver/routing_enums.pb.h>
#include <ortools/constraint_solver/routing_index_manager.h>
#include <ortools/constraint_solver/routing_parameters.h>
#include <google/protobuf/text_format.h>
#include <math.h>
#include <utility>
#include <vector>
#include <limits>
#include <stdexcept>
#include <string>
#include "fields2cover/route_planning/route_planner_base.h"

namespace f2c::rp {

namespace ortools = operations_research;

F2CRoute RoutePlannerBase::genRoute(
    const F2CCells& cells, const F2CSwathsByCells& swaths,
    bool show_log, double d_tol, bool redirect_swaths,
    long int time_limit_seconds, bool search_for_optimum) {
  F2CGraph2D shortest_graph = createShortestGraph(cells, swaths, d_tol);

  F2CGraph2D cov_graph = createCoverageGraph(
      cells, swaths, shortest_graph, d_tol, redirect_swaths);

  std::vector<long long int> v_route = computeBestRoute(
      cov_graph, show_log, time_limit_seconds, search_for_optimum, -1);
  return transformSolutionToRoute(
      v_route, swaths, cov_graph, shortest_graph);
}

void RoutePlannerBase::setStartAndEndPoint(const F2CPoint& p) {
  this->r_start_end = p;
}

F2CGraph2D RoutePlannerBase::createShortestGraph(
    const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
    double d_tol) const {
  F2CGraph2D g;
  // Add points from swaths that touches border
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      g.addEdge(s.startPoint(), cells.closestPointOnBorderTo(s.startPoint()));
      g.addEdge(s.endPoint(),   cells.closestPointOnBorderTo(s.endPoint()));
    }
  }

  // Add points in the border
  for (auto&& cell : cells) {
    for (auto&& ring : cell) {
      for (size_t i = 0; i < ring.size()-1; ++i) {
        g.addEdge(ring.getGeometry(i), ring.getGeometry(i+1));
      }
    }
  }

  // Add start and end point if they exists
  if (this->r_start_end) {
    g.addEdge(*r_start_end, cells.closestPointOnBorderTo(*r_start_end));
  }

  std::vector<F2CPoint> nodes = g.getNodes();

  // Connect nodes that are near other edges.
  for (int i = 0; i < 2; ++i) {
    auto edges = g.getEdges();
    for (auto&& edge : edges) {
      size_t from = edge.first;
      for (auto&& e : edge.second) {
        size_t to = e.first;

        F2CLineString border {g.indexToNode(from), g.indexToNode(to)};
        for (auto&& n : nodes) {
          if (n != border[0] && n != border[1] && n.distance(border) < d_tol) {
            g.addEdge(border[0], n);
            g.addEdge(n, border[1]);
            g.removeEdge(border[0], border[1]);
          }
        }
      }
    }
  }
  return g;
}


F2CGraph2D RoutePlannerBase::createCoverageGraph(
    const F2CCells& cells, const F2CSwathsByCells& swaths_by_cells,
    F2CGraph2D& shortest_graph,
    double d_tol, bool redirect_swaths) const {
  F2CGraph2D g;
  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      F2CPoint mid_p {(s.startPoint() + s.endPoint()) * 0.5};
      if (redirect_swaths) {
        g.addEdge(s.startPoint(), mid_p, 0);
        g.addEdge(s.endPoint(), mid_p, 0);
      } else {
        g.addDirectedEdge(s.startPoint(), mid_p, 0);
        g.addDirectedEdge(mid_p, s.endPoint(), 0);
      }
    }
  }

  for (const auto& swaths1 : swaths_by_cells) {
    for (const auto& s1 : swaths1) {
      auto s1_s = s1.startPoint();
      auto s1_e = s1.endPoint();
      for (const auto& swaths2 : swaths_by_cells) {
        for (const auto& s2 : swaths2) {
          auto s2_s = s2.startPoint();
          auto s2_e = s2.endPoint();
          if (redirect_swaths) {
            g.addEdge(s1_s, s2_s, shortest_graph);
            g.addEdge(s1_e, s2_e, shortest_graph);
            g.addEdge(s1_s, s2_e, shortest_graph);
            g.addEdge(s1_e, s2_s, shortest_graph);
          } else {
            g.addDirectedEdge(s1_e, s2_s, shortest_graph);
          }
        }
      }
    }
  }

  F2CPoint deposit(-1e8, -1e8);  // Arbitrary point
  if (this->r_start_end) {
    deposit = *this->r_start_end;
  }

  for (auto&& swaths : swaths_by_cells) {
    for (auto&& s : swaths) {
      if (this->r_start_end) {
        g.addEdge(s.startPoint(), deposit, shortest_graph);
        g.addEdge(s.endPoint(), deposit, shortest_graph);
      } else {
        g.addEdge(s.startPoint(), deposit, 0);
        g.addEdge(s.endPoint(), deposit, 0);
      }
    }
  }
  return g;
}

std::vector<long long int> RoutePlannerBase::computeBestRoute(
    const F2CGraph2D& cov_graph, bool show_log, long int time_limit_seconds,
    bool use_guided_local_search, long long int depot_id) const {
  if (depot_id < 0) {
    depot_id = static_cast<int>(cov_graph.numNodes()-1);
  }
  // Ensure depot_id is within valid int range
  if (depot_id > std::numeric_limits<int>::max()) {
    depot_id = std::numeric_limits<int>::max();
  }
  const ortools::RoutingIndexManager::NodeIndex depot{static_cast<int>(depot_id)};
  ortools::RoutingIndexManager manager(cov_graph.numNodes(), 1, depot);
  ortools::RoutingModel routing(manager);

  const int transit_callback_index = routing.RegisterTransitCallback(
      [&cov_graph, &manager] (long long int from, long long int to) -> long long int {
        auto from_node = manager.IndexToNode(from).value();
        auto to_node = manager.IndexToNode(to).value();
        return cov_graph.getCostFromEdge(from_node, to_node);
      });
  routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);
  ortools::RoutingSearchParameters searchParameters =
    ortools::DefaultRoutingSearchParameters();
  searchParameters.set_use_full_propagation(false);
  searchParameters.set_first_solution_strategy(
    ortools::FirstSolutionStrategy::AUTOMATIC);
  if (use_guided_local_search) {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
  } else {
    searchParameters.set_local_search_metaheuristic(
      ortools::LocalSearchMetaheuristic::AUTOMATIC);
  }
  searchParameters.mutable_time_limit()->set_seconds(time_limit_seconds);
  searchParameters.set_log_search(show_log);
  const ortools::Assignment* solution =
    routing.SolveWithParameters(searchParameters);

  long long int index = routing.Start(0);
  std::vector<long long int> v_id;

  index = solution->Value(routing.NextVar(index));

  while (!routing.IsEnd(index)) {
    v_id.emplace_back(manager.IndexToNode(index).value());
    index = solution->Value(routing.NextVar(index));
  }
  return v_id;
}

std::vector<long long int> solveRouting(
    const F2CGraph2D& cov_graph, const std::vector<int64_t>& cost,
    const ortools::RoutingSearchParameters& params, bool remove_missing_arcs,
    long long int stall_solution_limit, int depot) {
  const int64_t n = static_cast<int64_t>(cov_graph.numNodes());
  ortools::RoutingIndexManager manager(
      static_cast<int>(n), 1, ortools::RoutingIndexManager::NodeIndex{depot});
  ortools::RoutingModel routing(manager);
  const int transit_callback_index = routing.RegisterTransitCallback(
      [&cost, &manager, n] (int64_t from, int64_t to) -> int64_t {
        return cost[manager.IndexToNode(from).value() * n +
                    manager.IndexToNode(to).value()];
      });
  routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index);

  if (remove_missing_arcs) {
    routing.CloseModelWithParameters(params);
    for (int64_t idx = 0; idx < routing.Size(); ++idx) {
      const int64_t from = manager.IndexToNode(idx).value();
      std::vector<int64_t> allowed;
      for (const size_t to : cov_graph.getEdgesFrom(from)) {
        if (static_cast<int64_t>(to) == from ||
            cost[from * n + to] == std::numeric_limits<int64_t>::max()) {
          continue;
        } else if (static_cast<int>(to) == depot) {
          allowed.push_back(routing.End(0));
        } else {
          allowed.push_back(manager.NodeToIndex(
              ortools::RoutingIndexManager::NodeIndex{static_cast<int>(to)}));
        }
      }
      routing.NextVar(idx)->SetValues(allowed);
    }
  }

  int64_t best = std::numeric_limits<int64_t>::max();
  int64_t since_best = 0;
  if (stall_solution_limit > 0) {
    routing.AddAtSolutionCallback([&routing, &best, &since_best]() {
      const int64_t c = routing.CostVar()->Max();
      if (c < best) {
        best = c;
        since_best = 0;
      } else {
        ++since_best;
      }
    });
    routing.AddSearchMonitor(routing.solver()->MakeCustomLimit(
        [&since_best, stall_solution_limit]() {
          return since_best >= stall_solution_limit;
        }));
  }

  const ortools::Assignment* solution = routing.SolveWithParameters(params);
  std::vector<long long int> v_id;
  if (solution == nullptr) {
    return v_id;
  }
  int64_t index = solution->Value(routing.NextVar(routing.Start(0)));
  while (!routing.IsEnd(index)) {
    v_id.emplace_back(manager.IndexToNode(index).value());
    index = solution->Value(routing.NextVar(index));
  }
  return v_id;
}

std::vector<long long int> RoutePlannerBase::computeBestRouteWithParams(
    const F2CGraph2D& cov_graph, const std::string& search_params_text,
    bool remove_missing_arcs, long long int stall_solution_limit,
    long long int depot_id) const {
  const int64_t n = static_cast<int64_t>(cov_graph.numNodes());
  if (depot_id < 0) {
    depot_id = n - 1;
  }
  if (depot_id > std::numeric_limits<int>::max()) {
    depot_id = std::numeric_limits<int>::max();
  }
  ortools::RoutingSearchParameters params =
    ortools::DefaultRoutingSearchParameters();
  params.set_use_full_propagation(false);
  params.set_first_solution_strategy(
    ortools::FirstSolutionStrategy::AUTOMATIC);
  params.set_local_search_metaheuristic(
    ortools::LocalSearchMetaheuristic::GUIDED_LOCAL_SEARCH);
  params.mutable_time_limit()->set_seconds(10);
  if (!search_params_text.empty()) {
    ortools::RoutingSearchParameters overrides;
    if (!google::protobuf::TextFormat::ParseFromString(
          search_params_text, &overrides)) {
      throw std::invalid_argument(
          "computeBestRouteWithParams: cannot parse search_params_text");
    }
    params.MergeFrom(overrides);
    if (overrides.has_time_limit()) {
      *params.mutable_time_limit() = overrides.time_limit();
    }
  }

  // Dense copy of the arc costs: the transit callback is the solver's hot path.
  std::vector<int64_t> cost(n * n);
  for (int64_t i = 0; i < n; ++i) {
    for (int64_t j = 0; j < n; ++j) {
      cost[i * n + j] = cov_graph.getCostFromEdge(i, j);
    }
  }
  std::vector<long long int> route = solveRouting(
      cov_graph, cost, params, remove_missing_arcs, stall_solution_limit,
      static_cast<int>(depot_id));
  if (route.empty() && remove_missing_arcs) {
    route = solveRouting(cov_graph, cost, params, false, stall_solution_limit,
                         static_cast<int>(depot_id));
  }
  return route;
}

F2CRoute RoutePlannerBase::transformSolutionToRoute(
    const std::vector<long long int>& route_ids,
    const F2CSwathsByCells& swaths_by_cells,
    const F2CGraph2D& coverage_graph,
    F2CGraph2D& shortest_graph) const {
  F2CRoute route;
  F2CSwath swath;
  const size_t NS = swaths_by_cells.sizeTotal();
  for (int i = 0; i < route_ids.size()-2; ++i) {
    F2CPoint p_s = coverage_graph.indexToNode(route_ids[i]);
    F2CPoint p_e = coverage_graph.indexToNode(route_ids[i+2]);
    for (int j = 0; j < NS; ++j) {
      F2CSwath swath = swaths_by_cells.getSwath(j).clone();
      if (p_s == swath.startPoint() && p_e == swath.endPoint()) {
        if (route.isEmpty() && r_start_end) {
          route.addConnection(shortest_graph.shortestPath(
                *r_start_end, swath.startPoint()));
        }
        route.addSwath(swath, shortest_graph);
        break;
      } else if (p_e == swath.startPoint() && p_s == swath.endPoint()) {
        swath.reverse();
        if (route.isEmpty() && r_start_end) {
          route.addConnection(shortest_graph.shortestPath(
                *r_start_end, swath.startPoint()));
        }
        route.addSwath(swath, shortest_graph);
        break;
      }
    }
  }
  if (r_start_end) {
    route.addConnection(shortest_graph.shortestPath(
          route.endPoint(), *r_start_end));
  }
  return route;
}

}  // namespace f2c::rp


