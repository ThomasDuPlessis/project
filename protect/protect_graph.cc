#include <algorithm>
#include <iostream>
#include <vector>

#include "protect_graph.h"

using namespace std;

void print_paths(const vector<vector<int>> &paths) {
  for (const auto &path : paths) {
    for (const auto &area : path)
      cout << area << ", ";
    cout << endl;
  }
}


vector<vector<int>> cycles_length(const int base, const int length,
                                 const vector<vector<int>> &adjacency_list) {
  auto paths = paths_length(base, length, adjacency_list);

  // Filter out paths that do not end at base (not a cycle).
  paths.resize(std::remove_if(
                   paths.begin(), paths.end(),
                   [&base](const auto path) { return path.front() != base; }) -
               paths.begin());
  cout << "FILTERED PATHS NOT STARTING WITH" << base << endl;
  print_paths(paths);
  return paths;
}


// Return all paths of length, length starting at base.
vector<vector<int>> paths_length(const int base, const int length,
                                 const vector<vector<int>> &adjacency_list) {
  vector<vector<int>> result;
  if (length < 1) {
    result.push_back({ base });
    return result;
  }

  // Path if it doesn't move.
  for (auto &path : paths_length(base, length - 1, adjacency_list)) {
    path.push_back(base);
    result.push_back(path);
  }

  for (const auto &area : adjacency_list[base]) {
    for (auto &path : paths_length(area, length - 1, adjacency_list)) {
      path.push_back(base);
      result.push_back(path);
    }
  }
  return result;
}
