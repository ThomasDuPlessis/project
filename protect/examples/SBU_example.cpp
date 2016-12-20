#include <stdio.h>
#include <iostream>

#include "../protect.h"

int main(int argc, char *argv[]) {

  vector<PatrolArea> patrol_areas({1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                                   12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22});

  vector<vector<PatrolArea>> adjacency_list = {
    {                    },
    {2                   },
    {6, 9                },
    {4, 5                },
    {3, 5                },
    {3, 4, 8, 11, 12, 14 },
    {2, 7, 9             },
    {6, 8, 9             },
    {7, 9, 5             },
    {2, 6, 7, 8, 10, 11  },
    {9, 11               },
    {5, 9, 10, 12, 16    },
    {5, 11, 13, 14, 15   },
    {12, 14              },
    {5, 12, 13, 15       },
    {12, 14, 16          },
    {11, 15, 17          },
    {16, 20              },
    {19                  },
    {20, 21              },
    {17, 19              },
    {19, 22              },
    {21                  } 
  };

  vector<pair<int, int>> area_targets = {
      {1, 5},     {6, 10},    {11, 20},   {21, 30},   {31, 45},   {46, 55},
      {56, 60},   {61, 65},   {66, 75},   {76, 80},   {81, 90},   {91, 95},
      {96, 100},  {101, 105}, {106, 110}, {111, 120}, {121, 125}, {126, 130},
      {131, 135}, {135, 170}, {171, 175}, {176, 180}};

  // 2 actvities, full search with 3/4 effectiveness, and quick search with .5
  // effectiveness.
  vector<Activity> activities({{1, 3, .5}, {2, 5, .75}});

  vector<int> d_rewards  (181);
  vector<int> d_penalties(181);
  vector<int> a_rewards  (181);
  vector<int> a_penalties(181);

  // Catching parking ticket gets 30$, not catchin git is approx losing 15$ fee
  for (int i = 1; i < 181; i++) {
    d_rewards[i] = 30;
    d_penalties[i] = -15;
    a_penalties[i] = -30;
    a_rewards[i] = 15;
  }

  PatrolGraph graph(patrol_areas, adjacency_list, activities, d_rewards,
                    d_penalties, a_rewards, a_penalties, 181, area_targets);

  auto schedules = graph.generate_schedules(10, 10);

  graph.reduce_schedules(schedules);

  std::cout << "REDUCED:" << endl;

  print_schedules(schedules);

  const auto result = graph.create_strategy(schedules);
  cout <<  "strategy: ";
  for (const auto& r : result)
    cout << r << ",";

  cout << endl;

  return 0;
}
