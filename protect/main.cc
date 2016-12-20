#include <stdio.h>
#include <iostream>

#include "protect.h"

int main(int argc, char *argv[]) {
  vector<PatrolArea> patrol_areas = {{1, 2, 3}, {4, 5, 6},
                                     {7, 8, 9}};
  vector<Activity> activities({{1, 3, .5}, {2, 5, .75}});
  vector<int> d_rewards =   {0 , 30, 30, 30, 30, 30, 30, 30, 30, 30};
  vector<int> d_penalties = {0 , -15, -15, -15, -15, -15, -15, -15, -15, -15};
  vector<int> a_rewards =   {0 , 30, 30, 30, 30, 30, 30, 30, 30, 30};
  vector<int> a_penalties = {0 , -30, -30, -30, -30, -30, -30, -30, -30, -30};
  ProtectData data;
  data.PatrolAreas = patrol_areas;
  data.a_penalties = a_penalties;
  data.a_rewards = a_rewards;
  data.d_penalties = d_penalties;
  data.d_rewards = d_rewards;
  data.activities = activities;
  const auto compact_strats = generate_compact_strategies(10, data);
  const auto result = create_strategy(compact_strats, data);
  std::cout <<  "strategy: ";
  for (const auto& r : result)
    cout << r << ",";
  cout << endl;
  return 0;
}
