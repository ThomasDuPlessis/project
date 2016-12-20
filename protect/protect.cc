#include <algorithm>
#include <vector>

#include "PASAQ.h"
#include "protect.h"

void print_schedules(const std::vector<PatrolSchedule> &schedules) {
  int count = 0;
  for (const auto &schedule : schedules) {
    cout << "|" << count++ << "|";
    for (const auto &area_act : schedule)
      cout << "(" << area_act.area_num << ":k_" << area_act.activity.number << ")";
    cout << " |" << endl;
  }
}

/** 
 * Generate compact schedule up to Size N.
 */
std::vector<std::vector<int>>
generate_compact_schedules(const int n, const ProtectData& data) {
  std::vector<std::vector<int>> result;
  result.emplace_back();
  for (int i = 0; i < n; i++) {
    std::vector<std::vector<int>> result_extension;
    for (const auto& schedule : result) {
      result_extension.push_back(schedule);
      result_extension.back().push_back(i);
    }
    result.insert(result.end(), result_extension.begin(),
                  result_extension.end());
  }
  return result;
}

std::vector<PatrolSchedule>
create_compact_strategies(const std::vector<int> compact_schedule,
                          const ProtectData &data) {
  std::vector<PatrolSchedule> patrol_schedules(1);
  for (const int &area : compact_schedule) {
    std::vector<PatrolSchedule> new_patrols;
    for (const auto &patrol_schedule : patrol_schedules) {
      for (const auto &activity : data.activities) {
        new_patrols.push_back(patrol_schedule);
        new_patrols.back().emplace_back(area, activity);
      }
    }
    patrol_schedules = new_patrols;
  }
  return patrol_schedules;
}

/**
 * Given a set of compact schedules, create all possible strategies given
 * possible defensive activites. Schedules are compact if they are ordered and
 * do not repeat.
*/
std::vector<PatrolSchedule>
create_compact_strategies(const std::vector<std::vector<int>> compact_schedules,
                          const ProtectData &data) {
  std::vector<PatrolSchedule> strategies;

  for (const auto &schedule : compact_schedules) {
    const auto schedule_strategies =
      create_compact_strategies(schedule, data);
    strategies.insert(strategies.end(), schedule_strategies.begin(),
                      schedule_strategies.end());
  }
  return strategies;
}

std::vector<PatrolSchedule>
generate_compact_strategies(const int time, const ProtectData &data) {
  const auto &min_activity = std::min_element(
      data.activities.begin(), data.activities.end());
  int n_hat = time / min_activity->time;

  std::cout << "Longest possible schedule is " << n_hat << " stops long"
            << endl;

  auto schedules = generate_compact_schedules(n_hat, data);
  if (schedules[0].size() < 1) {
    schedules.erase(schedules.begin());
  }

  std::cout<< "Generated schedules, now creating strategies" << std::endl;
  const auto strategies = create_compact_strategies(schedules, data);

  return strategies;
}

// Reduce a schedule by removing repeat repeat nodes, keeping the one with the
// bigger payoff.
void reduce_schedule(PatrolSchedule &schedule) {
  for (size_t i = 0; i < schedule.size(); i++) {
    const auto &patrol = schedule[i];
    for (size_t j = 0; j < schedule.size(); j++) {
      if (i != j && patrol.area_num == schedule[j].area_num) {
        // The has a greater payoff, so delete the other.
        if (patrol.activity.effectiveness >=
            schedule[j].activity.effectiveness) {
          schedule.erase(schedule.begin() + j);
        } else {
          schedule.erase(schedule.begin() + i);
          break;
        }
      }
    }
  }
}

bool schedule_equals(const PatrolSchedule &s1, const PatrolSchedule &s2) {
  // Get max area.
  int max_area =
      max_element(s1.begin(), s1.end(), [](const Patrol &a, const Patrol &b) {
        return a.area_num < b.area_num;
      })->area_num;

  vector<int> activities(max_area + 1);
  for (const auto &patrol : s1)
    activities[patrol.area_num] = patrol.activity.number;
  for (const auto &patrol : s2)
    if (patrol.area_num > activities.size() ||
        activities[patrol.area_num] != patrol.activity.number)
      return false;

  return true;
}

void reduce_schedules(std::vector<PatrolSchedule> &schedules) {
  // filter out repeat areas
  for (PatrolSchedule &schedule : schedules)
    reduce_schedule(schedule);

  // Remove duplicates.
  for (size_t i = 0; i < schedules.size(); i++)
    for (size_t j = i; j < schedules.size(); j++)
      if (schedule_equals(schedules[i], schedules[j]))
        schedules.erase(schedules.begin() + j);
}

std::vector<double>
create_strategy(const std::vector<PatrolSchedule> &schedules,
                const ProtectData &data) {
  vector<vector<double>> A(1);

  const int num_targets = data.a_penalties.size();

  cout << "RUNNING PASAQ ON " << schedules.size() << " compact strategies, on "
       << num_targets << " targets" << endl;
  print_schedules(schedules);
  cout << "Effectiveness matrix size " << num_targets << "x" << schedules.size()
       << endl;

  // Initialize probability matrix.
  for (int i = 0; i < num_targets; i++)
    A.push_back(vector<double>(schedules.size(), 0));
  // cout << "A size: " << A.size() <<"x" << A[0].size() << endl;
  for (size_t j = 0; j < schedules.size(); j++) {
    const auto schedule = schedules[j];
    // cout << "schedule of size " << schedule.size() << endl;
    for (const auto &patrol : schedule) {
      // cout << "Targets: " << data.PatrolAreas[patrol.area_num].size() << endl;
      for (const auto target : data.PatrolAreas[patrol.area_num]) {
        A[target][j] += patrol.activity.effectiveness;
      }
    }
  }

  // Print out probability matrix.
  cout << "Effectiveness matrix: " << endl;
  for (const auto &row : A) {
    for (const auto &elem : row)
      cout << elem << "\t ";
    cout << endl;
  }
  
  PayoffMatrix Pm(data.a_rewards, data.a_penalties, data.d_rewards,
                  data.d_penalties);

  cout << "A size: " << A.size() <<"x" << A[0].size() << endl;
  cout << "Using Binary Search Method to Solve PASAQ" << endl;
  const auto result = BinarySearchMethod(0.5, 5, Pm, A, 0.5, 5);

  return result.second;
}
