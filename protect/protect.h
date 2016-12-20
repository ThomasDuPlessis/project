#ifndef PROTECT_H
#define PROTECT_H

#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;

struct Activity {
  int number; /** */
  int time;
  double effectiveness;
  /** activities are sorted based on their  */
  bool operator<(const Activity& other) const {
    return this->number < other.number;
  }
};

// PatrolArea is a set of targets.
typedef vector<int> PatrolArea;

struct Patrol {
  size_t area_num;
  Activity activity;
  Patrol(const int area, const Activity &activity)
      : area_num(area), activity(activity) {}
};

typedef std::vector<Patrol> PatrolSchedule;

/* 
 * Structure containing necessary data to create a strategy on defending
 * targets.
 *
 * NOTE: Each of the four reward/penalty vectors should be the same size, the
 * number of targets.
 */
struct ProtectData {
  vector<PatrolArea> PatrolAreas; // Set of patrol areas a defender can visit.
  vector<int> d_rewards;   // Reward for the defender to successfully defend
                           // each target.
  vector<int> d_penalties; // Penalty for the defender to fail to defend each
                           // target.
  vector<int> a_rewards;   // Reward for the attacker to successfully attack
                           // each target.
  vector<int> a_penalties; // Reward for the attacker to fail to attack each
                           // target.
  vector<Activity> activities; // Defender activities.
};

/** Enumerate all possible compact strategies, creating, essentially, the game
 * matrix.
*/
std::vector<PatrolSchedule>
generate_compact_strategies(const int time, const ProtectData &data);

/** Reduce a set of schedules to their compace representation.
 */
void reduce_schedules(std::vector<PatrolSchedule> &schedules);

std::vector<double>
create_strategy(const std::vector<PatrolSchedule> &schedules,
                const ProtectData &data);

void print_schedules(const std::vector<PatrolSchedule> &schedules);

#endif /* PROTECT_H */
