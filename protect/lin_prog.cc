#include "lin_prog.h"

#include <algorithm>
#include <stdexcept>
#include <exception>
#include <iostream>

using std::to_string;
using std::cout;
using std::endl;

lin_prog::lin_prog(string name) {
  this->num_vars = 1;
  this->name = name;
  this->rows.push_back(0);
  this->cols.push_back(0);
  this->vals.push_back(0);
  this->has_run = false;
  this->cur_row = 0;
  this->lp = glp_create_prob();
  glp_set_prob_name(lp, name.c_str());
  glp_set_obj_dir(lp, GLP_MIN); // default to minimize
}

lin_prog::~lin_prog() {
  glp_delete_prob(lp);  
}

bool lin_prog::has(string var) const {
  return std::find(variables.begin(), variables.end(), var) != variables.end();
}

size_t lin_prog::get_offset(string var) const {
  if (!has(var))
    throw std::invalid_argument(var + " does not exist.");
  for (size_t i = 0; i < variables.size(); i++) {
    if (variables[i] == var) {
      return i;
    }
  }
  return variables.size() + 1;
}

std::pair<size_t, size_t> lin_prog::get_bounds(string var) const {
  if (!has(var))
    throw std::invalid_argument(var + " does not exist.");
  for (size_t i = 0; i < variables.size() - 1; i++) {
    if (variables[i] == var) {
      return std::pair<size_t, size_t>(offsets[i], offsets[i+1]-1);
    }
  }
  return std::pair<size_t, size_t>(*offsets.rbegin(), num_vars - 1);
}

void lin_prog::declare_variables(string name, size_t num) {
  if (this->has(name))
    throw std::bad_alloc();
// #ifdef DEBUG
  std::cout << "declaring variable " + name << " with " << num << " indices"
            << endl;
// #endif
  glp_add_cols(lp, num);
  variables.push_back(name);
  offsets.push_back(num_vars);
  for (size_t i = 0; i < num; i++)
    glp_set_col_name(lp, num_vars + i, (name + std::to_string(i)).c_str());
  num_vars += num;
}

void lin_prog::add_constraint(string var, size_t index, double value) {
  const auto bounds = get_bounds(var);
  if (index < 1 || index - 1 + bounds.first > bounds.second) {
    string bounds_string = "(" + std::to_string(bounds.first) + ", " +
                           std::to_string(bounds.second) + ")";
    throw std::invalid_argument(std::to_string(index) + " is out of bounds " +
                                bounds_string + " for " + var);
  }
  if (cur_row < 1)
    throw std::invalid_argument("Must add row before adding constrains");
  rows.push_back(cur_row);
  cols.push_back(bounds.first + (index - 1));
  vals.push_back(value);
}

void lin_prog::set_row_bnd(int type, double lvalue, double rvalue) {
  glp_set_row_bnds(lp, cur_row, type, lvalue, rvalue);
}

void lin_prog::set_var_bnd(string var, size_t index, int type, double lvalue,
                           double rvalue) {
  const auto bounds = get_bounds(var);
  if (index < 1 || (index - 1) + bounds.first > bounds.second) {
    throw std::invalid_argument("[set_var_bnd] " + std::to_string(index) +
                                " is out of bounds for " + var);
  }
  glp_set_col_bnds(lp, bounds.first + (index - 1), type, lvalue, rvalue);
}

void lin_prog::set_var_kind(string var, size_t index, int type) {
  const auto bounds = get_bounds(var);
  if (index < 1 || (index - 1) + bounds.first > bounds.second) {
    throw std::invalid_argument("[set_var_kind] " + std::to_string(index) +
                                " is out of bounds for " + var);
  }
  glp_set_col_kind(lp, bounds.first + (index - 1), type);
}

void lin_prog::set_objective_var(string var, size_t index, double value) {
  const auto bounds = get_bounds(var);
  if (index < 1 || (index - 1) + bounds.first > bounds.second) {
    throw std::invalid_argument("[set_objective_var] " + std::to_string(index) +
                                " is out of bounds for " + var);
  }
  glp_set_obj_coef(lp, bounds.first + (index - 1), value);
}

void lin_prog::set_max() { glp_set_obj_dir(lp, GLP_MAX); }
void lin_prog::set_min() { glp_set_obj_dir(lp, GLP_MIN); }

void lin_prog::add_row() {
  glp_add_rows(lp, 1);
  ++cur_row;
}

void lin_prog::add_row(string name) {
  add_row();
  glp_set_row_name(lp, cur_row, name.c_str());
}

void lin_prog::apply_constraints() {
  glp_load_matrix(lp, (this->rows.size() - 1), &rows[0], &cols[0], &vals[0]);
}

int lin_prog::run(glp_iocp* parm) {
  if (parm != nullptr) {
    parm->presolve = GLP_ON;    
    glp_init_iocp(parm);
  }
  apply_constraints();
  has_run = true;
  return glp_intopt(lp, parm);
}

// return a string representation of this LP
void lin_prog::to_string() const {
  size_t row = 1;
  string new_line = "\n";
  string result = name + new_line + "\tst\t";
  new_line += "\t\t";
  for (size_t i = 0; i < rows.size(); i++) {
    if (rows[i] != row) {
      result += new_line;
      row = rows[i];
    }
    // search for variable name by checking highest offset below index.
    size_t var_index = 0;
    for (; var_index < variables.size() &&
           offsets[var_index] > cols[i]; var_index++) {}
    if (offsets[var_index] > cols[i])
      var_index--;

    result += (vals[i] >= 0 ? " + " : " - ") + std::to_string(vals[i]) +
              variables[var_index];
  }
}

double lin_prog::get_obj_val() const {
  if (!this->has_run)
    throw std::logic_error("LP has to be run before getting objective");
  return glp_get_obj_val(lp);
}


double lin_prog::get_var_val(string var, size_t index) const {
  if (!this->has_run)
    throw std::logic_error("LP has to be run before getting objective");
  const auto bounds = get_bounds(var);
  if (index < 1 || index - 1 + bounds.first > bounds.second)
    throw std::invalid_argument("[get_var_val] " + std::to_string(index) +
                                " is out of bounds for " + var);
  glp_get_col_prim(lp, bounds.first + (index - 1));
}
