#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <glpk.h>

using std::string;

class lin_prog {
private:
  size_t num_vars;
  size_t cur_row;
  std::vector<int> rows;
  std::vector<int> cols;
  std::vector<double> vals;
  std::vector<string> variables;
  std::vector<size_t> offsets;
  std::string name;
  bool has_run;
  glp_prob *lp;

  /** 
   * Apply constraints to linear program
   *
   */
  void apply_constraints();
  
  /** 
   * Return offset of the variable in our LP. this index should be the index of
   * var indexed at 0 in the gplk object.
   *
   * @param var name of variable
   *
   * @return offset of the column for variable
   */
  size_t get_offset(string var) const;

  /** 
   * Return the bounds of the variable in our LP. The first index should be the
   * index of var indexed at 0 in the gplk object and the second should be the
   * last index.
   *
   * @param var name of variable
   *
   * @return upper and lower bounds (inclusive) of columns var uses
   */
  std::pair<size_t, size_t> get_bounds(string var) const;

  /** 
   *  return true if variable name has been declared.
   *
   * @param var variable name
   *
   * @return bool if var has been declared
   */
  bool has(string var) const;

public:

  /** 
   * linear program object constructor.
   *
   * @param name name of the LP
   */
  lin_prog(string name);

  /** 
   * Deconstructor (deletes lp obj)
   *
   */
  ~lin_prog();

  /** 
   * declare a variable to be used in the LP
   *
   * @param name name of the variable 
   * @param num amount of sub variables for this variable.
   */
  void declare_variables(string name, size_t num);


  /** 
   * Adds constraint to current row, at index i.
   *
   * @param var string name of the argument
   * @param index index of the subvariable in var
   * @param value coefficient of the variable in current row constraint
   */
  void add_constraint(string var, size_t index, double value);

  /** 
   * Set the bounds for the curent row
   *
   * @param type  type of the bound (upper, lower fixed)
   * @param lvalue lower bound
   * @param rvalue upper bound
   */
  void set_row_bnd(int type, double lvalue, double rvalue);

  /** 
   * Set the bounds for variable var
   *
   * @param type  type of the bound (upper, lower fixed)
   * @param lvalue lower bound
   * @param rvalue upper bound
   */
  void set_var_bnd(string var, size_t index, int type, double lvalue, double rvalue);

  // Sets coefficient for variable var at index to value in objective function
  void set_objective_var(string var, size_t index, double value);

  /** 
   * Set the kind of a variable (for when variable needs to be something other
   * then a continuous real variable)
   *
   * @param var name of the variable being set
   * @param index index of the sub variable
   * @param type can be GLP_CV, GLP_IV, GLP_BV
   */
  void set_var_kind(string var, size_t index, int type);
  

  // run mixed integer optimization on the linear program.
  int run(glp_iocp* parm);

 /** 
  *  Add a new row of constraints
  *
  */
  void add_row();

  /** 
   * Add a new row, given name name.
   *
   * @param name name of the row 
   */
  void add_row(string name);

  // set objective to maximize.
  void set_max();

  // set object to minimize
  void set_min();

  // return a string representation of this LP
  void to_string() const;

  double get_obj_val() const;

  double get_var_val(string var, size_t index) const;

};
