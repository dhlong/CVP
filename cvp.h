#ifndef __CVP_H__
#define __CVP_H__

#include <ilcplex/ilocplex.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <cassert>
#include <iomanip>
#include "cvputility.h"
#include "function.h"
#include "network.h"

ILOSTLBEGIN

class CVP{
 private:
  void default_settings();
  
  // ILOG objects
  IloModel *proxy; // pointer to the proxy model
  IloObjective *proxy_obj; // pointer to the proxy objective
  IloCplex *cplex; // pointer to the cplex (CPLEX solver)

  // Helper functions
  Vector yy(const Vector &x, Real beta);
  bool is_optimal(const Vector &x, const Vector &y);
  
 protected:
  // Settings
  Real optimality_epsilon;
  Real initial_beta;
  Real beta_up_factor;
  Real beta_down_factor;
  bool to_reset_beta;
  bool to_golden_search;
  bool to_line_search;
  bool to_do_shortest_path;
  bool to_do_SOCP;
  int SP_iterations_per_SOCP;
  int SP_iterations;
  int line_search_iterations;

  IloObjective quad_proxy_obj(const Vector &y);
  IloObjective linear_proxy_obj(const Vector &x);

  // Actual objective function
  Function *obj;

  IloEnv env; // ILOG environment
  IloNumVarArray variables; // set of all variables
  IloRangeArray constraints; // set of all constraints

  virtual IloObjective initial_proxy_obj(); 
  virtual Vector initial_solution();
  
  Vector solve(const IloObjective &obj);
  
 public:
  // constructors
  CVP();
  CVP(Function *obj_);

  // destructor
  virtual ~CVP();

  // settings
  void read_settings(fstream &f);  
  
  // optimization procedures
  virtual Vector phase2(Vector *init=NULL);
  virtual Vector phase1(Vector *init=NULL);
  virtual Vector optimize();
};


// CVP_NF is a CVP class designed for single commodity network optimisation
// inherited from CVP class
// redefine the initial_proxy_obj
// adding function to generate flow constraints from network structure
class CVP_NF : public CVP {
private:
  Network net;
  void generate_flow_constraints();
  
protected:
  virtual IloObjective initial_proxy_obj();

public:
  CVP_NF(Function *obj, const Network &n);
  Vector solvelinear();
  Vector solvequad();
};

// CVP_MCNF is a CVP class designed for multi-commodity network optimisation
// inherited from CVP class
// redefine the initial_proxy_obj() and initial_solution()
// adding function to generate constraints from network structure
class CVP_MCNF : public CVP{
private:
  MultiCommoNetwork net;
  void generate_network_constraints();

protected:
  virtual IloObjective initial_proxy_obj();
  virtual Vector initial_solution();

public:
  CVP_MCNF(Function *obj_, MultiCommoNetwork net_);
  Vector solvelinear();
  Vector solve_by_dijkstra();
  Vector solve_by_dijkstra_and_SOCP();

  virtual Vector optimize();
};

#endif




