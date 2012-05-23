#include "cvp.h"
#include <set>
#include <list>
#include "dijkstra.h"

ILOSTLBEGIN

// report files, defined in main.cpp
extern fstream iteration_report;
extern fstream solve_report;

//////////////////////////////////////////////////////////////
///////// CVP class' function definitions
/////////
//////////////////////////////////////////////////////////////

// default constructor: no objective specified
CVP::CVP():
  env(), 
  variables(env), 
  constraints(env), 
  obj(NULL), 
  cplex(NULL), 
  proxy(NULL), 
  proxy_obj(NULL) { 
  default_settings(); 
}

// constructor from a potin to a objective function
CVP::CVP(Function *obj_) : 
  env(), 
  variables(env), 
  constraints(env),
  obj(obj_), 
  cplex(NULL), 
  proxy(NULL),
  proxy_obj(NULL) { 
  default_settings(); 
}

// default settings
void CVP::default_settings(){
  optimality_epsilon = 1e-7;
  beta_up_factor = 1.1;
  initial_beta = 10.0;
  beta_down_factor = 0.5;
  tau_multiplier = 4.0;
  to_reset_beta = false;
  to_line_search = true;
  to_golden_search = true;
  to_do_shortest_path = true;
  to_do_SOCP = true;
  SP_iterations_per_SOCP = 10;
  SP_iterations = 500;
  line_search_iterations = 20;
}

// Read new settings from file
// If file not exist, no changes made
void CVP::read_settings(fstream &f){
  if(!f){
    iteration_report<<"Default settings"<<endl;
    iteration_report<<"=================================="<<endl;
    return;
  }

  char line[300], *str;

  // parsing each line to get the setting names and values
  // the setting names and setting values are separated by '=' character
  while(!f.eof()){
    f.getline(line, 300);
    str = strrchr(line,'='); // delim by the '=' character
    if(str == NULL) continue; // no delim found, read next line
    stringstream ss(str+1, stringstream::in);

    if(strstr(line, "optimality epsilon") != NULL) 
      ss>>optimality_epsilon;
    else if(strstr(line, "beta up factor") != NULL) 
      ss>>beta_up_factor;
    else if(strstr(line, "beta down factor") != NULL) 
      ss>>beta_down_factor;
    else if(strstr(line, "line search iteration") != NULL) 
      ss>>line_search_iterations;
    else if(strstr(line, "initial beta") != NULL) 
      ss>>initial_beta;
    else if(strstr(line, "reset beta") != NULL)
      to_reset_beta = (strstr(line,"yes") != NULL);
    else if(strstr(line, "line search") != NULL)
      to_line_search = (strstr(line,"yes") != NULL);
    else if(strstr(line, "golden search") != NULL)
      to_golden_search = (strstr(line,"yes") != NULL);
    else if(strstr(line, "to do shortest path") != NULL)
      to_do_shortest_path = (strstr(line,"yes") != NULL);
    else if(strstr(line, "to do SOCP") != NULL)
      to_do_SOCP = (strstr(line,"yes") != NULL);
    else if(strstr(line, "SP iterations per SOCP") != NULL) 
      ss>>SP_iterations_per_SOCP;
    else if(strstr(line, "SP iterations") != NULL) 
      ss>>SP_iterations;
  }

  // report settings to iteration report
  iteration_report<<"Settings:"<<endl;
  iteration_report<<"\toptimality epsilon = "<<optimality_epsilon<<endl;
  iteration_report<<"\tinitial beta = "<<initial_beta<<endl;
  iteration_report<<(to_reset_beta?"\treset beta":"\tno reset beta")<<endl;
  iteration_report<<(to_golden_search?"\tgolden search":"\tnormal search")<<endl;
  iteration_report<<"\tline search iterations = "<<line_search_iterations<<endl;
  iteration_report<<(to_do_SOCP?"\tuse SOCP":"\tno use SOCP")<<endl;
  iteration_report<<(to_do_shortest_path?"\tuse shortest path":"\tno use shortest path")<<endl;
  if(to_do_shortest_path && to_do_SOCP)
    iteration_report<<"\tshortest path replications per iteration = "<<SP_iterations_per_SOCP<<endl;
  FOR(i, 51) iteration_report<<"="; iteration_report<<endl;
}


// default initial proxy objective
// projection the 0-point on the feasible set
IloObjective CVP::initial_proxy_obj(){
  Vector x(variables.getSize(), 0.0);
  return quad_proxy_obj(yy(x,1.0));
}

// default initial solution by solving the initial proxy objective
// can be (virtually) redefined in the inherited classes
Vector CVP::initial_solution(){
  return solve(this->initial_proxy_obj());
}

// return x - alpha * gradient(x)
Vector CVP::yy(const Vector &x, Real alpha){
  Vector z = obj->g(x);
  z *= (-alpha); z += x; // z = -alpha * g(x) + x
  return z;
}

// return the proxy objective min ||x - z||
IloObjective CVP::quad_proxy_obj(const Vector &z){
  return IloMinimize(env,
	       IloScalProd(variables, variables)
	       + (z*z)
	       - 2*IloScalProd(variables, Vector2IloNumArray(env,z)));
}

// return the proxy objective min gradient(x0) * x
IloObjective CVP::linear_proxy_obj(const Vector &x0){
  Vector g = obj->g(x0);
  Real norm = sqrt(g*g);
  g *= (1.0/norm);
  IloNumArray coef = Vector2IloNumArray(env, g);
  return IloMinimize(env,IloScalProd(coef, variables));
}

bool CVP::is_optimal(const Vector &x, const Vector &z){
  Vector d = obj->g(x), y = z;
  y -= x;
  return 1+(d*y)/sqrt((y*y)*(d*d)) <= optimality_epsilon;
}

Vector CVP::solve(const IloObjective &iloobj){
  Timer timer;
  double tic1, tic2, tic3, tic4, tic5, tic6;

  cout<<"New Solve"<<endl;
  tic1 = timer.elapsed();
  
  // If old objective exist, remove it from the proxy model
  if (proxy_obj != NULL) {
    if (proxy != NULL) proxy->remove(*proxy_obj);
    proxy_obj->end();
    delete proxy_obj;
  } 
  
  // create new proxy objective
  proxy_obj = new IloObjective(iloobj);   
  
  // If the proxy model has not been created, create it
  if (proxy == NULL){
    proxy = new IloModel(env);
    proxy->add(constraints);
  }
  tic2 = timer.elapsed();

  // add new proxy objective to the proxy model
  proxy->add(*proxy_obj);
  tic3 = timer.elapsed();

  // If cplex has not been created, create it
  if (cplex == NULL) cplex = new IloCplex(*proxy);
  tic4 = timer.elapsed();

  // solve
  if(!cplex->solve()) {
    cout<<cplex->getStatus()<<endl;
    solve_report<<cplex->getStatus()<<endl;
    return Vector();
  }

  tic5 = timer.elapsed();

  // obtain the optimal solution
  IloNumArray vals(env,variables.getSize());
  cplex->getValues(vals, variables);
  Vector x = IloNumArray2Vector(vals);

  // calculate new real objective (for reporting)
  Real obj_value = (obj==NULL ? -1 : obj->f(x));
  tic6 = timer.elapsed();
  
  // reporting to solve report
  solve_report << setprecision(10) << tic2 - tic1 << "\t"
                << setprecision(10) << tic3 - tic2 << "\t"
                << setprecision(10) << tic4 - tic3 << "\t"
                << setprecision(10) << tic5 - tic4 << "\t"
                << setprecision(10) << tic6 - tic5 << "\t"
                << setprecision(10) << obj_value << "\t"
		<< endl;

  return x;
}
  
CVP::~CVP(){
  env.end();
  //if(proxy) proxy->end();
  //if(cplex) cplex->end();
  if(cplex) delete cplex;
  if(proxy) delete proxy;
  if(proxy_obj) delete proxy_obj;
  if(obj) delete obj;
}

Vector CVP::phase1(Vector *init){
  Real pre_cost, cur_cost;
  Vector x, pre_x;
  Timer timer;
  int iteration = 0;
  double tic1, tic2;

  // Initialisation
  if(init!=NULL && init->size() == variables.getSize()) x = *init;
  else x = solve(this->initial_proxy_obj());
  cur_cost = obj->f(x);

  do{
    iteration++;
    tic1 = timer.elapsed();    
    pre_x = x, pre_cost = cur_cost;
    cur_cost = obj->f(x = solve(linear_proxy_obj(pre_x)));
    tic2 = timer.elapsed();
    
    iteration_report<<"Phase 1\t"<<iteration<<"\t"<<"\t"<<tic2-tic1<<endl;
    if(x==pre_x) return x;
  } while(cur_cost<=pre_cost);

  cout<<"Phase 1 failed"<<endl;
  return Vector();
}

Vector CVP::phase2(Vector *init){
  Real beta;
  int count = 0, iteration = 0;
  Timer timer;
  Real tic1, tic2, tic3;
  bool exit_flag = false;

  tic1 = timer.elapsed();

  // Initialisation by solving initial proxy objective
  Vector x0, x1 = solve(this->initial_proxy_obj());
  Real f0, f1 = obj->f(x1); 

  // Timing and Reporting
  tic2 = timer.elapsed();
  iteration_report<<"Initialisation: time = "<<tic2-tic1<<"s; obj = "<<f1<<endl;
  tic1 = timer.elapsed();

  beta = initial_beta * sqrt(x0*x0); // added by Hieu

  // header row of the iteration report
  FOR(i,71) iteration_report<<"-"; iteration_report<<endl;
  iteration_report << left << setw(6)  << "Iter"
		   << left << setw(8)  << "#solve"
		   << right << setw(8) << "time"
		   << right << setw(10) << "beta"
		   << right << setw(6)  << "ls?"
		   << right << setw(8) << "time_ls"
		   << right << setw(8) << "lambda*"
		   << right << setw(10) << "cosine"
		   << right << setw(17) << "obj"<<endl;
  FOR(i,71) iteration_report<<"-"; iteration_report<<endl;
  
  // Loops
  Vector g, z;
  while(!exit_flag) {
    if(to_reset_beta) beta = initial_beta;

    // Timing and Reporting
    iteration ++; count = 0; tic1 = timer.elapsed();
    solve_report<<"==================== New Iteration ======================"<<endl;

    // Previous best solution
    x0 = x1; f0 = f1;
    
    // Reduce beta and do projection until improvement
    g = obj->g(x0);
    Real normg = sqrt(g*g); // added by Hieu

    for(;;) {
      z = g; z *= (-beta)/normg; z += x0; // z = x0 - beta*g
      x1 = solve(quad_proxy_obj(z));
      f1 = obj->f(x1);
      count++; // counting number of solves (for reporting purpose)
      if(f1 < f0) break;
      beta *= beta_down_factor;
    }

    Real normdx_before_ls = sqrt((x1-x0)*(x1-x0)); // added by Hieu

    // Optimality check
    g = obj->g(x1); // g is now gradient at x1
    z -= x1; // z is now z - x1
    Real cosine = 1 + (z*g)/sqrt((z*z)*(g*g)); //added by Hieu
    exit_flag = cosine <= optimality_epsilon;
    
    tic2 = timer.elapsed();

    // Line Search
    Real lambda = 1.0;
    bool do_line_search = to_line_search && ((x0-x1)*g<0);
    if(do_line_search){
      if(to_golden_search) lambda = golden_search(x0, x1, obj, line_search_iterations);
      else lambda = line_search(x0, x1, obj, line_search_iterations);
      x1 -= x0; x1 *= lambda; x1 += x0; // x1 = x0 + betamin*(x1-x0)
      f1 = obj->f(x1);
    }

    Real normdx_after_ls = sqrt((x1-x0)*(x1-x0)); // added by Hieu
    
    // Timing and Reporting
    tic3 = timer.elapsed();
    iteration_report << left << setw(6)  << iteration
		     << left << setw(8)  << count
		     << right << setw(8) << setprecision(4) << fixed << tic3-tic1
		     << right << setw(10) << setprecision(5) << fixed << beta
		     << right << setw(6) << (do_line_search ? "YES":"NO")
		     << right << setw(8) << setprecision(4) << fixed << tic3-tic2
		     << right << setw(8) << setprecision(3) << fixed << lambda
		     << right << setw(10) << setprecision(5) << scientific << cosine
		     << right << setw(17) << setprecision(8) << scientific << f1 << endl;
  }
  
  // Reporting final results
  FOR(i,71) iteration_report<<"-"; iteration_report<<endl;
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;

  return x1;
}


// main function to call for optimisation
Vector CVP::optimize(){
  //Vector x = phase1();
  //if(x.size() == variables.getSize()) return x;
  return phase2();
}


//////////////////////////////////////////////////////////////
///////// CVP_NF class' function definitions
/////////
//////////////////////////////////////////////////////////////

// constructor
CVP_NF::CVP_NF(Function *obj, const Network &n): CVP(obj), net(n) {
  generate_flow_constraints();
}

// constraint generation
void CVP_NF::generate_flow_constraints(){
  int V = net.getNVertex(), A = net.arcs.size();
  IloExprArray node(env);
  FOR(j,V) node.add(IloExpr(env));

  FOR(i,A) {
    // Name the variable
    stringstream ssname;
    ssname<<"flow["<<net.arcs[i].head<<"->"<<net.arcs[i].tail<<"]";

    // Add the varibale to the variable set
    variables.add(IloNumVar(env, 0.0, net.arcs[i].cap, ssname.str().c_str()));

    // Expression for the left hand sides of the flow equations
    node[net.arcs[i].head] -= variables[i];
    node[net.arcs[i].tail] += variables[i];
  }
  
  // Set of vertex which is either supply or demand node
  set<Vertex> S;
  
  // Generate flow constraint for supply node and add to the constraint set
  FOR(i, net.supply.size()){
    constraints.add(node[net.supply[i].first] >= -net.supply[i].second);
    S.insert(net.supply[i].first);
  }

  // Generate flow constraint for demand node and add to the constraint set
  FOR(i, net.demand.size()){
    constraints.add(node[net.demand[i].first] >= net.demand[i].second);
    S.insert(net.demand[i].first);
  }
  
  // Generate constraint for normal node and add to the constraint set
  FOR(j, V) if(S.count(j)==0) constraints.add(node[j] == 0.0);
}

IloObjective CVP_NF::initial_proxy_obj(){
  IloExpr expr(env);
  Real coef;
  FOR(i, variables.getSize()){
    coef = 1./net.arcs[i].cap;
    coef *= coef*100;
    expr += coef * (variables[i]*variables[i]);

    coef = Real(net.arcs[i].head)/Real(net.arcs[i].tail);
    expr += coef * variables[i];
  }
  return IloMinimize(env, expr);
}

Vector CVP_NF::solvelinear(){    
  int A = net.arcs.size();
  IloExpr expr(env);
  FOR(a, A) expr += (net.arcs[a].cost * variables[a]);
  return solve(IloMinimize(env, expr));
}

Vector CVP_NF::solvequad(){    
  cout<<"Solve quad"<<endl;
  int A = net.arcs.size();
  IloExpr expr(env);
  FOR(a, A) expr += (net.arcs[a].cost * IloSquare(variables[a]));
  return solve(IloMinimize(env, expr));
}

//////////////////////////////////////////////////////////////
///////// CVP_MCNF class' function definitions
/////////
//////////////////////////////////////////////////////////////

// constructor
CVP_MCNF::CVP_MCNF(Function *obj_, MultiCommoNetwork net_)
  : CVP(obj_), net(net_) {
}

Vector CVP_MCNF::optimize(){
  if (!to_do_SOCP){
    env.end();
    return solve_by_dijkstra();
  }
  generate_network_constraints();  
  if (!to_do_shortest_path) return phase2();
  return solve_by_dijkstra_and_SOCP();  
}

// constraint generation
void CVP_MCNF::generate_network_constraints(){
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  IloExprArray node(env);
  FOR(i, K) FOR(n, V) node.add(IloExpr(env));

  FOR(a,A){
    IloExpr sum(env);
    FOR(i, K){
      // Name the variable
      stringstream ssname;
      ssname<<"flow["<<i<<","<<net.arcs[a].head<<">"<<net.arcs[a].tail<<"]";

      // Add the variable to the set of variables
      variables.add(IloNumVar(env, 0.0, IloInfinity, ssname.str().c_str()));
    }

    FOR(i, K){
      // Update flow equation for each commodity at each arc with the new var
      node[i*V + net.arcs[a].head] -= variables[a*K + i];
      node[i*V + net.arcs[a].tail] += variables[a*K + i];
    }
  }

  // Generate flow constraints and add them to the constraint set
  FOR(i,K) FOR(n, V){
    if(n == net.commoflows[i].origin)
      constraints.add(node[i*V + n] == - net.commoflows[i].demand);
    else if(n == net.commoflows[i].destination)
      constraints.add(node[i*V + n] == + net.commoflows[i].demand);
    else
      constraints.add(node[i*V + n] == 0);
  }
  
  // Reporting problem size to the iteration report
  iteration_report << endl << "Problem size:"<<endl;
  iteration_report << "\tNumber of network nodes  = "<<setw(10)<<right<<net.getNVertex()<<endl;
  iteration_report << "\tNumber of network arcs   = "<<setw(10)<<right<<net.arcs.size()<<endl;
  iteration_report << "\tNumber of commodities    = "<<setw(10)<<right<<net.commoflows.size()<<endl;
  iteration_report << "\tNumber of variables      = "<<setw(10)<<right<<variables.getSize()<<endl;
  iteration_report << "\tNumber of constraints    = "<<setw(10)<<right<<constraints.getSize()<<endl;
  iteration_report << endl;
}

// redefine the initial proxy objective 
// to be the linear cost objective
IloObjective CVP_MCNF::initial_proxy_obj(){
  int K = net.commoflows.size(), A = net.arcs.size();
  IloExpr expr(env);
  assert(variables.getSize() == K*A);
  FOR(a, A) FOR(i, K) expr += (net.arcs[a].cost * variables[a*K + i]);
  return IloMinimize(env, expr);
}

// redefine the initial solution to be the combination of solutions 
// when solving single-commodity network for each commodity one by one
Vector CVP_MCNF::initial_solution(){
  int K = net.commoflows.size(), A = net.arcs.size();
  Vector x(variables.getSize());
  FOR(k, K){
    Network newnet;
    newnet.arcs = net.arcs;
    newnet.supply.push_back(make_pair(net.commoflows[k].origin, net.commoflows[k].demand));
    newnet.demand.push_back(make_pair(net.commoflows[k].destination, net.commoflows[k].demand));
    CVP_NF cvp(NULL, newnet);
    Vector xx = cvp.solvequad();
    FOR(a, A) x[a*K + k] = xx[a];
  }
  return x;
}

Vector CVP_MCNF::solvelinear(){    
  int K = net.commoflows.size(), A = net.arcs.size();
  IloExpr expr(env);
  assert(variables.getSize() == K*A);
  FOR(a, A) FOR(i, K) expr += (net.arcs[a].cost * variables[a*K + i]);
  return solve(IloMinimize(env, expr));
}

Vector CVP_MCNF::solve_by_dijkstra(){
  int count = 0, iteration = 0;
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  Timer timer;
  Real tic1, tic2, tic3, tic4, start = timer.elapsed();
  bool exit_flag = false;
  int output_row_width = 57+12;
  
  //////////////////////////////////////////////////////////////////////////
  /////////// Adjacent List Initialisation
  ///////////
  AdjacentList adjl;
  vector< vector<arc_t> > indexarcl(V), indexadjl(V);
  vector< vertex_t * > trace(V);
  vector< char * > vb(V);
  vector< int > nv(V, 0);

  adjl.V = V; adjl.A = A;
  malloc_adjl(&adjl);
  FOR(i, V) trace[i] = MALLOC(vertex_t, V);
  FOR(i, V) vb[i] = MALLOC(char, V);

  FOR(i,V)  
    indexarcl[i] = vector<arc_t>(V, -1), 
    indexadjl[i] = vector<arc_t>(V, -1);
    
  FOR(i, V+1) adjl.n_arcs[i] = 0;
  FOR(a, A)   adjl.n_arcs[net.arcs[a].head]++;
  FOR(i, V)   adjl.n_arcs[i+1] += adjl.n_arcs[i];
    
  FOR(a, A) {
    int aa =  --adjl.n_arcs[net.arcs[a].head];
    int u = net.arcs[a].head, v = net.arcs[a].tail;
    indexarcl [u] [v] = a;
    indexadjl [u] [v] = aa;
    adjl.adjacent_vertices[aa] = net.arcs[a].tail;
  }
  
  FOR(i, V) FOR(j, V) vb[i][j]= 0;
  FOR(k, net.commoflows.size()) {
    vb[net.commoflows[k].origin][net.commoflows[k].destination] = 1;
    nv[net.commoflows[k].origin]++;
  }  
  ///////////
  /////////// end of Adjacent List initialisation
  ////////////////////////////////////////////////////////////////////////

  tic1 = timer.elapsed();

  // Initialisation by solving the intial network shortest paths
  FOR(a, A){
    NetworkArc arc = net.arcs[a];
    adjl.costs[indexadjl[arc.head][arc.tail]] = cost_t(arc.cost);
  }
  FOR(i, V) dijkstra(adjl, i, vb[i], nv[i], trace[i]);
  
  Vector x(A*K, 0.0);
  FOR(k, K){
    int u = net.commoflows[k].origin, v = net.commoflows[k].destination;
    Real demand = net.commoflows[k].demand;
    while(v>=0 && v!=u){
      x[indexarcl[trace[u][v]][v]*K + k] = demand;
      v = trace[u][v];
    }
  }

  // header row of the iteration report
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  iteration_report << left << setw(5)  << "Iter"
		   << right << setw(8) << "t_total"
		   << right << setw(8) << "t_SP"
		   << right << setw(8) << "t_LS"
		   << right << setw(8) << "tau"
		   << right << setw(20) << "obj"
		   << right << setw(12) << "t_elapsed"
		   << endl;
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  
  // Loops
  while(!exit_flag) {
    iteration++;
    cout<<"Iteration "<<iteration<<endl;
    tic1 = timer.elapsed();

    Vector g = obj->g(x);
    FOR(a, A){
      int u = net.arcs[a].head, v = net.arcs[a].tail;
      adjl.costs[indexadjl[u][v]] = cost_t(g[indexarcl[u][v]*K]);
    }
      
    FOR(i, V) if(nv[i]>0) dijkstra(adjl, i, vb[i], nv[i], trace[i]);
    
    Vector sp(x.size(), 0);
    FOR(k, K){
      int u = net.commoflows[k].origin, v = net.commoflows[k].destination;
      Real demand = net.commoflows[k].demand;
      while(v>=0 && v!=u) {
	sp[ indexarcl [trace[u][v]] [v] * K + k] = demand;
	v = trace[u][v];
      }
    }
      
    tic2 = timer.elapsed();
      
    Vector dx(sp); dx -= x;
    Vector gsp = obj->g(sp);
    Real f, fsp, tau = 0.0;
    if((gsp*dx)*(g*dx)<0) {
      tau = golden_search(x, sp, obj, line_search_iterations);
      x -= sp; x *= (1-tau); x += sp;
      f = obj->f(x);
    } 
    else if( (fsp = obj->f(sp)) < f ) x = sp, f = fsp, tau = 1.0;
    else break;

    tic3 = timer.elapsed();
  
    // Timing and Reporting
    tic4 = timer.elapsed();
    iteration_report << left << setw(5)  << iteration
		     << right << setw(8) << setprecision(4) << fixed << tic3-tic1
		     << right << setw(8) << setprecision(4) << fixed << tic2-tic1
		     << right << setw(8) << setprecision(4) << fixed << tic3-tic2
		     << right << setw(8) << setprecision(4) << fixed << tau
		     << right << setw(20) << setprecision(10) << scientific << f 
		     << right << setw(12) << setprecision(3) << fixed << timer.elapsed()-start
		     << endl;

    exit_flag = iteration >= SP_iterations;
  }
  
  // Reporting final results
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x)<<endl;
  
  cout<<"Release adjl"<<endl;
  free_adjl(&adjl);
  cout<<"Relase vb"<<endl;
  FOR(i, V) free(vb[i]);
  cout<<"Relase trace"<<endl;
  FOR(i, V) free(trace[i]);

  return x;
}

Vector CVP_MCNF::solve_by_dijkstra_and_SOCP(){
  Real beta;
  int count = 0, iteration = 0;
  Timer timer;
  Real tic1, tic2, tic3, tic4, start = timer.elapsed();
  bool exit_flag = false;

  //////////////////////////////////////////////////////////////////////////
  /////////// Adjacent List Initialisation
  ///////////
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  AdjacentList adjl;
  vector< vector<arc_t> > indexarcl(V);
  vector< vector<arc_t> > indexadjl(V);
  vector< vertex_t * > trace(V);
  vector< char * > vb(V);
  vector< int > nv(V, 0);
  int output_row_width = 117;

  adjl.V = V; adjl.A = A;
  malloc_adjl(&adjl);
  FOR(i, V) trace[i] = MALLOC(vertex_t, V);
  FOR(i, V) vb[i] = MALLOC(char, V);
  
  FOR(i,V)  
    indexarcl[i] = vector<arc_t>(V, -1), 
    indexadjl[i] = vector<arc_t>(V, -1);
  
  FOR(i, V+1) adjl.n_arcs[i] = 0;
  FOR(a, A)   adjl.n_arcs[net.arcs[a].head]++;
  FOR(i, V)   adjl.n_arcs[i+1] += adjl.n_arcs[i];
    
  FOR(a, A) {
    int aa =  --adjl.n_arcs[net.arcs[a].head];
    int u = net.arcs[a].head, v = net.arcs[a].tail;
    indexarcl [u] [v] = a;
    indexadjl [u] [v] = aa;
    adjl.adjacent_vertices[aa] = net.arcs[a].tail;
  }
  
  FOR(i, V) FOR(j, V) vb[i][j]= 0;
  FOR(k, net.commoflows.size()) {
    vb[net.commoflows[k].origin][net.commoflows[k].destination] = 1;
    nv[net.commoflows[k].origin]++;
  }  
  ///////////
  /////////// end of Adjacent List initialisation
  ////////////////////////////////////////////////////////////////////////

  tic1 = timer.elapsed();

  // Initialisation by solving the intial network shortest paths
  FOR(a, A){
    NetworkArc arc = net.arcs[a];
    adjl.costs[indexadjl[arc.head][arc.tail]] = cost_t(arc.cost);
  }
  FOR(i, V) dijkstra(adjl, i, vb[i], nv[i], trace[i]);
  
  Vector x0, x1(A*K, 0.0);
  FOR(k, K){
    int u = net.commoflows[k].origin, v = net.commoflows[k].destination;
    Real demand = net.commoflows[k].demand;
    while(v>=0 && v!=u){
      x1[indexarcl[trace[u][v]][v]*K + k] = demand;
      v = trace[u][v];
    }
  }
  Real f0, f1 = obj->f(x1); 
  beta = initial_beta * sqrt(x1*x1);

  // Timing and Reporting
  tic2 = timer.elapsed();
  iteration_report<<"Initialisation: time = "<<tic2-tic1<<"s; obj = "<<f1<<endl;
  tic1 = timer.elapsed();

  // header row of the iteration report
  output_row_width = 130;
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  iteration_report << left << setw(4)  << "Iter"
		   << right << setw(7) << "#solve"
		   << right << setw(8) << "t_total"
		   << right << setw(12) << "beta"
		   << right << setw(5) << "ls?"
		   << right << setw(8) << "lambda*"
		   << right << setw(8) << "tau*0"
		   << right << setw(8) << "tau*n"
		   << right << setw(8) << "t_SP"
		   << right << setw(20) << "obj_ls"
		   << right << setw(20) << "obj_final"
		   << right << setw(10) << "cosine"
		   << right << setw(12) << "t_elapsed"
		   <<endl;
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  
  // Loops
  Vector g, z;
  Real taubound = -1, taustar = 0.5/2;
  while(!exit_flag) {
    if(to_reset_beta) beta = initial_beta;

    // Timing and Reporting
    iteration ++; count = 0; tic1 = timer.elapsed();
    solve_report<<"==================== New Iteration ======================"<<endl;

    // Previous best solution
    x0 = x1; f0 = f1;
    
    // Reduce beta and do projection until improvement
    g = obj->g(x0);
    //Vector gg = obj->gg(x0);
    //FOR(i, g.size()) if(gg[i]>1e-6) g[i] /= gg[i]; else g[i] = 0.0;
    Real normg = sqrt(g*g);
    cout<<normg;
    g *= (1/normg);

    for(;;) {
      z = g; z *= (-beta); z += x0; // z = x0 - beta*g
      x1 = solve(quad_proxy_obj(z));
      f1 = obj->f(x1);
      count++; // counting number of solves (for reporting purpose)
      if(f1 < f0) break;
      beta *= beta_down_factor;
    }

    // Optimality check
    g = obj->g(x1); // g is now gradient at x1
    z -= x1; // z is now z - x1
    Real cosine = 1 + (z*g)/sqrt((z*z)*(g*g));
    exit_flag = cosine  <= optimality_epsilon; 
    
    tic2 = timer.elapsed();

    // Line Search
    Real lambda = 1.0;
    bool do_line_search = to_line_search && ((x0-x1)*g<0);
    if(do_line_search){
      if(to_golden_search) lambda = golden_search(x0, x1, obj, line_search_iterations);
      else lambda = line_search(x0, x1, obj, line_search_iterations);
      x1 -= x0; x1 *= lambda; x1 += x0; // x1 = x0 + betamin*(x1-x0)
      f1 = obj->f(x1);
    }
    
    tic3 = timer.elapsed();
    
    Real f_ls = f1, fsp;
    Real taustar0 = 0.0;
    FOR(iter, SP_iterations_per_SOCP) {
      double tick_gs_start = timer.elapsed(), tick_gs_end;
      g = obj->g(x1);
      //iteration_report<<"1: "<<timer.elapsed() - tick_gs_start<<endl;

      FOR(a, A){  
        int u = net.arcs[a].head, v = net.arcs[a].tail;
        adjl.costs[indexadjl[u][v]] = cost_t(g[indexarcl[u][v]*K]);
      }
      //iteration_report<<"2: "<<timer.elapsed() - tick_gs_start<<endl;
        
      FOR(i, V) if(nv[i]>0) dijkstra(adjl, i, vb[i], nv[i], trace[i]);
      //iteration_report<<"3: "<<timer.elapsed() - tick_gs_start<<endl;
      
      Vector sp(x1.size(), 0);
      FOR(k, K){
        int u = net.commoflows[k].origin, v = net.commoflows[k].destination;
        Real demand = net.commoflows[k].demand;
        while(v>=0 && v!=u) {
	  sp[ indexarcl [trace[u][v]] [v] * K + k] = demand;
	  v = trace[u][v];
        }
      }
      //iteration_report<<"4: "<<timer.elapsed() - tick_gs_start<<endl;

      //Vector gsp, dx(sp); dx -= x1;
      //Real gxdx = (dx*g);

      if (taustar <=0) taubound = 1;
      else taubound = taustar * 4, sp -= x1, sp *= taubound, sp += x1;// dx *= taubound;
      
      int ii = 0;
      /*
      FOR(i, 20) {
	gsp = obj->g(sp);
	ii ++;
	if((gsp*dx)*gxdx<0) break;
	taubound *= 2;
	sp += dx;
	dx *= 2.0;
      }
      */
      //iteration_report<<"5: "<<ii<<" "<<timer.elapsed() - tick_gs_start<<endl;

      taustar = golden_search(x1, sp, obj, line_search_iterations);
      //iteration_report<<"6: "<<timer.elapsed() - tick_gs_start<<endl;

      x1 -= sp; x1 *= (1-taustar); x1 += sp;
      taustar *= taubound;
      f1 = obj->f(x1);

      //iteration_report<<"7: "<<timer.elapsed() - tick_gs_start<<endl;
      if(iter == 0) taustar0 = taustar;
    }
    
    // Timing and Reporting
    tic4 = timer.elapsed();
    iteration_report << left << setw(4)  << iteration
		     << right << setw(7)  << count
		     << right << setw(8) << setprecision(3) << fixed << tic4-tic1
		     << right << setw(12) << setprecision(2) << fixed << beta
		     << right << setw(5) << setprecision(3) << (do_line_search?"YES":"NO")
		     << right << setw(8) << setprecision(5) << fixed << lambda
		     << right << setw(8) << setprecision(5) << fixed << taustar0
		     << right << setw(8) << setprecision(5) << fixed << taustar
		     << right << setw(8) << setprecision(3) << fixed << tic4 - tic3
		     << right << setw(20) << setprecision(10) << scientific << f_ls
		     << right << setw(20) << setprecision(10) << scientific << f1 
		     << right << setw(10) << setprecision(1) << scientific << cosine
		     << right << setw(12) << setprecision(3) << fixed << timer.elapsed()-start
		     << endl;
    if(taustar == 0.0) taustar = 1.0;
  }
  
  // Reporting final results
  FOR(i,output_row_width) iteration_report<<"-"; iteration_report<<endl;
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;

  return x1;
}
