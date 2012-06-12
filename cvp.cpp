#include "cvp.h"
#include <set>
#include <list>
#include "dijkstra.h"

ILOSTLBEGIN

// report files, defined in main.cpp
extern fstream iteration_report;
extern fstream solve_report;
extern SettingMapper settings;

#include <Eigen/Dense>

MatrixXd projection_matrix(const Graph &net){
  int V = net.getNVertex(), A = net.arcs.size();
  MatrixXd m = MatrixXd::Zero(V-1,A);
  FOR(a, A){
    if(net.arcs[a].head != V-1) m(net.arcs[a].head, a) = 1.0;
    if(net.arcs[a].tail != V-1) m(net.arcs[a].tail, a) = -1.0;
  }
  
  MatrixXd M = (-1.0)*m.transpose()*(m*m.transpose()).inverse()*m;
  FOR(a, A) M(a, a) += 1.0;
  
  return M;
}

Vector projection(const MultiCommoNetwork &net, const Vector &x, const MatrixXd &M){
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  
  MatrixXd X = MatrixXd::Zero(A,K);
  ITER(x, itx) X(itx.index()/K, itx.index()%K) = itx.value();
  
  MatrixXd P = M*X;

  Vector p(A*K);
  FOR(a, A) FOR(k, K) p.insert(a*K+k) = P(a,k);

  return p;
}

bool check_nonnegative(const Vector &x){
  double max_deviation = 0.0;
  ITER(x, it)  updatemin(max_deviation, it.value());
  //cout<<"Non negative deviation = "<<max_deviation<<endl;
  return max_deviation > -1.0e-6;
}

bool check_conservation(const MultiCommoNetwork &net, const Vector &x){
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  vector< vector<double> > netflow(K);
  double max_deviation = 0.0;
  FOR(k, K) netflow[k] = vector<double>(V, 0.0);

  ITER(x, itx){
    int a = itx.index()/K, k = itx.index() % K;
    netflow[k][net.arcs[a].head] -= itx.value();
    netflow[k][net.arcs[a].tail] += itx.value();
  }

  FOR(k, K) FOR(v, V)
    if(v == net.commoflows[k].origin)
      updatemax(max_deviation, fabs(netflow[k][v]+net.commoflows[k].demand));
    else if(v == net.commoflows[k].destination)
      updatemax(max_deviation, fabs(netflow[k][v]-net.commoflows[k].demand));
    else updatemax(max_deviation, fabs(netflow[k][v]));

  //cout<<"Max deviation from conservation = "<<max_deviation<<endl;

  return max_deviation < 1e-6;
}

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
  proxy_obj(NULL) { }

// constructor from a potin to a objective function
CVP::CVP(Function *obj_) : 
  env(), 
  variables(env), 
  constraints(env),
  obj(obj_), 
  cplex(NULL), 
  proxy(NULL),
  proxy_obj(NULL) { }

// default initial proxy objective
// projection the 0-point on the feasible set
IloObjective CVP::initial_proxy_obj(){
  Vector x(variables.getSize());
  return quad_proxy_obj(x);
}

// default initial solution by solving the initial proxy objective
// can be (virtually) redefined in the inherited classes
Vector CVP::initial_solution(){
  return solve(this->initial_proxy_obj());
}

// return the proxy objective min ||x - z||
IloObjective CVP::quad_proxy_obj(const Vector &z){
  IloExpr sum(env);
  ITER(z, itz) sum += (2*itz.value())*variables[itz.index()];
  
  IloObjective obj = IloMinimize(env,
				 IloScalProd(variables, variables)
				 + (z.dot(z))
				 - sum);
  sum.end();
  return obj;
}

// return the proxy objective min gradient(x0) * x
IloObjective CVP::linear_proxy_obj(const Vector &x0){
  IloExpr sum(env);
  Vector g(obj->g(x0));
  g *= (1.0/sqrt(g.dot(g)));
  ITER(g, itg) sum += itg.value()*variables[itg.index()];
  IloObjective obj = IloMinimize(env, sum);
  sum.end();
  return obj;
}

int mark = true;

bool CVP::is_optimal(const Vector &x, const Vector &z){
  Vector d(obj->g(x)), y(z);
  y -= x;
  return 1+(d.dot(y))/sqrt((y.dot(y))*(d.dot(d))) <= settings.getr("optimality epsilon");
}

Vector CVP::solve(const IloObjective &iloobj){
  Timer *timer = new CPUTimer;
  TableReport tr("%15.5f%15.5f%15.5f%15.5f%15.5f%20.10e");
  if(proxy_obj == NULL) // if this is first solve
    tr.print_header(solve_report,"t1", "t2", "t3", "t4", "t5", "obj");

  cout<<"New Solve"<<endl;
  timer->record();
  
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
  timer->record(); //timing

  // add new proxy objective to the proxy model
  proxy->add(*proxy_obj);
  timer->record(); // timing
  
  // If cplex has not been created, create it
  if (cplex == NULL){
    cplex = new IloCplex(*proxy);
    //cplex->exportModel("mypro2.rlp");
  }
  timer->record(); // timing

  cplex->setParam(IloCplex::Threads, 1);
  cplex->setOut(env.getNullStream());
  // solve
  if(!cplex->solve()) {
    cout<<cplex->getStatus()<<endl;
    solve_report<<cplex->getStatus()<<endl;
    return Vector();
  }

  timer->record(); // timing

  // obtain the optimal solution
  IloNumArray vals(env,variables.getSize());
  cplex->getValues(vals, variables);

  Vector x(vals.getSize());
  FOR(i, vals.getSize()) if(vals[i] > 1e-10) x.coeffRef(i) = Real(vals[i]);

  vals.end();

  // calculate new real objective (for reporting)
  Real obj_value = (obj==NULL ? -1 : obj->f(x));
  timer->record();
  
  // reporting to solve report
  tr.print_row(solve_report,
	       timer->elapsed(-6,-5), timer->elapsed(-5,-4),
	       timer->elapsed(-4,-3), timer->elapsed(-3,-2),
	       timer->elapsed(-2,-1), obj_value);
  delete timer;
  return x;
}
  
CVP::~CVP(){
  if(cplex) delete cplex;
  if(proxy) delete proxy;
  if(proxy_obj) delete proxy_obj;
  if(obj) delete obj;
}

Vector CVP::phase1(Vector *init){
  Real pre_cost, cur_cost;
  Vector x, pre_x;
  Timer *timer = new CPUTimer;
  int iteration = 0;
  double tic1, tic2;

  // Initialisation
  if(init!=NULL && init->size() == variables.getSize()) x = *init;
  else x = solve(this->initial_proxy_obj());
  cur_cost = obj->f(x);

  do{
    iteration++;
    timer->record();    
    pre_x = x, pre_cost = cur_cost;
    cur_cost = obj->f(x = solve(linear_proxy_obj(pre_x)));
    timer->record();
    
    iteration_report<<"Phase 1\t"<<iteration<<"\t"<<"\t"<<timer->elapsed(-1,-2)<<endl;
    if(x==pre_x){
      delete timer;
      return x;
    }
  } while(cur_cost<=pre_cost);
  
  delete timer;
  return Vector();
}

Vector CVP::phase2(Vector *init){
  Real beta;
  int count = 0, iteration = 0;
  Timer *timer = new CPUTimer();
  bool exit_flag = false;

  timer->record();

  // Initialisation by solving initial proxy objective
  Vector x0, x1 = solve(this->initial_proxy_obj());
  Real f0, f1 = obj->f(x1); 

  // Timing and Reporting
  timer->record();
  iteration_report<<"Initialisation: time = "<<timer->elapsed()<<"s; obj = "<<f1<<endl;
  timer->record();

  beta = settings.getr("initial beta") * sqrt(x0.dot(x0)); // added by Hieu

  // header row of the iteration report
  TableReport tr("%-6d%-8d%8.4f%10.5f%6s%8.4f%8.3f%10.5e%17.8e");
  tr.print_header(iteration_report,
		  "Iter", "#solve", "time", "beta", "ls?",
		  "time_ls", "lambda*", "cosine", "obj");
  
  // Loops
  Vector g, z;
  while(!exit_flag) {
    if(settings.getb("to reset beta")) 
      beta = settings.getr("initial beta") * sqrt(x1.dot(x1));

    // Timing and Reporting
    ++iteration; count = 0; timer->record();

    // Previous best solution
    x0 = x1; f0 = f1;
    
    // Reduce beta and do projection until improvement
    g = obj->g(x0);
    Real normg = sqrt(g.dot(g)); // added by Hieu

    for(;;) {
      z = g; z *= (-beta/normg); z += x0; // z = x0 - beta*g
      x1 = solve(quad_proxy_obj(z));
      f1 = obj->f(x1);
      count++; // counting number of solves (for reporting purpose)
      if(f1 < f0) break;
      beta *= settings.getr("beta down factor");
    }

    //Real normdx_before_ls = sqrt((x1-x0)*(x1-x0)); // added by Hieu

    // Optimality check
    g = obj->g(x1); // g is now gradient at x1
    z -= x1; // z is now z - x1
    Real cosine = 1 + (z.dot(g))/sqrt((z.dot(z))*(g.dot(g))); //added by Hieu
    exit_flag = cosine <= settings.getr("optimality epsilon");
    
    timer->record();

    // Line Search
    Real lambda = 1.0;
    bool do_line_search = settings.getb("to do line search") && (g.dot(x0-x1)<0);
    if(do_line_search){
      if(settings.getb("to do golden search")) 
	lambda = section_search(x0, x1, obj, settings.geti("line search iterations"));
      else 
	lambda = line_search(x0, x1, obj, settings.geti("line search iterations"));
      x1 -= x0; x1 *= lambda; x1 += x0; // x1 = x0 + betamin*(x1-x0)
      f1 = obj->f(x1);
    }

    //Real normdx_after_ls = sqrt((x1-x0)*(x1-x0)); // added by Hieu
    
    // Timing and Reporting
    timer->record();
    tr.print_row(iteration_report, 
		 iteration, count, timer->elapsed(-1,-3), beta, 
		 do_line_search ? "YES" : "NO",
		 timer->elapsed(-2, -1), lambda, cosine, f1);
  }
  
  // Reporting final results
  tr.print_line(iteration_report);
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;
  delete timer;
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
    ssname<<"flow["<<net.arcs[i].head<<">"<<net.arcs[i].tail<<"]";

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

  node.end();
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
  IloObjective obj = IloMinimize(env, expr);
  expr.end();
  return obj;
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
CVP_MCNF::CVP_MCNF(Function *obj_, const MultiCommoNetwork &net_)
  : CVP(obj_), net(net_) {
}

Vector CVP_MCNF::optimize(){
  if (!settings.getb("to do SOCP")){
    env.end();
    return solve_by_dijkstra();
  }
  this->generate_network_constraints();  
  if (!settings.getb("to do shortest path")) return phase2();
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
      if(settings.getb("to include capacity constraints")) 
	sum += variables[a*K +i];
    }

    FOR(i, K){
      // Update flow equation for each commodity at each arc with the new var
      node[i*V + net.arcs[a].head] -= variables[a*K + i];
      node[i*V + net.arcs[a].tail] += variables[a*K + i];
    }
    
    // Capacity constraints
    if(settings.getb("to include capacity constraints")) 
      constraints.add(sum <= net.arcs[a].cap);
    sum.end();
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
  IloObjective iloobj(IloMinimize(env, expr));
  expr.end();
  return iloobj;
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
    ITER(xx, itxx) x.coeffRef(itxx.index()*K+k) = itxx.value();
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
  return solve_by_dijkstra_only(net,obj,settings.geti("SP iterations"));
}

Vector CVP_MCNF::solve_by_dijkstra_and_SOCP(){
  Real beta;
  int count = 0, iteration = 0;
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  Timer* timer = new CPUTimer;
  bool exit_flag = false;
  Vector x0(A*K), x1(A*K), sp(A*K);

  timer->record();

  // Initialisation by solving the intial network shortest paths
  ShortestPathOracle DA(net);
  FOR(a, A){
    NetworkArc arc = net.arcs[a];
    DA.set_cost(arc.head, arc.tail, cost_t(arc.cost));
  }
  DA.get_flows(x1);

  Real f0, f1 = obj->f(x1); 
  beta = settings.getr("initial beta") * sqrt(x1.dot(x1));
  
  // Timing and Reporting
  timer->record();
  iteration_report<<"Initialisation: time = "<<timer->elapsed()<<"s; obj = "<<f1<<endl;

  // format and header of the iteration report
  TableReport tr("%-4d  %7d     %8.3f   %12.2f"
		 "%5.3s %5.3s "
		 "%8.5f   %8.5f   %8.5f"
		 "%8.3f %20.10e %20.10e %10.1e  %12.3f");
  tr.print_header(iteration_report,
		  "Iter", "#solve",  "t_total",   "beta",
		  "ls?",  "proj?", 
		  "lambda*", "tau*0",     "tau*n",
		  "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");

  Vector g, z, y1, y0;
  ReducableFunction *cobj = dynamic_cast<ReducableFunction*>(obj);
  Function *robj = (cobj)? cobj->reduced_function() : NULL;

  MatrixXd M = projection_matrix(net);
  bool use_analytical_projection = false;

  // Loops
  Real taubound = -1, taustar = 0.5/5;
  while(!exit_flag) {
    if(settings.getb("to reset beta")) beta = settings.getr("initial beta") * sqrt(x1.dot(x1));

    // Timing and Reporting
    ++iteration; count = 0; timer->record();

    // Previous best solution
    x0 = x1; f0 = f1; y0 = y1;
    
    // normalized gradient
    g = obj->g(x0);
    g *= (1/sqrt(g.dot(g)));

    cout<<"Objective before SOCP = "<<f1<<endl;

    for(;;) {
      z = g; z *= (-beta); z += x0; // z = x0 - beta*g

      // Try analytical solution first
      x1 = projection(net, z-x0, M) + x0;
      assert(check_conservation(net, x1));

      if(check_nonnegative(x1)){
	use_analytical_projection = true;
	break;
      }

      use_analytical_projection = false;
      x1 = solve(quad_proxy_obj(z));

      f1 = obj->f(x1);
      count++; // counting number of solves (for reporting)
      if(f1 < f0) break;
      beta *= settings.getr("beta down factor");
    }

    cout<<"Objective after SOCP = "<<f1<<endl;

    // Optimality check
    g = obj->g(x1); // g is now gradient at x1
    z -= x1; // z is now z - x1
    Real cosine = 1 + (z.dot(g))/sqrt((z.dot(z))*(g.dot(g)));
    exit_flag = cosine  <= settings.getr("optimality epsilon"); 
    
    timer->record(); // for timing

    // Line Search
    Real lambda = 1.0;
    bool do_line_search = settings.getb("to do line search") && (g.dot(x0-x1)<0);
    if(do_line_search){
      if(settings.getb("to do golden search")) 
	lambda = section_search(x0, x1, obj, settings.geti("line search iterations"));
      else lambda = line_search(x0, x1, obj, settings.geti("line search iterations"));
      x1 -= x0; x1 *= lambda; x1 += x0; // x1 = x0 + betamin*(x1-x0)
      f1 = obj->f(x1);
    }
    
    timer->record(); // for timing
    
    Real f_ls = f1, fsp, taustar0 = 0.0; // for reporting

    // If objective function is reducable, use the reduced function
    // for more efficient section search
    if(cobj){
      Vector y1(cobj->reduced_variable(x1));
      FOR(iter, settings.geti("SP iterations per SOCP")) {
	Vector gy(robj->g(y1));
	DA.reset_cost();
	ITER(gy, itgy)
	  DA.set_cost(net.arcs[itgy.index()].head,
		      net.arcs[itgy.index()].tail,
		      cost_t(itgy.value()));
	DA.get_flows(sp);
	Vector ysp(cobj->reduced_variable(sp));
	timer->record();
	taustar = section_search(y1, ysp, robj, settings.geti("line search iterations"));
	timer->record();
	x1 -= sp;  x1 *= (1-taustar); x1 += sp; 
	y1 -= ysp; y1 *= (1-taustar); y1 += ysp;
	timer->record();
	cout<<"Search time = "<<timer->elapsed(-2,-3)<<"; Cal iter "<< settings.geti("SP iterations per SOCP")<<"; Search = "<<timer->elapsed()<<endl;
      //cout<<"SP time = "<<timer->elapsed(-2,-3)<<"; Search = "<<timer->elapsed()<<endl;
	if(iter == 0) taustar0 = taustar;
      }
    }
    else 
      FOR(iter, settings.geti("SP iterations per SOCP")) {
	g = obj->g(x1);
	FOR(a, A) DA.set_cost(net.arcs[a].head, net.arcs[a].tail, cost_t(g.coeff(a*K)));
	DA.get_flows(sp);
	taustar = section_search ( x1, sp, obj, settings.geti("line search iterations"));
	x1 -= sp; x1 *= (1-taustar); x1 += sp;
	if(iter == 0) taustar0 = taustar; // for reporting
      }
    
    f1 = obj->f(x1);
    // Timing and Reporting
    timer->record();
    tr.print_row(iteration_report,
		 iteration, count, timer->elapsed(-1,-4), beta,
		 (do_line_search?"YES":"NO"), 
		 (use_analytical_projection?"YES":"NO"),
		 lambda, taustar0, taustar,
		 timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));
    if(taustar == 0.0) taustar = 1.0; 

    //assert(check_conservation(net, x1));
    //assert(check_nonnegative(x1));
  }
  
  // Reporting final results
  tr.print_line(iteration_report);
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;
  delete timer;
  return x1;
}

Vector solve_by_dijkstra_only(const MultiCommoNetwork &net, Function *obj, int iterations)
{
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  ShortestPathOracle DA(net);
  Timer *timer = new CPUTimer();
  Vector x(A*K), sp(A*K);

  timer->record();

  // Initialisation by solving the intial network shortest paths
  FOR(a, A) 
    DA.set_cost(net.arcs[a].head, 
		net.arcs[a].tail, 
		cost_t(net.arcs[a].cost));
  DA.get_flows(x);

  // header row of the iteration report
  TableReport tr("%-5d%8.4f%8.4f%8.4f%8.4f%20.10e%12.3f");
  tr.print_header(iteration_report, 
		  "Iter", "t_total", "t_SP", "t_LS", 
		  "tau", "obj", "t_elapsed");
  
  // If obj is reducable, a more efficient algorithm is used
  ReducableFunction *cobj = dynamic_cast<ReducableFunction*>(obj);
  Function * robj = NULL;
  Vector y;
  Real tau;
  if(cobj) robj = cobj->reduced_function(), y = cobj->reduced_variable(x);

  FOR(iteration, iterations) {
    cout<<"Iteration "<<iteration<<endl;
    timer->record();
    
    if(cobj){
      Vector g(robj->g(y));
      ITER(g, itg) 
	DA.set_cost(net.arcs[itg.index()].head, 
		    net.arcs[itg.index()].tail, 
		    cost_t(itg.value()));
      DA.get_flows(sp);
      timer->record();

      Vector ysp(cobj->reduced_variable(sp));      
      if(4*tau >= 1.0) tau = section_search(y, ysp, robj);
      else tau = section_search(y, ysp, robj, 
				settings.geti("line search iterations"), 
				4*tau*(1-PHI), 4*tau*PHI);
      x -= sp;  x *= (1-tau); x += sp;
      y -= ysp; y *= (1-tau); y += ysp;
    }
    else{
      Vector g(obj->g(x));
      FOR(a, A) 
	DA.set_cost(net.arcs[a].head, 
		    net.arcs[a].tail, 
		    cost_t(g.coeff(a*K)));
      DA.get_flows(sp);
      timer->record();
      
      if(4*tau >= 1.0) tau = section_search(x, sp, robj);
      else tau = section_search(x, sp, obj, 20, 4*tau*(1-PHI), 4*tau*PHI);
      x -= sp; x *= (1-tau); x += sp;
    }

    timer->record();
      
    // Timing and Reporting
    if(iteration%50 == 0)
      tr.print_row(iteration_report, 
		   iteration, timer->elapsed(-1,-3), timer->elapsed(-2,-3),
		   timer->elapsed(-2,-1), tau, obj->f(x), timer->elapsed(0,-1));
  }
  
  if(robj) delete robj;
  
  // Reporting final results
  tr.print_line(iteration_report);
  iteration_report<<"Optimal objective: "
		  <<scientific<<setprecision(12)<<obj->f(x)<<endl;
  delete timer;
  return x;  
}


CVP_MCNF_KL::CVP_MCNF_KL(const MultiCommoNetwork &n): 
  CVP_MCNF((Function *) new KleinrockFunction(n), n),
  capconstraints(env)
{
}

void CVP_MCNF_KL::generate_network_constraints(){
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
      sum += variables[a*K + i];
    }

    FOR(i, K){
      // Update flow equation for each commodity at each arc with the new var
      node[i*V + net.arcs[a].head] -= variables[a*K + i];
      node[i*V + net.arcs[a].tail] += variables[a*K + i];
    }
    
    // Capacity constraints
    capconstraints.add(sum <= net.arcs[a].cap*0.9999);
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


Vector CVP_MCNF_KL::solve_by_dijkstra_and_SOCP(){
  int count = 0, iteration = 0;
  int V = net.getNVertex(), A = net.arcs.size(), K = net.commoflows.size();
  Timer *timer = new CPUTimer;
  bool exit_flag = false;
  KleinrockFunction *kl = dynamic_cast<KleinrockFunction*>(obj); 
  assert(kl!=NULL);
  Function *rkl = kl->reduced_function();

  timer->record();
  ShortestPathOracle DA(net);

  timer->record();

  // Initialisation by solving the intial network shortest paths
  constraints.add(capconstraints);
  Vector x1 = solve(this->initial_proxy_obj()), x0(A*K), sp(A*K);
  proxy->remove(capconstraints);
  Real f0, f1 = obj->f(x1); 
  Real beta = settings.getr("initial beta") * sqrt(x1.dot(x1));
  
  // Timing and Reporting
  timer->record();
  iteration_report<<"Initialisation time = "<<timer->elapsed()<<"s; obj = "<<f1<<endl;

  // header row of the iteration report
  // header row of the iteration report
  TableReport tr("%-4d  %7d     %8.3f   %12.2f"
		 "%5.3s %8.5f   %8.5f   %8.5f"
		 "%8.3f %20.10e %20.10e %10.1e  %12.3f");
  tr.print_header(iteration_report,
		  "Iter", "#solve",  "t_total",   "beta",
		  "ls?",  "lambda*", "tau*0",     "tau*n",
		  "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");

  // Loops
  Vector g, z;
  Real taubound = -1, taustar = 0.5/5;
  while(!exit_flag) {
    // Timing and Reporting
    iteration ++; count = 0; timer->record();
    cout<<"Mark 0"<<endl;

    if(settings.getb("to reset beta")) 
      beta = settings.getr("initial beta") * sqrt(x1.dot(x1));

    // Previous best solution
    x0 = x1; f0 = f1;
    
    // normalized gradient
    g = obj->g(x0);
    g *= (1.0/sqrt(g.dot(g)));

    // reduce beta and do projection until feasibility and improvement
    for(;;) {
      z = g; z *= (-beta); z += x0; // z = x0 - beta*g
      x1 = solve(quad_proxy_obj(z));
      count++; // counting number of solves (for reporting)

      bool feasible = true;
      Vector y1(kl->reduced_variable(x1));
      ITER(y1, ity1) if (ity1.value() >= net.arcs[ity1.index()].cap){
	feasible = false;
	break;
      }

      if(feasible){
	f1 = rkl->f(y1);
	if(f1<f0) break;
      }
      beta *= settings.getr("beta down factor");
    }

    // Optimality check
    g = obj->g(x1); // g is now gradient at x1
    z -= x1;        // z is now z - x1
    Real cosine = 1 + (z.dot(g))/sqrt((z.dot(z))*(g.dot(g)));
    exit_flag = cosine  <= settings.getr("optimality epsilon"); 
    
    timer->record(); // for timing

    // Line Search
    Real lambda = 1.0;
    bool do_line_search = settings.getb("to do line search") && (g.dot(x0-x1)<0);
    if(do_line_search){
      if(settings.getb("to do golden search")) 
	lambda = section_search(x0, x1, obj, settings.geti("line search iterations"));
      else lambda = line_search(x0, x1, obj, settings.geti("line search iterations"));
      x1 -= x0; x1 *= lambda; x1 += x0; // x1 = x0 + betamin*(x1-x0)
      f1 = obj->f(x1);
    }
    
    timer->record(); // for timing
    
    Real f_ls = f1, fsp, taustar0 = 0.0; // for reporting
    Vector y1(kl->reduced_variable(x1));
    FOR(iter, settings.geti("SP iterations per SOCP")) {
      Vector gy(rkl->g(y1));
      ITER(gy, itgy)
	DA.set_cost(net.arcs[itgy.index()].head, 
		    net.arcs[itgy.index()].tail, 
		    cost_t(itgy.value()));
      DA.get_flows(sp);
      Vector ysp(kl->reduced_variable(sp));

      // feasibility search
      Real alpha = 1.0;
      ITER(ysp, itysp){
	int a = itysp.index();
	Real yspa = itysp.value(), y1a = y1.coeff(a);
	if(yspa > y1a) 
	  updatemin(alpha, (net.arcs[a].cap*0.9999 - y1a)/(yspa-y1a));
      }
      assert(alpha>0.0); // Hopefully we can find a new feasible solution
      sp  -= x1; sp  *= alpha; sp  += x1;
      ysp -= y1; ysp *= alpha; ysp += y1;
      
      taustar = section_search ( y1, ysp, rkl, settings.geti("line search iterations"));
      x1 -= sp;  x1 *= (1-taustar); x1 += sp; 
      y1 -= ysp; y1 *= (1-taustar); y1 += ysp; 
      f1 = rkl->f(y1);

      if(iter == 0) taustar0 = taustar; // for reporting
    }
    
    f1 = obj->f(x1);

    // Timing and Reporting
    timer->record();
    tr.print_row(iteration_report,
		 iteration, count, timer->elapsed(-1,-4), beta,
		 (do_line_search?"YES":"NO"), lambda, taustar0, taustar,
		 timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));
    if(taustar == 0.0) taustar = 1.0; 
  }
  
  // Reporting final results
  tr.print_line(iteration_report);
  iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;
  delete timer;
  return x1;
}
