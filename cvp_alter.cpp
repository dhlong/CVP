#include "network.h"
#include "function.h"
#include <list>
#include <ilcplex/cplex.h>

#define FREE(p)	  \
	do{ \
		if(p) free(p); \
		p = NULL; \
	} while(0)

int A, V, K;

CPXENVptr     env = NULL;
CPXLPptr      lp  = NULL;

int numcols, numrows, numnz, numqnz;

char     *probname = NULL;  
double   *obj = NULL;
double   *rhs = NULL;
char     *sense = NULL;
int      *matbeg = NULL;
int      *matcnt = NULL;
int      *matind = NULL;
double   *matval = NULL;
double   *lb = NULL;
double   *ub = NULL;
int      *qmatbeg = NULL;
int      *qmatcnt = NULL;
int      *qmatind = NULL;
double   *qmatval = NULL;
int      status = 0;

double **z, **p;

fstream iteration_report; // globally-accessed iteration report
fstream solve_report;     // globally-accessed solve report
SettingMapper settings;   // globally-accessed setting mapper

bool check_nonnegative(const Vector &x){
	double max_deviation = 0.0;
	ITER(x, it)  updatemin(max_deviation, it.value());
	//cout<<"Non negative deviation = "<<max_deviation<<endl;
	return max_deviation > -1.0e-6;
}

bool check_conservation(const MultiCommoNetwork &net,const Vector &x){
	vector< vector<double> > netflow(K);
	double max_deviation = 0.0;
	FOR(k, K) netflow[k] = vector<double>(V, 0.0);
	
	ITER(x, itx){
		int a = itx.index()/K, k = itx.index() % K;
		netflow[k][net.arcs[a].head] -= itx.value();
		netflow[k][net.arcs[a].tail] += itx.value();
	}
	
	FOR(k, K) FOR(v, V)
		if(v == net.commoflows[k].origin){
			if(fabs(netflow[k][v]+net.commoflows[k].demand)>1e-6)
				cout<<"Violation at commo "<<k<<" node "<<v<<" (origin)."<<endl;
			updatemax(max_deviation, fabs(netflow[k][v]+net.commoflows[k].demand));
		}
		else if(v == net.commoflows[k].destination){
			if(fabs(netflow[k][v]-net.commoflows[k].demand)>1e-6)
				cout<<"Violation at commo "<<k<<" node "<<v<<" (destin)."<<endl;
			updatemax(max_deviation, fabs(netflow[k][v]-net.commoflows[k].demand));
		}
		else{
			if(fabs(netflow[k][v])>1e-6)
				cout<<"Violation at commo "<<k<<" node "<<v<<" (interm)."<<endl;
			updatemax(max_deviation, fabs(netflow[k][v]));
		}
	
	//cout<<"Max deviation from conservation = "<<max_deviation<<endl;
	
	return max_deviation < 1e-6;
}

bool check_capacity(const MultiCommoNetwork &net, const Vector &x){
	vector<double> flows(A, 0.0);
	ITER(x, itx) flows[itx.index()/K] += itx.value();
	FOR(a,A) if(flows[a] > net.arcs[a].cap){
		cout<<"Violate of cap constraint at arc "<<a<<endl;
		return false;
	}
	return true;
}

void allocate() {
	probname = (char   *) malloc (90      * sizeof(char));   
	obj      = (double *) malloc (numcols * sizeof(double));
	rhs      = (double *) malloc (numrows * sizeof(double));
	sense    = (char   *) malloc (numrows * sizeof(char)); 
	matbeg   = (int    *) malloc (numcols * sizeof(int));   
	matcnt   = (int    *) malloc (numcols * sizeof(int));   
	matind   = (int    *) malloc (numnz   * sizeof(int));   
	matval   = (double *) malloc (numnz   * sizeof(double));
	lb       = (double *) malloc (numcols * sizeof(double));
	ub       = (double *) malloc (numcols * sizeof(double));
	qmatbeg  = (int    *) malloc (numcols * sizeof(int)); 
	qmatcnt  = (int    *) malloc (numcols * sizeof(int)); 
	qmatind  = (int    *) malloc (numqnz  * sizeof(int)); 
	qmatval  = (double *) malloc (numqnz  * sizeof(double)); 
}

void init(const MultiCommoNetwork &net){
	V = net.getNVertex();
	A = net.arcs.size();
	K = net.commoflows.size();

	numcols = A; numrows = V; numnz = 2*A; numqnz = A;
	allocate();
  
	strcpy (probname, "example");

	FOR(v, V) 
		sense[v] = 'E', 
		rhs[v] = 0.0;

	FOR(a, A) 
		obj[a]        = 0.0,
		lb[a]         = 0.0, 
		ub[a]         = CPX_INFBOUND,
		matcnt[a]     = 2, 
		matbeg[a]     = 2*a,
		matind[2*a]   = net.arcs[a].head,
		matind[2*a+1] = net.arcs[a].tail,
		matval[2*a]   = -1.0,
		matval[2*a+1] = 1.0,
		qmatbeg[a]    = a,
		qmatind[a]    = a, 
		qmatcnt[a]    = 1, 
		qmatval[a]    = 2.0;

	assert((env = CPXopenCPLEX (&status)) != NULL);
	assert((lp  = CPXcreateprob (env, &status, probname)) != NULL);
    
	assert(!CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF));
	assert(!CPXsetintparam (env, CPX_PARAM_THREADS, 1));

	assert(!CPXcopylp (env, lp, numcols, numrows, CPX_MIN, obj, rhs, 
	                   sense, matbeg, matcnt, matind, matval,
	                   lb, ub, NULL));
  
	assert(!CPXcopyquad (env, lp, qmatbeg, qmatcnt, qmatind, qmatval));

	z = (double**) malloc(K*sizeof(double*));
	p = (double**) malloc(K*sizeof(double*));
	FOR(k, K) z[k] = (double*) malloc(A*sizeof(double));
	FOR(k, K) p[k] = (double*) malloc(A*sizeof(double));
}

void socp(const MultiCommoNetwork &net, const Vector &x0, const Vector &g, Real beta, Vector &p_){
	typedef pair<int, Real> PAIRIR;
	typedef list<PAIRIR> LPAIRIR;
	int *collist = (int*) malloc(A*sizeof(int));
	double *z = (double*) malloc(A*sizeof(double));
	double *p = (double*) malloc(A*sizeof(double));
	vector<LPAIRIR> x(K);
	vector<LPAIRIR> tmp(A);
	double rhsval[2];
	int rhsind[2];

	FOR(a, A) z[a] = 0.0;
	ITER(g, itg) z[itg.index()] = 2*beta*itg.value();

	ITER(x0, itx0)
		x[itx0.index()%K].push_back(make_pair(itx0.index()/K, itx0.value()));

	//ITER(z_, itz) z[itz.index()%K][itz.index()/K] = -2*itz.value();

	FOR(a, A) collist[a] = a;
  
	FOR(k, K){
		for(LPAIRIR::iterator it = x[k].begin(); it != x[k].end(); ++it)
			z[(*it).first] -= 2*(*it).second;

		assert(!CPXchgobj(env, lp, A, collist, z));
    
		rhsind[0] =   net.commoflows[k].origin;
		rhsind[1] =   net.commoflows[k].destination;
		rhsval[0] = - net.commoflows[k].demand;
		rhsval[1] =   net.commoflows[k].demand;
		assert(!CPXchgrhs(env, lp, 2, rhsind, rhsval));

		assert(!CPXqpopt (env, lp));
		assert(!CPXgetx(env, lp, p, 0, A-1));

		rhsval[0] = rhsval[1] = 0.0;
		assert(!CPXchgrhs(env, lp, 2, rhsind, rhsval));    

		for(LPAIRIR::iterator it = x[k].begin(); it != x[k].end(); ++it)
			z[(*it).first] += 2*(*it).second;

		FOR(a, A) if(p[a] > 1e-10) tmp[a].push_back(make_pair(k, p[a]));
	}
  
	p_ = Vector(A*K);
	FOR(a, A)
		for(LPAIRIR::iterator it = tmp[a].begin(); it != tmp[a].end(); ++it)
			p_.insert(a*K+(*it).first) = (*it).second;

	FREE(z);
	FREE(p);
	FREE(collist);
}

void release(){
	FOR(k, K) FREE(z[k]);
	FOR(k, K) FREE(p[k]);
	FREE(z);
	FREE(p);
	FREE(probname);
	//cout<<"Mark 2.1"<<endl;
	//FREE(obj);
	//cout<<"Mark 2.2"<<endl;
	FREE(rhs);
	FREE(sense);
	FREE(matbeg);
	FREE(matcnt);
	FREE(matind);
	FREE(matval);
	FREE(lb);
	FREE(ub);
	FREE(qmatbeg);
	FREE(qmatcnt);
	FREE(qmatind);
	FREE(qmatval);
	CPXcloseCPLEX(&env);
}

void allocate2() {
	probname = (char   *) malloc (90      * sizeof(char));   
	obj      = (double *) malloc (numcols * sizeof(double));
	rhs      = (double *) malloc (numrows * sizeof(double));
	sense    = (char   *) malloc (numrows * sizeof(char)); 
	matbeg   = (int    *) malloc (numcols * sizeof(int));   
	matcnt   = (int    *) malloc (numcols * sizeof(int));   
	matind   = (int    *) malloc (numnz   * sizeof(int));   
	matval   = (double *) malloc (numnz   * sizeof(double));
	lb       = (double *) malloc (numcols * sizeof(double));
	ub       = (double *) malloc (numcols * sizeof(double));
}

Vector init2(const MultiCommoNetwork &net){
	V = net.getNVertex();
	A = net.arcs.size();
	K = net.commoflows.size();

	numcols = A*K+1; numrows = V*K+A; numnz = 3*A*K+2*K;
	allocate2();
  
	strcpy (probname, "Concurrent flow");

	FOR(i, V*K)
		sense[i] = 'E',
		rhs[i] = 0.0;
	FOR(a, A)
		sense[a+V*K] = 'L',
		rhs[a+V*K] = net.arcs[a].cap;
    
	FOR(i, A*K)
		obj[i]        = 0.0,
		lb[i]         = 0.0, 
		ub[i]         = CPX_INFBOUND,
		matcnt[i]     = 3, 
		matbeg[i]     = 3*i,
		matind[3*i]   = net.arcs[i/K].head*K + i%K,
		matind[3*i+1] = net.arcs[i/K].tail*K + i%K,
		matind[3*i+2] = V*K+i/K,
		matval[3*i]   = -1.0,
		matval[3*i+1] = 1.0,
		matval[3*i+2] = 1.0;

	int i = A*K;
	obj[i] = 1.0;
	lb[i] = 0.0;
	ub[i] = CPX_INFBOUND;
	matcnt[i] = 2*K;
	matbeg[i] = 3*i;
	FOR(k, K)
		matind[3*i + 2*k]     = net.commoflows[k].origin*K + k,
		matind[3*i + 2*k + 1] = net.commoflows[k].destination*K + k,
		matval[3*i + 2*k]     = net.commoflows[k].demand,
		matval[3*i + 2*k + 1] = -net.commoflows[k].demand;

	assert((env = CPXopenCPLEX (&status)) != NULL);
	assert((lp  = CPXcreateprob (env, &status, probname)) != NULL);
    
	assert(!CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_ON));
	assert(!CPXsetintparam (env, CPX_PARAM_THREADS, 1));

	assert(!CPXcopylp (env, lp, numcols, numrows, CPX_MAX, obj, rhs, 
	                   sense, matbeg, matcnt, matind, matval,
	                   lb, ub, NULL));
  
	double* p = (double *) malloc(A*K*sizeof(double ));
	double lambda = 1.0;

	assert(!CPXbaropt (env, lp));  
	CPXwriteprob(env, lp, "concur.rlp", NULL);
	int statind = CPXgetstat (env, lp);
	char buffer[100];
	char* ptr = CPXgetstatstring (env, statind, buffer);
	cout<<ptr<<endl;

	assert(!CPXgetx(env, lp, p, 0, A*K-1));
	assert(!CPXgetx(env, lp, &lambda, A*K, A*K));
	assert(lambda >= 1.0);

	Vector x0(A*K);
	FOR(i, A*K) if(p[i] > 1e-8) x0.insert(i) = p[i]/lambda;
	FREE(p);
	return x0;
}

void release2(){
	FREE(probname);
	FREE(obj);
	FREE(rhs);
	FREE(sense);
	FREE(matbeg);
	FREE(matcnt);
	FREE(matind);
	FREE(matval);
	FREE(lb);
	FREE(ub);
	CPXcloseCPLEX(&env);
	env = NULL; lp = NULL;
}

void solve(const MultiCommoNetwork &net, ReducableFunction *obj){
	Real beta;
	int count = 0;
	Timer* timer = new CPUTimer;
	bool exit_flag = false;
	Vector x0(A*K), x1(A*K), sp(A*K);

	timer->record();

	// Initialisation by solving the intial network shortest paths
	ShortestPathOracle DA(net);
	FOR(a, A)
		DA.set_cost(net.arcs[a].head, 
		            net.arcs[a].tail, 
		            cost_t(net.arcs[a].cost));

	cout<<"Before dijkstra"<<endl;
	DA.get_flows(x1);
	cout<<"After dijkstra"<<endl;


	Real f0, f1 = obj->f(x1); 
	beta = settings.getr("initial beta") * sqrt(x1.dot(x1));
  
	// Timing and Reporting
	timer->record();
	iteration_report << "Initialisation: time = " << timer->elapsed() 
	                 << "s; obj = " << f1
	                 << endl;

	// format and header of the iteration report
	TableReport tr("%-4d  %7d     %8.3f   %12.2f"
	               "%5.3s"
	               "%8.5f   %8.5f   %8.5f"
	               "%8.3f %20.10e %20.10e %10.1e  %12.3f");

	tr.print_header(&iteration_report,
	                "Iter", "#solve",  "t_total",   "beta",
	                "ls?", 
	                "lambda*", "tau*0",     "tau*n",
	                "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");
	tr.print_header(&cout,
	                "Iter", "#solve",  "t_total",   "beta",
	                "ls?", 
	                "lambda*", "tau*0",     "tau*n",
	                "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");

	Vector g, z, y0;
	Vector y1(obj->reduced_variable(x1));
	Function *robj = obj->reduced_function();

	// Loops
	Real taubound = -1, taustar = 0.5/5, taustar0 = 1.0;
	for(int iteration = 1; !exit_flag; iteration++) {
		timer->record(); // Timming

		if(settings.getb("to reset beta")) 
			beta = settings.getr("initial beta") * sqrt(x1.dot(x1));

		// Previous best solution
		x0 = x1; f0 = f1; y0 = y1;
    
		// normalized gradient
		Vector g0 = robj->g(y0);
		g0 *= (1/sqrt(g0.dot(g0)*K));

		for(count = 1;; count++) {
			//z = g; z *= (-beta); z += x0;
			//z = x0 - beta*g;
			socp(net, x0, g0, beta, x1);
			f1 = obj->f(x1);
			if(f1 < f0) break;
			//if((x0-x1).dot(obj->g(x1)) < 0) break;
			beta *= settings.getr("beta down factor");
		}

		y1 = obj->reduced_variable(x1);
		Vector g1 = robj->g(y1); // g is now gradient at x1

		// Optimality check
		Vector dx = x0 - x1, dy = y0 - y1;
		Real g1dx = dy.dot(g1), g1g1 = K*g1.dot(g1), dxdx = dx.dot(dx);
		Real g0dx = dy.dot(g0), g0g1 = K*g0.dot(g1), g0g0 = K*g0.dot(g0);
		Real cosine = 1 + ( (g1dx - beta*g0g1) /
		                    sqrt((dxdx -2*beta*g0dx + beta*beta*g0g0)*g1g1) );
		exit_flag = cosine  <= settings.getr("optimality epsilon"); 
    
		timer->record(); // for timing


		// Line Search
		Real lambda = 1.0;
		bool do_line_search = ( settings.getb("to do line search") && 
		                        g1dx < 0 );
		if(do_line_search){
			lambda = section_search (y0, y1, robj, 
			                         settings.geti("line search iterations"));
			x1 = x0 + lambda*(x1-x0);
			y1 = y0 + lambda*(y1-y0);
			//x1 -= x0; x1 *= lambda; x1 += x0;
			f1 = robj->f(y1);
		}
    
		timer->record(); // for timing
		Real f_ls = f1; // for reporting

		double tau = taustar0*20;
		updatemin(tau, 1.0);
		FOR(iter, settings.geti("SP iterations per SOCP")) {
			g1 = robj->g(y1);
			DA.reset_cost();
			ITER(g1, itg1)
				DA.set_cost ( net.arcs[itg1.index()].head,
				              net.arcs[itg1.index()].tail,
				              cost_t(itg1.value()));
			DA.get_flows(sp);
			Vector ysp(obj->reduced_variable(sp));

			taustar = section_search(y1, ysp, robj, 
			                         settings.geti("line search iterations"),
			                         false, tau*(1-PHI), tau*PHI);

			x1 += taustar * (sp  - x1);  //x1 -= sp;  x1 *= (1-taustar); x1 += sp; 
			y1 += taustar * (ysp - y1); //y1 -= ysp; y1 *= (1-taustar); y1 += ysp;

			if(iter == 0) taustar0 = taustar;
		}
    
		f1 = robj->f(y1);
		// Timing and Reporting
		timer->record();
		tr.print_row (&iteration_report,
		              iteration, count, timer->elapsed(-1,-4), beta,
		              (do_line_search?"YES":"NO"), 
		              lambda, taustar0, taustar,
		              timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));

		tr.print_row (&cout,
		              iteration, count, timer->elapsed(-1,-4), beta,
		              (do_line_search?"YES":"NO"), 
		              lambda, taustar0, taustar,
		              timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));
		if(taustar == 0.0) taustar = 1.0; 
      
	}

	// Check feasibility of the final solution
	assert(check_conservation(net,x1));
	assert(check_nonnegative(x1));

	// Reporting final results
	tr.print_line(iteration_report);
	tr.print_line(cout);
	iteration_report << "Optimal objective = " 
	                 << scientific << setprecision(12) << obj->f(x1)
	                 << endl;
	delete timer;
}

void solve_KL(const MultiCommoNetwork &net){
	KleinrockFunction *obj = new KleinrockFunction(net);
	Real beta;
	int count = 0;
	Timer* timer = new CPUTimer;
	bool exit_flag = false;
	Vector x0(A*K), x1(A*K), sp(A*K);
	ShortestPathOracle DA(net);

	timer->record();

	// Initialisation by solving the concurrent flow problem
	x1 = init2(net);
	release2();

	Real f0, f1 = obj->f(x1); 
	beta = settings.getr("initial beta") * sqrt(x1.dot(x1));
  
	// Timing and Reporting
	timer->record();
	iteration_report << "Initialisation: time = " << timer->elapsed()
	                 << "s; obj = " << f1 
	                 << endl;

	// format and header of the iteration report
	TableReport tr("%-4d  %7d     %8.3f   %12.2f"
	               "%5.3s"
	               "%8.5f   %8.5f   %8.5f"
	               "%8.3f %20.10e %20.10e %10.1e  %12.3f");

	tr.print_header(&iteration_report,
	                "Iter", "#solve",  "t_total",   "beta",
	                "ls?", 
	                "lambda*", "tau*0",     "tau*n",
	                "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");
	tr.print_header(&cout,
	                "Iter", "#solve",  "t_total",   "beta",
	                "ls?", 
	                "lambda*", "tau*0",     "tau*n",
	                "t_SP", "obj_ls",  "obj_final", "cosine", "t_elapsed");

	Vector g, z, y0;
	Vector y1(obj->reduced_variable(x1));
	Function *robj = obj->reduced_function();

	init(net);

	// Loops
	Real taubound = -1, taustar = 0.5/5, taustar0 = 1.0;
	for(int iteration = 1; !exit_flag; iteration++) {
		if(settings.getb("to reset beta")) 
			beta = settings.getr("initial beta") * sqrt(x1.dot(x1));

		// Timing and Reporting
		timer->record();

		// Previous best solution
		x0 = x1; f0 = f1; y0 = y1;
    
		// normalized gradient
		g = obj->g(x0);
		g *= (1/sqrt(g.dot(g)));

		for(count = 1;;count++) {
			//z = g; z *= (-beta); z += x0;
			z = x0 - beta*g;
			//socp(net, z, x1);
			if(check_capacity(net, x1)){
				f1 = obj->f(x1);
				if(f1 < f0) break;
				if((x0-x1).dot(obj->g(x1)) < 0) break;
			}
			beta *= settings.getr("beta down factor");
		}

		// Optimality check
		g = obj->g(x1); // g is now gradient at x1
		z -= x1; // z is now z - x1
		Real cosine = 1 + (z.dot(g))/sqrt((z.dot(z))*(g.dot(g)));
		exit_flag = cosine  <= settings.getr("optimality epsilon"); 
    
		timer->record(); // for timing

		y1 = obj->reduced_variable(x1);

		// Line Search
		Real lambda = 1.0;
		bool do_line_search = ( settings.getb("to do line search") && 
		                        g.dot(x0-x1)<0 );
		if(do_line_search){
			lambda = section_search (y0, y1, robj, 
			                         settings.geti("line search iterations"));
			x1 = x0 + lambda*(x1-x0);
			y1 = y0 + lambda*(y1-y0);
			//x1 -= x0; x1 *= lambda; x1 += x0;
			f1 = robj->f(y1);
		}
    
		timer->record(); // for timing
    
		Real f_ls = f1; // for reporting

		double tau = taustar0*20;
		updatemin(tau, 1.0);
		FOR(iter, settings.geti("SP iterations per SOCP")) {
			Vector gy(robj->g(y1));
			DA.reset_cost();
			ITER(gy, itgy)
				DA.set_cost(net.arcs[itgy.index()].head,
				            net.arcs[itgy.index()].tail,
				            cost_t(itgy.value()));
			DA.get_flows(sp);

			Vector ysp(obj->reduced_variable(sp));
			taustar = section_search(y1, ysp, robj, 
			                         settings.geti("line search iterations"),
			                         false, tau*(1-PHI), tau*PHI);

			x1 += taustar*(sp-x1);  //x1 -= sp;  x1 *= (1-taustar); x1 += sp; 
			y1 += taustar*(ysp-y1); //y1 -= ysp; y1 *= (1-taustar); y1 += ysp;

			if(iter == 0) taustar0 = taustar;
		}
    
		f1 = robj->f(y1);

		// Timing and Reporting
		timer->record();
		tr.print_row(&iteration_report,
		             iteration, count, timer->elapsed(-1,-4), beta,
		             (do_line_search?"YES":"NO"), 
		             lambda, taustar0, taustar,
		             timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));
		tr.print_row(&cout,
		             iteration, count, timer->elapsed(-1,-4), beta,
		             (do_line_search?"YES":"NO"), 
		             lambda, taustar0, taustar,
		             timer->elapsed(), f_ls, f1, cosine, timer->elapsed(0,-1));
		if(taustar == 0.0) taustar = 1.0;     
	}
  
	// check feasibility of the final solution
	assert(check_conservation(net,x1));
	assert(check_nonnegative(x1));
	assert(check_capacity(net,x1));

	// Reporting final results
	tr.print_line(iteration_report);
	tr.print_line(cout);
	iteration_report<<"Optimal objective: "<<scientific<<setprecision(12)<<obj->f(x1)<<endl;
	delete timer;
	delete obj;
}

int main(){
	Timer *timer = new CPUTimer;
  
	// get the system date and time
	// append date and time to the end of the report file name
	// so that the reports won't be overwritten after each run
	time_t rawtime;
	struct tm * timeinfo;
	char timestr[100];
	time(&rawtime );
	timeinfo = localtime (&rawtime);
	sprintf(timestr,"%04d%02d%02d%02d%02d",
	        timeinfo->tm_year+1900, 
	        timeinfo->tm_mon+1, 
	        timeinfo->tm_mday, 
	        timeinfo->tm_hour, 
	        timeinfo->tm_min);

	// start the timer
	timer->record();

	// Reading settings from file "CVP.ini"
	string ininame = "CVP.ini";
	fstream fini(ininame.c_str(), fstream::in);
	settings.read(fini);
	fini.close();
  
	// generate input and report file names
	string inputname  = settings.gets("input file");
	string outputname = inputname + "_iterations_" + string(timestr) + ".txt";
	string timername  = inputname + "_solves_" + string(timestr) + ".txt";
	string pajekname  = inputname + "pajek.net";
 
	// prepare files for reporting
	iteration_report.open(outputname.c_str(), fstream::out);
	solve_report.open(timername.c_str(), fstream::out);
  
	iteration_report<<"Solving problem \""<<inputname<<"\""<<endl<<endl;

	// Report settings
	iteration_report<<"Settings:"<<endl;
	settings.report(iteration_report);
  
	// check format of the input based on the filename given
	// if file name contains "grid" or "planar" --> genflot format
	// else tntp format
	FileFormat format = GENFLOT;
	if ((strstr(inputname.c_str(),"grid")   == NULL ) &&
	    (strstr(inputname.c_str(),"planar") == NULL))
		format = TNTP;
  
	MultiCommoNetwork net(inputname.c_str(), format);
  
	cout<<"Solving"<<endl;
	timer->record();

	if(settings.gets("Function") == "bpr"){
		init(net);
		BPRFunction *obj = new BPRFunction(net);
		solve(net,obj);
		release();
		delete obj;
	}
	else{
		solve_KL(net);
		release();
	}
	timer->record();

	iteration_report << "Toal CPU time = "
	                 << left << setw(10) << setprecision(5) << fixed
	                 << timer->elapsed(0,-1) << "s"
	                 << endl;

	iteration_report << "Peak memory   = " 
	                 << left << setw(10) << setprecision(2) << fixed
	                 << memory_usage()/1e3 << "kB" <<endl;

	iteration_report.close();
	solve_report.close();

	delete timer;
	delete obj;
	return 0;
}

