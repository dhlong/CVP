#include "function.h"

#ifndef INFINITY
#define INFINITY 1e200
#endif


/*****************************************************************
 * Delay Function definition
 *
 ****************************************************************/
DelayFunction::DelayFunction(const MultiCommoNetwork &n, 
                             Real a, Real b, Real g, Real l) :
	net(n), bpr(n,a,b), gamma(g), lambda(l)
{
}

Real DelayFunction::f(Vector &x) const{
	int V = net.getNVertex();
	vector<Real> sumnode(V, 0.0);
	Real sum = 0.0;
  
	ITER(x, itx)
		sumnode[net.arcs[itx.index()].head] += itx.value()/net.arcs[itx.index()].cap;

	FOR(v, V) sum += pow(sumnode[v], lambda);
	return bpr.f(x) + gamma*sum;
}

Vector DelayFunction::g(Vector &x) const {
	int V = net.getNVertex();
	vector<Real> sumnode(V, 0.0);
  
	ITER(x, itx)
		sumnode[net.arcs[itx.index()].head] += itx.value()/net.arcs[itx.index()].cap;

	Vector d = bpr.g(x);
	ITER(d, itd) 
		itd.value()+= gamma*lambda*pow(sumnode[net.arcs[itd.index()].head], lambda-1);

	return d;
}

Vector DelayFunction::gg(Vector &x) const{
	int V = net.getNVertex();
	vector<Real> sumnode(V, 0.0);
  
	ITER(x, itx)
		sumnode[net.arcs[itx.index()].head] += itx.value()/net.arcs[itx.index()].cap;

	Vector d = bpr.gg(x);
	ITER(d, itd) 
		itd.value() += ( gamma * lambda * (lambda-1) *
		                 pow(sumnode[net.arcs[itd.index()].head], lambda-2) );

	return d;
}

/*****************************************************************
 * BPR Function definition
 *
 ****************************************************************/

BPRFunction::BPRFunction(const MultiCommoNetwork &n, Real a, Real b): 
	net(n), alpha(a), beta(b) {}

Real BPRFunction::f(Vector &x) const {
	int A = net.arcs.size();
	assert(A == x.size()); // debug
	Real sum = 0.0, ya, ca, ta;
	ITER(x, itx) {
		int a = itx.index();
		ya = itx.value();
		ca = net.arcs[a].cap;
		ta = net.arcs[a].cost;
		sum += ta*ya*(1 + alpha/(beta+1)*pow(ya/ca,beta));
	}
	return sum;
}

Vector BPRFunction::g(Vector &x) const {
	int A = net.arcs.size(), a = 0, pa = 0;
	assert(A == x.size()); // debug
	Vector d(A);
	Real ya, ca, ta;
	ITER(x, itx){
		pa = a;
		a = itx.index();
		for(int aa = pa+1; aa<a; aa++) d.insert(aa) = net.arcs[aa].cost;
		ya = itx.value(); 
		ca = net.arcs[a].cap; 
		ta = net.arcs[a].cost;
		d.insert(a) = ta + ta*alpha*pow(ya/ca,beta);
	}
	for(int aa = a+1; aa<A; aa++) d.insert(aa) = net.arcs[aa].cost;
	return d;
}

Vector BPRFunction::gg(Vector &x) const {
	int A = net.arcs.size();
	assert(A == x.size()); // debug
	Vector d(A);
	Real ya, ca, ta;
	ITER(x, itx) {
		int a = itx.index();
		ya = itx.value(); ca = net.arcs[a].cap; ta = net.arcs[a].cost;
		assert(ca>0.0);
		d.insert(a) = ta*alpha*beta*pow(ya/ca,beta-1)/ca;
	}
	return d;
}

/*****************************************************************
 * Kleinrock Function definition
 *
 ****************************************************************/

Real KleinrockFunction::f(Vector &x) const {
	int A = net.arcs.size();
	assert(A == x.size()); // debug
	Real sum = 0.0, ya, ca;
	ITER(x, itx){
		int a = itx.index();
		ya = itx.value();
		ca = net.arcs[a].cap;
		//assert(ca-ya>0);
		if(ca>ya) sum += ya/(ca-ya); else return INFINITY;
	}
	return sum;
}

Vector KleinrockFunction::g(Vector &x) const{
	int A = net.arcs.size(), preva = 0;
	assert(A == x.size()); // debug
	Vector d(A);
	Real ya, ca, dd;
	ITER(x, itx){
		int a = itx.index();
		for(int aa = preva; aa < a; aa++) d.insert(aa) = 1.0/net.arcs[aa].cap;
		ya = itx.value(); ca = net.arcs[a].cap;
		//assert(ca-ya>0);
		dd = ca - ya;
		if(dd>0) dd = ca/(dd*dd);
		else dd = INFINITY;
		d.insert(a) = dd;
		preva = a+1;
	}
	for(int aa = preva; aa < A; aa++) d.insert(aa) = 1.0/net.arcs[aa].cap;
	return d;
}

Vector KleinrockFunction::gg(Vector &x) const{
	int A = net.arcs.size();
	assert(A == x.size()); // debug
	Vector d(A);
	Real ya, ca, dd;
	ITER(x, itx){
		int a = itx.index();
		ya = itx.value(); ca = net.arcs[a].cap;
		dd = ca - ya;
		if(dd>0) dd = 2*ca/(dd*dd*dd);
		else dd = INFINITY;
		d.insert(a) = dd;
	}

	return d;
}

KleinrockFunction::KleinrockFunction(const MultiCommoNetwork &n) : net(n){
}

/*****************************************************************
 * ReducableFunction
 *
 ****************************************************************/

Vector compressed_vector(Vector &x, int K){
	assert(x.size() % K == 0);
	Vector y(x.size()/K);
	Vector::iterator itx = x.get_iterator();
	while(!itx.end()){
		int a = itx.index()/K;
		Real ya = 0.0;
		while(!itx.end() && itx.index()/K == a) ya += itx.value(), ++itx;
		if(fabs(ya) > 1e-10) y.insert(a) = ya;
	}
	return y;
}

template <class F>
MCFReducableFunction<F>::MCFReducableFunction(const MultiCommoNetwork &n,
                                              const F &f) :
	net(n), func(f)
{
}

template <class F>
Function* MCFReducableFunction<F>::reduced_function() const{
	return dynamic_cast<Function*>(new F(func));
}

template <class F>
Vector MCFReducableFunction<F>::reduced_variable(Vector &x) const{
	assert(x.size() == net.arcs.size() * net.commoflows.size());
	return compressed_vector(x, net.commoflows.size());
}

template <class F>
Real MCFReducableFunction<F>::f(Vector &x) const{
	Vector y = this->reduced_variable(x);
	return func.f(y);
}

template <class F>
Vector MCFReducableFunction<F>::g(Vector &x) const{
	int A = net.arcs.size(), K = net.commoflows.size();
	Vector y = this->reduced_variable(x);
	Vector dy = func.g(y);
	Vector dx(A*K);
	ITER(dy, itdy) FOR(k, K)
		dx.insert(itdy.index()*K + k) = itdy.value();
	return dx;
}

template <class F>
Vector MCFReducableFunction<F>::gg(Vector &x) const {
	int A = net.arcs.size(), K = net.commoflows.size();
	Vector y = this->reduced_variable(x);
	Vector dy = func.gg(y);
	Vector dx(A*K);
	ITER(dy, itdy) FOR(k, K)
		dx.insert(itdy.index()*K + k) = itdy.value();
	return dx;
}

/*****************************************************************
 * Some ReducableFunctions that are based on DelayFunction,
 * BPRFunction, KleinrockFunction
 *
 ****************************************************************/

MCFDelayFunction::MCFDelayFunction(const MultiCommoNetwork &n, 
                                   Real a, Real b, Real g, Real l) : 
	MCFReducableFunction<DelayFunction>(n, DelayFunction(n,a,b,g,l))
{
}

MCFBPRFunction::MCFBPRFunction(const MultiCommoNetwork &n, 
                               Real a, Real b) :
	MCFReducableFunction<BPRFunction>(n, BPRFunction(n,a,b))
{
}

MCFKleinrockFunction::MCFKleinrockFunction(const MultiCommoNetwork &n) :
	MCFReducableFunction<KleinrockFunction>(n, KleinrockFunction(n))
{
}

/*****************************************************************
 * Section search & golden section search
 *
 ****************************************************************/


#define PHI 0.6180339887498948482045868343656

Real golden_section_search ( Vector &A,
                             Vector &B,
                             Function *obj,
                             int iterations)
{
	Vector *x1 = new Vector(A), *x4 = new Vector(B);
	Vector *x2 = new Vector(B), *x3 = new Vector(B), *xtmp;
  
	(*x2) -= A; (*x2) *= (1-PHI); (*x2) += A;
	(*x3) += A; (*x3) -= (*x2);

	Real f1 = obj->f(A), f4 = obj->f(B);
	Real f2 = obj->f(*x2), f3 = obj->f(*x3), fm = f1;
	Real b1 = 0, b4 = 1, b2 = 1-PHI, b3 = PHI, bm = 0;

	if(fm > f2) fm = f2, bm = b2;
	if(fm > f3) fm = f3, bm = b3;
	if(fm > f4) fm = f4, bm = b4;

	FOR(i, iterations){
		if(f2>f3){
			xtmp = x1; x1 = x2; x2 = x3; x3 = xtmp;
			f1 = f2; f2 = f3;
			b1 = b2; b2 = b3;

			(*x3) = (*x1); (*x3) += (*x4); (*x3) -= (*x2);
			f3 = obj->f(*x3);
			b3 = b1 + b4 - b2;      

			if(fm > f3) fm = f3, bm = b3;
		}
		else{
			xtmp = x4; x4 = x3; x3 = x2; x2 = xtmp;
			f4 = f3; f3 = f2;
			b4 = b3; b3 = b2;

			(*x2) = (*x1); (*x2) += (*x4); (*x2) -= (*x3);
			f2 = obj->f(*x2);
			b2 = b1 + b4 - b3;

			if(fm > f2) fm = f2, bm = b2;
		}
	}

	delete x1; delete x2; delete x3; delete x4;
	return bm;
}

extern fstream iteration_report;

Real general_section_search ( Vector &A,
                              Vector &B,
                              Function *obj,
                              int iterations,
                              Real b2, Real b3)
{
	//iteration_report << "general section search" << endl;
	Vector *x1 = new Vector(A), *x4 = new Vector(B);
	Vector *x2 = new Vector(B), *x3 = new Vector(B), *xtmp;
  
	(*x2) -= A; (*x2) *= b2; (*x2) += A;
	(*x3) -= A; (*x3) *= b3; (*x3) += A;

	Real f1 = obj->f(A), f4 = obj->f(B);
	Real f2 = obj->f(*x2), f3 = obj->f(*x3), fm = f1;
	Real fbound=0, fnewbound;
	Real b1 = 0, b4 = 1, bm = 0, tmp;

	if(fm > f2) fm = f2, bm = b2;
	if(fm > f3) fm = f3, bm = b3;
	if(fm > f4) fm = f4, bm = b4;

	FOR(i, iterations){
		// make sure that the order is x1 -> x2 -> x3 -> x4
		if(b2>b3){
			xtmp = x2; x2 = x3; x3 = xtmp;
			tmp  = b2; b2 = b3; b3 = tmp;
			tmp  = f2; f2 = f3; f3 = tmp;
		}

		// safe guard the case where x2 and x3 are too close to each other
		if(fabs(b2-b3)<1e-12 && fabs(f2-f3) < 1e-6){
			iteration_report << "x2 and x3 too close at iteration "<<i<<endl;
			(*x3) -= (*x4); (*x3) *= (PHI-1); (*x3) += (*x2);
			f3 = obj->f(*x3);
			b3 -= b4; b3 *= (PHI-1); b3 += b2;
			if(fm > f3) fm = f3, bm = b3;
		}

		if(f2>f3){
			fnewbound = (b4-b1)/(b2-b1)*(f2-f1)+f1;
			if(fnewbound > fbound) fbound = fnewbound;

			xtmp = x1; x1 = x2; x2 = x3; x3 = xtmp;
			f1 = f2; f2 = f3;
			b1 = b2; b2 = b3;

			(*x3) = (*x1); (*x3) += (*x4); (*x3) -= (*x2);
			f3 = obj->f(*x3);
			b3 = b1 + b4 - b2;      

			if(fm > f3) fm = f3, bm = b3;
		}
		else{
			fnewbound = (b1-b4)/(b3-b4)*(f3-f4)+f4;
			if(fnewbound > fbound) fbound = fnewbound;
	    
			xtmp = x4; x4 = x3; x3 = x2; x2 = xtmp;
			f4 = f3; f3 = f2;
			b4 = b3; b3 = b2;

			(*x2) = (*x1); (*x2) += (*x4); (*x2) -= (*x3);
			f2 = obj->f(*x2);
			b2 = b1 + b4 - b3;

			if(fm > f2) fm = f2, bm = b2;
		}
		if(fabs(fbound-fm)/(fbound+fm) < 0.5*1e-9) break;
	}  

	delete x1; delete x2; delete x3; delete x4;
	return bm;
}


// golden search between A and B
Real section_search ( Vector &A, 
                      Vector &B, 
                      Function *obj, 
                      int iterations,
                      bool to_use_golden_ratio,
                      Real b1, Real b2)
{
	// Check whether the function can be reduced (a reducable function)
	// by casting it to ReducableFunction class
	ReducableFunction *casted_obj = dynamic_cast<ReducableFunction*> (obj);

	if(casted_obj != NULL){
		// If it can be casted --> it is a reducable function
		// then do the search with the reduced function and variables instead
		Function* reduced_obj = casted_obj->reduced_function();
		Vector A_ = casted_obj->reduced_variable(A);
		Vector B_ = casted_obj->reduced_variable(B);
		Real lambda = section_search( A_, B_,
		                              reduced_obj,
		                              iterations,
		                              to_use_golden_ratio,
		                              b1, b2);
		delete reduced_obj;
		return lambda;
	}

	if(to_use_golden_ratio) return golden_section_search(A,B,obj,iterations);
	return general_section_search (A,B,obj,iterations,b1,b2);
}

// Naive line search between A and B
Real line_search (Vector &A, Vector &B, Function *obj, int niteration){
	Vector x(A), dx(B);
	Real fmin = obj->f(x), imin = 0.0, f;
	dx -= A; dx *= (1.0/niteration);
	FOR(i, niteration){
		x += dx;
		f = obj->f(x);
		if(f < fmin) fmin = f, imin = i+1;
	}
	return imin / niteration;
}
