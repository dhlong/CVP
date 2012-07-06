#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include "network.h"

#define PHI 0.6180339887498948482045868343656

using namespace std;

// Prototype of a differientiable function
class Function{
 public:
	// return value of the function at point x
	virtual Real f(Vector &x) const = 0;

	// return gradient of the function at point x
	virtual Vector g(Vector &x) const = 0;

	// Diagonal of Hessian matrix
	virtual Vector gg(Vector &x) const = 0;

	// destructor
	virtual ~Function() {} ;
};

class BPRFunction: public Function {
 private:
	MultiCommoNetwork net;
	Real alpha, beta;
	
 public:
	virtual Real f(Vector &x) const;
	virtual Vector g(Vector &x) const;
	virtual Vector gg(Vector &x) const;
	BPRFunction(const MultiCommoNetwork &n, const Real a=0.15, Real b=4);
};

class DelayFunction : public Function {
private:
	MultiCommoNetwork net;
	BPRFunction bpr;
	Real gamma, lambda;

 public:
	DelayFunction(const MultiCommoNetwork &n, Real a = 0.15, Real b = 4, Real g = 0.01, Real l = 10);
	virtual Real f(Vector &v) const;
	virtual Vector g(Vector &v) const;
	virtual Vector gg(Vector &v) const;
};

class KleinrockFunction : public Function {
 private:
	MultiCommoNetwork net;

 public:
	virtual Real f(Vector &x) const;
	virtual Vector g(Vector &x) const;
	virtual Vector gg(Vector &x) const;
	KleinrockFunction(const MultiCommoNetwork &n);  
};

class ReducableFunction : public Function {
 public:
	virtual Function* reduced_function() const = 0;
	virtual Vector reduced_variable(Vector &) const = 0;
};

template <class F>
class MCFReducableFunction : public ReducableFunction {
private:
	MultiCommoNetwork net;
	F func;
	
public:
	MCFReducableFunction(const MultiCommoNetwork &n, const F &f);
	virtual Function* reduced_function() const;
	virtual Vector reduced_variable(Vector &) const;

	virtual Real f(Vector &x) const;
	virtual Vector g(Vector &x) const;
	virtual Vector gg(Vector &x) const;
};

// quartic function with delay propagation
class MCFDelayFunction : public MCFReducableFunction<DelayFunction> {
public:
	MCFDelayFunction(const MultiCommoNetwork &n, Real a, Real b, Real g, Real l);
};

// BPR Function on a multi-commodity network
class MCFBPRFunction: public MCFReducableFunction<BPRFunction> {
public:
	MCFBPRFunction(const MultiCommoNetwork &n, Real a=0.15, Real b=4);
};

class MCFKleinrockFunction: public MCFReducableFunction<KleinrockFunction>{
public:
	MCFKleinrockFunction(const MultiCommoNetwork &n);  
};


Real section_search ( Vector &x0, 
                      Vector &x1, 
                      Function *obj, 
                      int iterations = 20,
                      bool to_use_golden_ratio = true,
                      Real b1 = 1-PHI, 
                      Real b2 = PHI);

Real line_search (Vector &x0, Vector &x1, Function *obj, int niteration = 20);

#endif

