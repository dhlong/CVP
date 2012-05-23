#ifndef __FUNCTION_H__
#define __FUNCTION_H__

#include "network.h"

#define PHI 0.6180339887498948482045868343656

using namespace std;

// Prototype of a differientiable function
class Function{
public:
  // return value of the function at point x
  virtual Real f(const Vector &x) const = 0;

  // return gradient of the function at point x
  virtual Vector g(const Vector &x) const = 0;

  // Diagonal of Hessian matrix
  virtual Vector gg(const Vector &x) const = 0;

  // destructor
  virtual ~Function() {} ;
};

// quartic function with delay propagation
class QuarticFunction : public Function{
private:
  Network net;

  // to[i] containts the list of arcs that link to arc i
  vector< vector<int> > to;

public:
  QuarticFunction(const Network &n);
  virtual Real f(const Vector &v) const;
  virtual Vector g(const Vector &v) const;
  virtual Vector gg(const Vector &v) const;
};

// BPR Function on a multi-commodity network
class BPRFunction: public Function{
private:
  MultiCommoNetwork net;
  Real alpha, beta;

public:
  virtual Real f(const Vector &x) const;
  virtual Vector g(const Vector &x) const;
  virtual Vector gg(const Vector &x) const;
  BPRFunction(const MultiCommoNetwork &n, Real a=0.15, Real b=4);
  MultiCommoNetwork getNetwork() const{
    return net;
  }
};

class ReducedBPRFunction: public Function{
 private:
  MultiCommoNetwork net;
  Real alpha, beta;

 public:
  virtual Real f(const Vector &x) const;
  virtual Vector g(const Vector &x) const;
  virtual Vector gg(const Vector &x) const;
  ReducedBPRFunction(const MultiCommoNetwork &n, const Real a=0.15, Real b=4);
};

class KleinrockFunction : public Function {
 private:
  MultiCommoNetwork net;
  Real alpha, beta;

 public:
  virtual Real f(const Vector &x) const;
  virtual Vector g(const Vector &x) const;
  virtual Vector gg(const Vector &x) const;
  KleinrockFunction(const MultiCommoNetwork &n, const Real a=0.15, Real b=4);  
};

double golden_search ( const Vector &x0, 
		       const Vector &x1, 
		       Function *obj, 
		       int niteration = 20,
		       double ratio = PHI);
		       
double line_search (const Vector &x0, const Vector &x1, Function *obj, int niteration = 20);

#endif

