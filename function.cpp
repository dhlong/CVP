#include "function.h"

QuarticFunction::QuarticFunction(const Network &n) : 
  net(n), to(n.getNVertex()) {
  FOR(i,n.arcs.size()) to[net.arcs[i].head].push_back(i);
}

Real QuarticFunction::f(const Vector &v) const{
  Real sum = 0.0, x, c;
  int arc;
  
  assert(v.size()==net.arcs.size());

  FOR(i, v.size()){
    FOR(j, to[net.arcs[i].tail].size()){
      arc = to[net.arcs[i].tail][j];
      x = (v[arc]/net.arcs[arc].cap);
      sum += 20*x*x*x*x;
    }
    x = (v[i]/net.arcs[i].cap);
    sum += 100*x*x;
    
    c = Real(net.arcs[i].head+1.0)/Real(net.arcs[i].tail+1.0);
    sum += c*v[i];
  }
  return sum;
}

Vector QuarticFunction::g(const Vector &v) const{
  Vector d(v.size(),0.0);
  Real x, c, sum;
  int a;

  FOR(i, v.size()) {
    sum = 0;
    FOR(j, to[net.arcs[i].tail].size()){
      a = to[net.arcs[i].tail][j];
      x = v[a]; c = net.arcs[a].cap;
      d[a] += 80*x*x*x/(c*c*c*c);
    }

    x = v[i]; c = net.arcs[i].cap;
    d[i] += 200*x/(c*c);

    c = Real(net.arcs[i].head+1.0)/Real(net.arcs[i].tail+1.0);
    d[i] += c;
  }
  return d;
}

//!!!! have not debugged
// do not use yet
Vector QuarticFunction::gg(const Vector &v) const{
  Vector d(v.size(),0.0);
  Real x, c, sum;
  int a;

  FOR(i, v.size()) {
    sum = 0;
    FOR(j, to[net.arcs[i].tail].size()){
      a = to[net.arcs[i].tail][j];
      x = v[a]; c = net.arcs[a].cap;
      d[a] += 240*x*x/(c*c*c*c);
    }

    x = v[i]; c = net.arcs[i].cap;
    d[i] += 200/(c*c);

    c = Real(net.arcs[i].head+1.0)/Real(net.arcs[i].tail+1.0);
    d[i] += c;
  }
  return d;
}

BPRFunction::BPRFunction(const MultiCommoNetwork &n, Real a, Real b): 
  net(n), alpha(a), beta(b) {}

Real BPRFunction::f(const Vector &x) const {
  int K = net.commoflows.size(), A = net.arcs.size();
  assert(K*A == x.size()); // debug
  Real sum = 0.0, ya, ca, ta;
  FOR(a, A){
    ya = 0.0;
    FOR(i, K) ya += x[a*K + i];
    ca = net.arcs[a].cap;
    ta = net.arcs[a].cost;
    sum += ta*ya*(1 + alpha/(beta+1)*pow(ya/ca,beta));
  }
  return sum;
}

Vector BPRFunction::g(const Vector &x) const {
  int K = net.commoflows.size(), A = net.arcs.size();
  assert(K*A == x.size()); // debug
  Vector d(K*A);
  Real ya, ca, dd, ta;
  FOR(a, A){
    ya = 0.0; ca = net.arcs[a].cap; ta = net.arcs[a].cost;
    FOR(i, K) ya += x[a*K + i];
    dd = ta + ta*alpha*pow(ya/ca,beta);
    FOR(i,K) d[a*K + i] = dd;
  }
  return d;
}

Vector BPRFunction::gg(const Vector &x) const {
  int K = net.commoflows.size(), A = net.arcs.size();
  assert(K*A == x.size()); // debug
  Vector d(K*A);
  Real ya, ca, dd, ta;
  FOR(a, A){
    ya = 0.0; ca = net.arcs[a].cap; ta = net.arcs[a].cost;
    FOR(i, K) ya += x[a*K + i];
    assert(ca>0.0);
    dd = ta*alpha*beta*pow(ya/ca,beta-1)/ca;
    FOR(i,K) d[a*K + i] = dd;
  }
  return d;
}

ReducedBPRFunction::ReducedBPRFunction(const MultiCommoNetwork &n, Real a, Real b): 
  net(n), alpha(a), beta(b) {}

Real ReducedBPRFunction::f(const Vector &x) const {
  int A = net.arcs.size();
  assert(A == x.size()); // debug
  Real sum = 0.0, ya, ca, ta;
  FOR(a, A){
    ya = x[a];
    ca = net.arcs[a].cap;
    ta = net.arcs[a].cost;
    sum += ta*ya*(1 + alpha/(beta+1)*pow(ya/ca,beta));
  }
  return sum;
}

Vector ReducedBPRFunction::g(const Vector &x) const {
  int A = net.arcs.size();
  assert(A == x.size()); // debug
  Vector d(A);
  Real ya, ca, ta;
  FOR(a, A){
    ya = x[a]; ca = net.arcs[a].cap; ta = net.arcs[a].cost;
    d[a] = ta + ta*alpha*pow(ya/ca,beta);
  }
  return d;
}

Vector ReducedBPRFunction::gg(const Vector &x) const {
  int A = net.arcs.size();
  assert(A == x.size()); // debug
  Vector d(A);
  Real ya, ca, ta;
  FOR(a, A){
    ya = x[a]; ca = net.arcs[a].cap; ta = net.arcs[a].cost;
    assert(ca>0.0);
    d[a] = ta*alpha*beta*pow(ya/ca,beta-1)/ca;
  }
  return d;
}

#define PHI 0.6180339887498948482045868343656

//
// recuresively do golden search
// x1 and x2 are the two middle points in the search segment AB
// e.g    A----x1----x2----B
// f1, f2, b1, b2 are the values of the function and the betas
// correspoding to x1 and x2
// fm and bm are to save the minimum value and beta while searching
//
void golden_search_recursive (Vector &x1, Vector &x2, 
			      Real f1, Real f2, Real &fm,
			      Real b1, Real b2, Real &bm,
			      Function *obj,
			      int count,
			      double ratio) 
{
  if(ratio > 0.5) ratio = 1 - ratio;
  if(fabs(3*ratio - 1) < 1e-6) {
    Vector dx(x2); dx -= x1;
    Real dphi = (1-2*ratio)*((1-PHI)-ratio);
    Real db = dphi*(b2-b1);
    dx *= dphi;
    x1 += dx; x2 -= dx;
    b1 += db; b2 -= db;
    f1 = obj->f(x1); f2 = obj->f(x2);
    ratio += dphi/(1-2*ratio);
  }
  if (f1 > f2) {
    x1 -= x2; x1 *= ( (1 - 3*ratio)/(1 - 2*ratio) ); x1 += x2;
    f1 = obj->f(x1);
    b1 = b1 + ratio/(1-2*ratio)*(b2-b1);
    if(f1<fm) fm=f1, bm=b1;
  } 
  else {
    x2 -= x1; x2 *= ( (1 - 3*ratio)/(1 - 2*ratio) ); x2 += x1;
    f2 = obj->f(x2);
    b2 = b2 + ratio/(1-2*ratio)*(b1-b2);
    if(f2<fm) fm=f2, bm=b2;
  }
  
  if(count == 0) return;
  //  if(fabs(b2-b1) < 1e-3*0.5*(fabs(b2)+fabs(b1))) return;
  golden_search_recursive (x2,x1,f2,f1,fm,b2,b1,bm,obj,count-1,ratio/(1-ratio));
}

// golden search between A and B
double golden_search ( const Vector &A, 
		       const Vector &B, 
		       Function *obj, 
		       int niteration,
		       double ratio)
{
  
  BPRFunction *bpr = dynamic_cast<BPRFunction*> (obj);
  if(bpr != NULL && dynamic_cast<ReducedBPRFunction*> (obj) == NULL){
    MultiCommoNetwork net = bpr->getNetwork();
    ReducedBPRFunction rbpr(net);
    int nA = net.arcs.size(), nK = net.commoflows.size();
    Vector AA(nA), BB(nA);
    FOR(a, nA){
      AA[a] = BB[a] = 0.0;
      FOR(k, nK) AA[a] += A[a*nK + k], BB[a] += B[a*nK + k];
    }
    //cout<<setprecision(20)<<obj->f(A)<<" "<<setprecision(20)<<rbpr.f(AA)<<endl;
    assert(fabs(obj->f(A)-rbpr.f(AA))<1e-2);
    assert(fabs(obj->f(B)-rbpr.f(BB))<1e-2);
    return golden_search(AA, BB, &rbpr, niteration, ratio);
  }
  

  Vector x1(A), x2(A);
  if(ratio < 0.5) ratio = 1-ratio;

  x1 -= B; x1*=ratio; x1+=B; // x1 = PHI*A + (1-PHI)*B = PHI*(A-B) - B;
  x2 += B; x2 -= x1; // x2 = (1-PHI)*A + PHI*B = A + B - x1

  Real f1 = obj->f(x1), f2 = obj->f(x2), b1 = 1-ratio, b2 = ratio;
  Real fm = f1, bm = b1;

  golden_search_recursive (x1,x2,f1,f2,fm,b1,b2,bm,obj,niteration,1-ratio);
  return bm;
}

// Naive line search between A and B
double line_search (const Vector &A, const Vector &B, Function *obj, int niteration){
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
