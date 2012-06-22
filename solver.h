#ifndef __SOLVER__H__
#define __SOLVER__H

#include <iostream>
#include <cstdio>
#include <ilcplex/cplex.h>

extern "C"{
#include <gurobi_c.h>
}


using namespace std;

class Solver {
public:
	virtual int copylp(int cols, int rows, int objsense, double* obj, double *rhs,
	                   char *sense, int *matbeg, int *metcnt, int *matind, double *matval,
	                   double *lb, double *ub) = 0;
	virtual int copyquad(int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval) = 0;
	virtual int chgobj(int cnt, int *ind, double *obj) = 0;
	virtual int chgrhs(int cnt, int *ind, double *rhs) = 0;
	virtual int solve() = 0;
	virtual int getx(double *x, int beg, int end) = 0;
	virtual ~Solver() {}
};

class CPXSolver : public Solver{
private:
	CPXENVptr env;
	CPXLPptr  lp;

public:
	CPXSolver();

	virtual int copylp(int cols, int rows, int objsense, double* obj, double *rhs,
	                   char *sense, int *matbeg, int *matcnt, int *matind, double *matval,
	                   double *lb, double *ub);
	virtual int copyquad(int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval);
	virtual int chgobj(int cnt, int *ind, double *obj);
	virtual int chgrhs(int cnt, int *ind, double *rhs);
	virtual int solve();
	virtual int getx(double *x, int beg, int end);
	virtual ~CPXSolver();
};

class GRBSolver : public Solver{
private:
	GRBenv *env;
	GRBmodel *model;

public:
	GRBSolver();

	virtual int copylp(int cols, int rows, int objsense, double* obj, double *rhs,
	                   char *sense, int *matbeg, int *metcnt, int *matind, double *matval,
	                   double *lb, double *ub);
	virtual int copyquad(int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval);
	virtual int chgobj(int cnt, int *ind, double *obj);
	virtual int chgrhs(int cnt, int *ind, double *rhs);
	virtual int solve();
	virtual int getx(double *x, int beg, int end);
	virtual ~GRBSolver();
};

#endif
