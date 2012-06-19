#include "solver.h"
#include <cassert>

CPXSolver::CPXSolver(){
	int status;
	env = CPXopenCPLEX (&status);
	assert(env != NULL);
  
	lp  = CPXcreateprob (env, &status, "new problem");
	assert(lp != NULL);

	assert(!CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF));
	assert(!CPXsetintparam (env, CPX_PARAM_THREADS, 1));
}

CPXSolver::~CPXSolver(){
	CPXcloseCPLEX(&env);
	env = NULL;
	lp = NULL;
}

int CPXSolver::copylp(int cols, int rows, int objsense, double* obj, double *rhs,
                      char *sense, int *matbeg, int *matcnt, int *matind, double *matval,
                      double *lb, double *ub){
	int error = CPXcopylp(env, lp, cols, rows, objsense, obj, rhs,
	                      sense, matbeg, matcnt, matind, matval,
	                      lb, ub, NULL);
	assert(error == 0);
	return error;
}

int CPXSolver::copyquad(int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval){
	int error = CPXcopyquad(env, lp, qmatbeg, qmatcnt, qmatind, qmatval);
	assert(error == 0);
	return error;
}

int CPXSolver::chgobj(int cnt, int *ind, double *obj){
	int error = CPXchgobj(env, lp, cnt, ind, obj);
	assert(error == 0);
	return error;
}

int CPXSolver::chgrhs(int cnt, int *ind, double *rhs){
	int error = CPXchgrhs(env, lp, cnt, ind, rhs);
	assert(error == 0);
	return error;
}

int CPXSolver::solve(){
	int error  = CPXbaropt(env, lp);
	assert(error == 0);
	return error;
}

int CPXSolver::getx(double *x, int beg, int end){
	int error = CPXgetx(env, lp, x, beg, end);
	assert(error == 0);
	return error;
}



