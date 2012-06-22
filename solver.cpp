#include "solver.h"
#include <cassert>
#include <cstdio>
#include <cstdlib>

#define FOR(i,n) for(int i = 0, n_ = (n); i < n_; i++)
#define FREE(a) do{ if(a!=NULL) free(a); a = NULL; } while(0)

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

GRBSolver::GRBSolver(){
	int error = GRBloadenv(&env, "grb_solver.log");
	GRBsetintparam(env, "Threads", 1);
	GRBsetintparam(env, "OutputFlag", 0);
	assert(error == 0 && env != NULL);

	error = GRBnewmodel(env, &model, "grb model", 0, 
	                    NULL, NULL, NULL, NULL, NULL);
	assert(error == 0 && model != NULL);

}

GRBSolver::~GRBSolver(){
	GRBfreemodel(model);
	GRBfreeenv(env);
	model = NULL; env = NULL;
}

int GRBSolver::copylp(int cols, int rows, int objsense, double* obj, double *rhs,
                      char *sense, int *matbeg, int *matcnt, int *matind, double *matval,
                      double *lb, double *ub){
	int *rbeg = (int*) malloc(rows*sizeof(int));
	FOR(i, rows) rbeg[i] = 0;

	int error = GRBaddconstrs(model, rows, 0, rbeg, NULL, NULL, sense, rhs, NULL);
	assert(error == 0);
	
	int nz = matbeg[cols-1] + matcnt[cols-1];
	cout<<"N Constraints = "<<rows<<"; N vars = "<<cols<<"; NZeros = "<<nz<<endl;
	GRBgetintattr(model, "NumConstrs", &rows);

	GRBupdatemodel(model);
	              

	error = GRBaddvars(model, cols, matbeg[cols-1] + matcnt[cols-1],
	                   matbeg, matind, matval,
	                   obj, lb, ub, NULL, NULL);
	GRBupdatemodel(model);

	if(error != 0)
		cout<<"Error code: "<<error<<"; Error: "<<GRBgeterrormsg(env)<<endl;
	assert(error == 0);
	FREE(rbeg);
	return error;
}

int GRBSolver::copyquad(int *qmatbeg, int *qmatcnt, int *qmatind, double *qmatval){
	int cols = 0, idx = 0;
	GRBgetintattr(model, "NumVars", &cols);

	int qnz = qmatbeg[cols-1] + qmatcnt[cols-1];
	int *colind = (int*) malloc(qnz*sizeof(int));
	FOR(col, cols) FOR(row, qmatcnt[col]) colind[idx++] = col;
	FOR(i, qnz) qmatval[i] *= 0.5;

	int error = GRBaddqpterms(model, qnz, colind, qmatind, qmatval);
	assert(error == 0);
	return error;
}

int GRBSolver::chgobj(int cnt, int *ind, double *obj){
	int error = GRBsetdblattrlist(model, "Obj", cnt, ind, obj);
	assert(error == 0);
	return error;
}

int GRBSolver::chgrhs(int cnt, int *ind, double *rhs){
	int error = GRBsetdblattrlist(model, "RHS", cnt, ind, rhs);
	assert(error == 0);
	return error;
}

int GRBSolver::solve(){
	int error = GRBoptimize(model);
	assert(error == 0);
	return error;
}

int GRBSolver::getx(double *x, int beg, int end){
	int cols = GRBgetintattr(model, "NumVars", &cols);
	int *ind = (int*) malloc((end-beg+1)*sizeof(int));
	FOR(i, end-beg+1) ind[i] = beg+i;
	int error = GRBgetdblattrlist(model, "X", end-beg+1, ind, x);
	assert(error == 0);
	FREE(ind);
	return error;
}

