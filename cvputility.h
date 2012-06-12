#ifndef __CVPUTILITY_H__
#define __CVPUTILITY_H__

#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <cstdarg>
#include <map>
#include <typeinfo>
#include <cctype>
#include <algorithm>
#include <utility>

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <Eigen/Sparse>

#define FOR(i,n) for(int i=0, n_=(n); i<n_; i++)
#define MALLOC(type, size) (type*) malloc((size)*sizeof(type))
#define ITER(vec, it) for(Vector::InnerIterator it(vec); it; ++it)
#define FREE(var)				\
  do {						\
    if(var) free(var);				\
    (var) = NULL;				\
  } while(0)

using namespace std;
using namespace Eigen;

typedef double Real;
typedef SparseVector<Real> Vector;
//typedef vector<Real> Vector;
typedef int Vertex;

// report error message and terminate the program 
void error_handle(const string &);

class TableReport {
 private:
  string row_format, header_format, line;

 public:
  TableReport(const string &format);
  ostream& print_header(ostream &f, ...);
  ostream& print_row(ostream &f, ...);
  ostream& print_line(ostream &f);
};

class Timer{
 private:
  vector<double> records;

 protected:
  virtual double get_clock() const;

 public:
  void record();
  double elapsed(int i=-1, int j=-2);
};

class CPUTimer : public Timer{
 protected:
  virtual double get_clock() const;
};

double memory_usage(int peak = 1);

class SettingMapper{
 private:
  map<string, int> int_params;
  map<string, double> real_params;
  map<string, string> str_params;
  map<string, bool> bool_params;

 public:
  
  void setb(const string &key, bool value);
  void seti(const string &key, int value);
  void setr(const string &key, double value);
  void sets(const string &key, const string &value);
  void set_by_str(const string &key, const string &value_str);

  bool& getb(const string &key);
  int& geti(const string &key);
  double& getr(const string &key);
  string& gets(const string &key);

  void read(fstream &f);
  void report(fstream &f);
};


/**
 * Some utility functions to work with vector
 */

/*
inline Vector operator + (const Vector &x, const Vector &y){
  Vector z(x);
  FOR(i, x.size()) z[i] += y[i];
  return z;
}

inline Vector operator - (const Vector &x, const Vector &y){
  Vector z(x);
  FOR(i, x.size()) z[i] -= y[i];
  return z;
}

inline Vector & operator += (Vector &x, const Vector &y){
  FOR(i, x.size()) x[i] += y[i];
  return x;
}

inline Vector & operator -= (Vector &x, const Vector &y){
  FOR(i, x.size()) x[i] -= y[i];
  return x;
}

inline Vector & operator *= (Vector &x, Real alpha){
  FOR(i, x.size()) x[i] *= alpha;
  return x;
}

inline Real operator * (const Vector &x, const Vector &y){
  Real sum = 0.0;
  FOR(i, x.size()) sum += x[i]*y[i];
  return sum;
}

inline Vector operator *(Real alpha, const Vector &x){
  Vector z(x);
  z*=alpha;
  return z;
}

*/

bool operator == (const Vector &x, const Vector &y);

template <typename T> inline bool updatemin(T &a, T b){
  if(a>b) return a=b, true;
  return false;
}

template <typename T> inline bool updatemax(T &a, T b){
  if(a<b) return a=b, true;
  return false;
}

template<typename T>
bool read_param(string fininame, const string &alias, T &t){
  fstream f(fininame.c_str(), fstream::in);
  char line[305], *str;

  if(!f) error_handle("No CVP.ini file found.");
  while(!f.eof()){
    f.getline(line, 300);
    str = strrchr(line,'#'); // comments are marked by '#' character
    if(str != NULL) continue;
    if(strstr(line, alias.c_str()) == NULL) continue;
    str = strrchr(line,'='); // value delim by the '=' character
    if(str == NULL) continue; // no delim found, read next line

    stringstream ss(str+1, stringstream::in);
    ss>>t;
    f.close();
    return true;
  }

  f.close();
  return false;
}

#endif

