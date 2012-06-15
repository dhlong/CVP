#ifndef __MY_SPARSE_VECTOR__
#define __MY_SPARSE_VECTOR__

#include <list>
#include <iostream>
#include <utility>
#include <cstdio>
#include <iterator>

using namespace std;

class MySparseVector{
 private:
  list< pair<int, double> > val;
  bool sorted;
	int size_;

  void sort();

 public:
	typedef pair<int, double> IVPair;
	typedef list<IVPair> IVPairL;
	class iterator : public IVPairL::iterator {
		MySparseVector *v;
	public:
		iterator(MySparseVector *v_): v(v_), IVPairL::iterator(v_->val.begin()){}
		inline bool end(){ return (*this) == v->val.end();}
		inline int index(){ return (*this)->first; }
		inline double value(){ return (*this)->second; }
	};
	
	iterator get_iterator();

  MySparseVector(int n);

  MySparseVector & operator += (MySparseVector &);
  MySparseVector & operator -= (MySparseVector &);
  MySparseVector & operator *= (double);
	bool operator == (MySparseVector &);

  double & insert(int);
  double & operator[](int);
	double coeff(int);
	double & coeffRef(int);

  double dot(MySparseVector &);
  double squaredNorm();
  double norm();

	int size() const;
	int nonZeros() const;
	bool check();

  void output();
  
};

#endif
