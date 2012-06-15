#include <cassert>
#include <cmath>
#include <algorithm>
#include "my_sparse_vector.h"

inline double sqr(double a){
	return a*a;
}

MySparseVector::iterator MySparseVector::get_iterator(){
	sort();
	return iterator(this);
}

void MySparseVector::sort(){
	if(!sorted) val.sort();
  sorted = true;
}

MySparseVector::MySparseVector(int n) : size_(n), val(), sorted(true){
}

MySparseVector& MySparseVector::operator += (MySparseVector &x){
  assert(x.size_ == size_);
  sort(); x.sort();
  list< pair<int, double> >::iterator it = val.begin(), itx = x.val.begin();
  
  while(it != val.end() && itx != x.val.end()){
    if(it->first == itx->first){
      it->second += itx->second;
      ++itx;
      if(fabs(it->second) < 1e-10) it = val.erase(it);
      else ++it;
    }
    else if(it->first > itx->first){
      val.insert(it, *itx);
      ++itx;
    }
    else ++it;
  }

  while(itx != x.val.end()){
    val.push_back(*itx);
    ++itx;
  }
  return *this;
}

MySparseVector& MySparseVector::operator -= (MySparseVector &x){
  assert(x.size_ == size_);
  sort(); x.sort();
  list< pair<int, double> >::iterator it = val.begin(), itx = x.val.begin();

  
  while(it != val.end() && itx != x.val.end()){
    if(it->first == itx->first){
      it->second -= itx->second;
      ++itx;
      if(fabs(it->second) < 1e-10) it = val.erase(it);
      else ++it;
    }
    else if(it->first > itx->first){
      val.insert(it, make_pair(itx->first, -itx->second));
      ++itx;
    }
    else ++it;
  }

  while(itx != x.val.end()){
    val.push_back(make_pair(itx->first, -itx->second));
    ++itx;
  }
  return *this;
}

MySparseVector& MySparseVector::operator *= (double alpha){
  if(alpha == 0.0){
    val.clear();
    return *this;
  }

  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
    it->second *= alpha;
  return *this;
}

double& MySparseVector::insert(int i){
	assert(i<size_);
	if(val.empty()){
		val.push_back(make_pair(i, 0.0));
		return val.back().second;
	}
  if(sorted){
    if(val.front().first > i){
      val.push_front(make_pair(i, 0.0));
      return val.front().second;
    }
    if(val.back().first < i){
      val.push_back(make_pair(i, 0.0));
      return val.back().second;
    }
  }
  val.push_back(make_pair(i, 0.0));
  sorted = false;
  return val.back().second;
}

double& MySparseVector::operator[](int i){
	assert(i<size_);
  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
    if(it->first == i) return it->second;
  return insert(i);
}

double MySparseVector::dot(MySparseVector & x){
  assert(x.size_ == size_);
  sort(); x.sort();
  list< pair<int, double> >::iterator it = val.begin(), itx = x.val.begin();
  double sum = 0.0;
  
  while(it != val.end() && itx != x.val.end()){
    if(it->first == itx->first){
      sum += ((it->second) * (itx->second));
      ++it; ++itx;
    }
    else if(it->first > itx->first) ++itx;
    else ++it;
  }
  return sum;
}

bool MySparseVector::operator ==(MySparseVector & x){
  sort(); x.sort();
  list< pair<int, double> >::iterator it = val.begin(), itx = x.val.begin();
  double sum = 0.0;
  
  while(it != val.end() && itx != x.val.end()){
    if(it->first == itx->first){
	    if(fabs(it->second - itx->second) > 1e-8) return false;
      ++it; ++itx;
    }
    else if(it->first > itx->first) 
	    if(fabs(itx->second) > 1e-8) return false;
	    else ++itx;
    else 
	    if(fabs(it->second) > 1e-8) return false;
	    else ++it;
  }

  while(it != val.end()) 
	  if(fabs(it->second) > 1e-8) return false;
	  else ++it;

  while(itx != x.val.end()) 
	  if(fabs(it->second) > 1e-8) return false;
	  else ++itx;
		                              
  return true;
}

double MySparseVector::squaredNorm() {
  double sum = 0.0;
  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
	  sum += sqr(it->second);
  return sum;
}

double MySparseVector::norm(){
  return sqrt(squaredNorm());
}

void MySparseVector::output(){
  sort();
  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
    cout<<"["<<it->first<<"]="<<it->second<<endl;
}

int MySparseVector::size() const{
	return size_;
}

int MySparseVector::nonZeros() const{
	return val.size();
}

bool MySparseVector::check(){
  sort();
  int i = -1;
  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
    if(it->first <= i) return false;
    else i = it->first;
  return true;
  
}

double MySparseVector::coeff(int i){
  for(list< pair<int, double> >::iterator it = val.begin();
      it != val.end(); ++it)
	  if(it->first == i) return it->second;
  return 0.0;
}

double& MySparseVector::coeffRef(int i){
	return (*this)[i];
}

