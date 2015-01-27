// gcore_Exclusion.cc

#include "Exclusion.h"
#include <vector>
#include <algorithm>

using gcore::Exclusion;
using namespace std;

class gcore::Exclusion_i{
  friend class gcore::Exclusion;
  vector<int> d_excl;
  Exclusion_i(): d_excl(){}
  Exclusion_i(const Exclusion_i &e):d_excl(e.d_excl){}
};

gcore::Exclusion::Exclusion(): d_this(new Exclusion_i()){}

Exclusion::Exclusion(const Exclusion &e): d_this (new Exclusion_i(*e.d_this)){}

Exclusion::~Exclusion(){delete d_this;}

Exclusion &Exclusion::operator=(const Exclusion &e){
  if (this != &e){
    this->~Exclusion();
    new(this) Exclusion(e);
  }
  return *this;
}

void Exclusion::insert(int i){
  vector<int>::iterator iter=d_this->d_excl.begin(), 
    to=d_this->d_excl.end();
  while(iter!=to){
    if(*iter == i) return;
    if(*iter > i ){
      d_this->d_excl.insert(iter,i);
      return;
    }
    ++iter;
  }
  d_this->d_excl.push_back(i);
}

void Exclusion::erase(int i){
  vector<int>::iterator iter=d_this->d_excl.begin(), 
    to=d_this->d_excl.end();
  while(iter!=to){
    if(*iter == i){
      d_this->d_excl.erase(iter);
      return;
    }
    ++iter;
  }
}


int Exclusion::size()const{
  return d_this->d_excl.size();
}

int Exclusion::atom(int i)const{
  return d_this->d_excl[i];
}

bool Exclusion::contains(int i)const{
  unsigned int x = 0;
  for(; x < d_this->d_excl.size(); ++x)
    if (d_this->d_excl[x] == i)
      break;
  
  return x != d_this->d_excl.size();
}
