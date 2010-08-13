// gcore_LJException.cc

#include <new>
#include <iostream>
#include <set>
#include "LJException.h"

using gcore::LJException;

LJException::LJException(int a, int b){
  if(a<b){
    d_a[0]=a;
    d_a[1]=b;
  }
  else{
    d_a[0]=b;
    d_a[1]=a;
    std::cerr << "NOTE: order of atoms changed in LJException:\n";
    std::cerr << "      " << a+1 << "," << b+1 << " -> " << b+1 << "," << a+1 << std::endl;
  }
  d_type=-1;
}

LJException::LJException(const LJException &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_type=a.d_type;
  d_cond=a.d_cond;
  d_ind=a.d_ind;
}

LJException &LJException::operator=(const LJException &b){
  if(this != &b){
    this->~LJException();
    new(this) LJException(b);
  }
  return *this;
}

int gcore::operator<(const LJException &a, const LJException &b){
  if (a[0]<b[0])return 1;
  else if((a[0]==b[0])&&(a[1]<b[1]))return 1;
  else if((a[0]==b[0])&&(a[1]==b[1])&&(a.type()<b.type()))return 1;
  return 0;
}
