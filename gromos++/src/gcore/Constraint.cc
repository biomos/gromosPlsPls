// gcore_Constraint.cc

#include "Constraint.h"
#include <new>

using gcore::Constraint;

Constraint::Constraint(int a, int b){
  if(a<b){
    d_a[0]=a;
    d_a[1]=b;
  }
  else{
    d_a[0]=b;
    d_a[1]=a;
  }
  d_dist=-1;
}

Constraint::Constraint(const Constraint &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_dist=a.d_dist;
}

Constraint &Constraint::operator=(const Constraint &b){
  if(this != &b){
    this->~Constraint();
    new(this) Constraint(b);
  }
  return *this;
}

int gcore::operator<(const Constraint &a, const Constraint &b){
  if (a[0]<b[0])return 1;
  else if((a[0]==b[0])&&(a[1]<b[1]))return 1;
  return 0;
}
