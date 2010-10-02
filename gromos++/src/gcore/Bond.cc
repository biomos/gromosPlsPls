// gcore_Bond.cc

#include "Bond.h"
#include <new>
#include <iostream>

using gcore::Bond;

Bond::Bond(int a, int b, bool warn){
  if(a<b){
    d_a[0]=a;
    d_a[1]=b;
  }
  else{
    d_a[0]=b;
    d_a[1]=a;
    if(warn) {
      std::cerr << "NOTE: order of atoms changed in bond:\n";
      std::cerr << "      " << a+1 << "," << b+1 << " -> " << b+1 << "," << a+1 << std::endl;
    }
  }
  d_type=-1;
}

Bond::Bond(const Bond &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_type=a.d_type;
}

Bond &Bond::operator=(const Bond &b){
  if(this != &b){
    this->~Bond();
    new(this) Bond(b);
  }
  return *this;
}

int gcore::operator<(const Bond &a, const Bond &b){
  if (a[0]<b[0])return 1;
  else if((a[0]==b[0])&&(a[1]<b[1]))return 1;
  else if((a[0]==b[0])&&(a[1]==b[1])&&(a.type()<b.type()))return 1;
  return 0;
}
