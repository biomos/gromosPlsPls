// gcore_Angle.cc

#include "Angle.h"
#include <new>
#include <iostream>

using gcore::Angle;

Angle::Angle(int a, int b, int c, bool warn){
  d_a[1]=b;
  if(a<c){
    d_a[0]=a;
    d_a[2]=c;
  }
  else{
    d_a[0]=c;
    d_a[2]=a;
    if (warn) {
      std::cerr << "NOTE: order of atoms changed in bond angle:\n";
      std::cerr << "      " << a + 1 << "," << b + 1 << "," << c + 1 << " -> "
              << c + 1 << "," << b + 1 << "," << a + 1 << std::endl;
    }
  }
  d_type=-1;
}

Angle::Angle(const Angle &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_a[2]=a.d_a[2];
  d_type=a.d_type;
}

Angle &Angle::operator=(const Angle &b){
  if(this != &b){
    this->Angle::~Angle();
    new(this) Angle(b);
  }
  return *this;
}

int gcore::operator<(const Angle &a, const Angle &b){
  if (a[1]<b[1])return 1;
  else if (a[1]>b[1])return 0;
  else if (a[0]<b[0])return 1;
  else if (a[0]>b[0])return 0;
  else if (a[2]<b[2])return 1;
  else if (a[2]==b[2]&&a.type()<b.type()) return 1;
  return 0;
}
