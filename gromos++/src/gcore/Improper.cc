// gcore_Improper.cc

#include "Improper.h"
#include <new>
#include<iostream>

using gcore::Improper;

Improper::Improper(int a, int b, int c, int d, bool warn){
  if(b<c){
    d_a[0]=a; 
    d_a[1]=b; 
    d_a[2]=c; 
    d_a[3]=d;
  } 
  else{
    d_a[0]=d;
    d_a[1]=c;
    d_a[2]=b;
    d_a[3] = a;
    if (warn) {
      std::cerr << "NOTE: order of atoms changed in improper dihedral:\n";
      std::cerr << "      " << a + 1 << "," << b + 1 << "," << c + 1 << "," << d + 1 << " -> "
              << d + 1 << "," << c + 1 << "," << b + 1 << "," << a + 1 << std::endl;
    }
  }
  d_type = -1;
}

Improper::Improper(const Improper &a){
  d_a[0]=a.d_a[0];
  d_a[1]=a.d_a[1];
  d_a[2]=a.d_a[2];
  d_a[3]=a.d_a[3];
  d_type=a.d_type;
}

Improper &Improper::operator=(const Improper &b){
  if(this != &b){
    this->Improper::~Improper();
    new(this) Improper(b);
  }
  return *this;
}

int gcore::operator<(const Improper &a, const Improper &b){
return (a[0]<b[0]||(a[0]==b[0]&&a[1]<b[1])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]<b[2])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]<b[3])
	||(a[0]==b[0]&&a[1]==b[1]&&a[2]==b[2]&&a[3]==b[3]&&a.type()<b.type()));
}
