// gcore_AtomPair.cc

#include "AtomPair.h"
#include <new>

using gcore::AtomPair;

AtomPair::AtomPair(int a, int b){
  if(a<b){
    d_a[0]=a;
    d_a[1]=b;
  }
  else{
    d_a[0]=b;
    d_a[1]=a;
  }
}

AtomPair::AtomPair(const AtomPair &a){
    d_a[0]=a.d_a[0];
    d_a[1]=a.d_a[1];
}

int gcore::operator<(const AtomPair &a, const AtomPair &b){
  if (a[0]<b[0])return 1;
  else if((a[0]==b[0])&&(a[1]<b[1]))return 1;
  return 0;
}
