// gmath_Vec.t.cc

#include "Vec.h"
#include <iostream>

using gmath::Vec;
using namespace std;

ostream &operator<<(ostream &o, const Vec &v){
  o << '(' << v[0] << ' '
    << v[1] << ' '
    << v[2] << ')';
  return o;
}

int main(){
  Vec v(1,2,3);
  Vec w(4,5,6);
  cout << "v = "<< v << " w = " << w << endl;
  cout << "v + w = "<< v+w << endl;
  cout << "v - w = "<< v-w << endl;
  cout << "v * 3 = "<< v*3 << endl;
  cout << "3 * v = "<< 3*v << endl;
  cout << "v / 2 = " << v/2 << endl;
  w-=v;
  cout << "w-=v: " << w << endl;
  w+=v;
  cout << "w+=v: " << w << endl;
  cout << "v.dot(w) = " << v.dot(w) << endl;
  cout << "v.cross(w) = " << v.cross(w) << endl;
  Vec t;
  t=w+v;
  cout << "t=w+v: " << t << endl;
  cout << "v.abs2() = " << v.abs2() << endl;
  cout << "v.abs() = " << v.abs() << endl;
  return 0;
}

