// bound_Triclinic.cc

#include <cmath>
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include "Triclinic.h"

using bound::Triclinic;
using gmath::Vec;
using gcore::Box;

Vec Triclinic::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const{
  Vec P = r2 - r1;
  int k,l,m;
  k = int(rint(box.cross_K_L_M()[0].dot(P)));
  l = int(rint(box.cross_K_L_M()[1].dot(P)));
  m = int(rint(box.cross_K_L_M()[2].dot(P)));
  
  P += box.K() * k + box.L() * l + box.M() * m;
  return r1 + P;
}
