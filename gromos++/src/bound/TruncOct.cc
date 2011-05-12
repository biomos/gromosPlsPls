// bound_TruncOct.cc

#include <cmath>
#include <iostream>
#include <sstream>
#include <set>
#include "TruncOct.h"
#include "../gmath/Vec.h"
#include "../gcore/System.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Bond.h"
#include "../gcore/Box.h"

using bound::TruncOct;
using gmath::Vec;
using gcore::Box;

Vec TruncOct::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const{
  Vec diff=r2-r1;
  Vec a;

  const double kabs = box.K().abs();
  a[0] = diff[0] - kabs * rint(diff[0]/kabs);
  a[1] = diff[1] - kabs * rint(diff[1]/kabs);
  a[2] = diff[2] - kabs * rint(diff[2]/kabs);

  if ( (0.75*kabs - fabs(a[0]) - fabs(a[1]) - fabs(a[2])) < 0.0) {
    const double half_kabs = 0.5 * kabs;
    a[0] = a[0] - a[0]/fabs(a[0])*half_kabs;
    a[1] = a[1] - a[1]/fabs(a[1])*half_kabs;
    a[2] = a[2] - a[2]/fabs(a[2])*half_kabs;

  }

  return r1 + a;
}
