// bound_RectBox.cc

#include <cmath>
#include "../gmath/Vec.h"
#include "../gcore/Box.h"
#include "RectBox.h"

using bound::RectBox;
using gmath::Vec;
using gcore::Box;

Vec RectBox::nearestImage(const Vec &r1, const Vec &r2, const Box &box)const {
  Vec diff = r2 - r1;
  Vec a;
  const double kabs = box.K().abs();
  a[0] = diff[0] - kabs * rint(diff[0] / kabs);
  const double labs = box.L().abs();
  a[1] = diff[1] - labs * rint(diff[1] / labs);
  const double mabs = box.M().abs();
  a[2] = diff[2] - mabs * rint(diff[2] / mabs);

  Vec rec = r1 + a;
  return rec;
}


