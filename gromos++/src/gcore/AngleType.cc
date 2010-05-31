// gcore_AngleType.cc
#include <iostream>
#include <cmath>
#include "AngleType.h"
#include <new>
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../gmath/Physics.h"

using gcore::AngleType;

AngleType &AngleType::operator=(const AngleType &b){
  if(this!=&b){
    this->~AngleType();
    new(this) AngleType(b);
  }
  return *this;
}

AngleType::AngleType(int c, double fc, double l) : d_code(c), d_t0(l),d_fc(fc),
        d_afc(0.0) {
  // we only got the harmonic-in-cosine force constant and have to calculate the
  // harmonic-in-angle force constant.
  // we do this as discribed in volume 3
  
  // first let's make sure it is not simply zero.
  if (d_fc == 0.0) {
    return;
  }
  
  double t_0 = d_t0 * gmath::physConst.get_degree2radian();
  
  // this is kT at 300K - a reasonable choice for biomolecular simulations.
  const double kT = gmath::physConst.get_boltzmann() * 300.0;
  
  const double X = cos(t_0) * cos(t_0);
  const double Y = (d_fc - kT) / d_fc;
  const double Z = X - sin(t_0) * sin(t_0);
  
  const double under_sqrt = X*X - Y*Z;
  if (under_sqrt < 0 || Z == 0.0)
    throw gromos::Exception("AngleType", "Cannot automatically convert "
            "force constant: discriminant negative.");

  const double cos_angle = (X - sqrt(under_sqrt)) / Z;
  if (cos_angle < -1.0 || cos_angle > 1.0)
    throw gromos::Exception("AngleType", "Cannot automatically convert "
            "force constant: cosine out of range.");
  
  const double angle = acos(cos_angle);
  if (angle == 0.0)
    throw gromos::Exception("AngleType", "Cannot automatically convert "
            "force constant: angle is zero.");
  
  d_afc = kT / (angle * angle);
  
  // do a warning as this may be a bit tricky
  if (!args::Arguments::outG96) {
    std::cerr << "Warning for bond angle type " << d_code + 1 << ": Harmonic force "
            "constant was calculated at 300K." << std::endl;
    std::cerr << "WARNING: check the calculated value before running the simulation!!!\n";
  }
}
