// gmath_Vec.cc

#include <iomanip>
#include <sstream>

#include "Vec.h"

std::string gmath::v2s(gmath::Vec const &v)
{
  std::ostringstream os;
  os << std::setw(20) << v[0]
     << std::setw(20) << v[1]
     << std::setw(20) << v[2];
  return os.str();
}



