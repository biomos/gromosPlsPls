// gcore_AngleType.cc
#include "AngleType.h"
#include <new>

using gcore::AngleType;

AngleType &AngleType::operator=(const AngleType &b){
  if(this!=&b){
    this->~AngleType();
    new(this) AngleType(b);
  }
  return *this;
}
