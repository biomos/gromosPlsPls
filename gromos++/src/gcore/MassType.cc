// gcore_MassType.cc
#include "MassType.h"

#include <new>

using gcore::MassType;

MassType &MassType::operator=(const MassType &b){
  if(this!=&b){
    this->~MassType();
    new(this) MassType(b);
  }
  return *this;
}

