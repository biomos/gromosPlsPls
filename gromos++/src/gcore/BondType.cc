// gcore_BondType.cc
#include "BondType.h"

#include <new>

using gcore::BondType;

BondType &BondType::operator=(const BondType &b){
  if(this!=&b){
    this->~BondType();
    new(this) BondType(b);
  }
  return *this;
}

