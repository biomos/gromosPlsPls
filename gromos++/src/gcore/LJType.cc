// gcore_LJType.cc

#include "LJType.h"
#include <new>

using gcore::LJType;

LJType &LJType::operator=(const LJType &l){
  if(this!=&l){
    this->~LJType();
    new(this) LJType(l);
  }
  return *this;
}


