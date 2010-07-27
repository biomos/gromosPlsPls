// gcore_LJExcType.cc

#include "LJExcType.h"
#include <new>

using gcore::LJExcType;

LJExcType &LJExcType::operator=(const LJExcType &l){
  if(this!=&l){
    this->~LJExcType();
    new(this) LJExcType(l);
  }
  return *this;
}


