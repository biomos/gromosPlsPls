// gcore_ImproperType.cc
#include "ImproperType.h"

#include <new>

using gcore::ImproperType;

ImproperType &ImproperType::operator=(const ImproperType &b){
  if(this!=&b){
    this->~ImproperType();
    new(this) ImproperType(b);
  }
  return *this;
}
