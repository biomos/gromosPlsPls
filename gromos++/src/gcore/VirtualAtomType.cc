// gcore_VirtualAtomType.cc
#include "VirtualAtomType.h"

#include <new>

using gcore::VirtualAtomType;

VirtualAtomType &VirtualAtomType::operator=(const VirtualAtomType &va) {
  if (this != &va) {
    this->~VirtualAtomType();
    new(this) VirtualAtomType(va);
  }
  return *this;
}


