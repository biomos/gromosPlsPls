// gcore_DihedralType.cc
#include "DihedralType.h"

#include <new>

using gcore::DihedralType;

DihedralType &DihedralType::operator=(const DihedralType &b){
  if(this!=&b){
    this->~DihedralType();
    new(this) DihedralType(b);
  }
  return *this;
}
