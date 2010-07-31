// gcore_LJExceptionType.cc
#include "LJExceptionType.h"

#include <new>

using gcore::LJExceptionType;

LJExceptionType &LJExceptionType::operator=(const LJExceptionType &b) {
  if (this != &b) {
    this->~LJExceptionType();
    new(this) LJExceptionType(b);
  }
  return *this;
}
