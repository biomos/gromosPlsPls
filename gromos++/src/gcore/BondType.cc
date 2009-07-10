// gcore_BondType.cc
#include "BondType.h"

#include <new>

using gcore::BondType;

BondType &BondType::operator=(const BondType &b) {
  if (this != &b) {
    this->~BondType();
    new(this) BondType(b);
  }
  return *this;
}

BondType::BondType(int c, double fc, double l) : d_code(c), d_fc(fc), d_b0(l) {
  d_hfc = 2.0 * d_b0 * d_b0 * d_fc;
}
