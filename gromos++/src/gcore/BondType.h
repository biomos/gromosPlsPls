// BondType.h
#ifndef INCLUDED_BONDTYPE
#define INCLUDED_BONDTYPE

namespace gcore{

class BondType
{
  double d_fc;
  double d_b0;
 public:
  BondType(double fc=0, double l=0): d_fc(fc), d_b0(l){}
  BondType(const BondType& b):d_fc(b.d_fc), d_b0(b.d_b0){}
  BondType &operator=(const BondType &b);
  ~BondType(){}
  double b0()const{return d_b0;}
  double fc()const{return d_fc;}
};

}
#endif
