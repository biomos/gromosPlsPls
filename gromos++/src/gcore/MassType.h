// MassType.h
#ifndef INCLUDED_MASSTYPE
#define INCLUDED_MASSTYPE

namespace gcore{

class MassType
{
  int d_n;
  double d_am;
 public:
  MassType(int n=0, double l=0): d_n(n), d_am(l){}
  MassType(const MassType& b):d_n(b.d_n), d_am(b.d_am){}
  MassType &operator=(const MassType &b);
  ~MassType(){}
  double n()const{return d_n;}
  double am()const{return d_am;}
};

}
#endif
