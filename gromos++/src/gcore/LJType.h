// LJType.h
#ifndef INCLUDED_LJTYPE
#define INCLUDED_LJTYPE

namespace gcore{

class LJType
{
  double d_c12, d_c6, d_cs12, d_cs6;
  public:
  LJType(double d1=0, double d2=0, double d3=0, double d4=0):
    d_c12(d1), d_c6(d2), d_cs12(d3), d_cs6(d4){}
  LJType(const LJType &l): d_c12(l.d_c12), d_c6(l.d_c6), 
    d_cs12(l.d_cs12), d_cs6(l.d_cs6){}
  LJType &operator=(const LJType &l);
  ~LJType(){}
  double c12()const{return d_c12;}
  double c6()const {return d_c6;}
  double cs12()const {return d_cs12;}
  double cs6()const {return d_cs6;}
};

}
#endif
