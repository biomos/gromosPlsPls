// AngleType.h
#ifndef INCLUDED_ANGLETYPE
#define INCLUDED_ANGLETYPE

namespace gcore{

class AngleType
{
  double d_t0;
  double d_fc;
 public:
  AngleType(double fc=0, double l=0): d_t0(l), d_fc(fc){}
  AngleType(const AngleType& b):d_t0(b.d_t0), d_fc(b.d_fc){}
  AngleType &operator=(const AngleType &b);
  ~AngleType(){}
  double t0()const{return d_t0;}
  double fc()const{return d_fc;}
};

}
#endif
