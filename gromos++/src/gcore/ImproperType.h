// ImproperType.h
#ifndef INCLUDED_IMPROPERTYPE
#define INCLUDED_IMPROPERTYPE

namespace gcore{

class ImproperType{
  double d_q0;
  double d_fc;
 public:
  ImproperType(double fc=0, double l=0): d_q0(l), d_fc(fc){}
  ImproperType(const ImproperType& b):d_q0(b.d_q0), d_fc(b.d_fc){}
  ~ImproperType(){}
  ImproperType &operator=(const ImproperType &);
  double q0()const{return d_q0;}
  double fc()const{return d_fc;}
};

}
#endif
