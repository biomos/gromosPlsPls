// gcore_Angle.h

#ifndef INCLUDED_GCORE_ANGLE
#define INCLUDED_GCORE_ANGLE

namespace gcore{

class Angle{
  int d_a[3];
  int d_type;
  // not implemented
  Angle();
 public:
  Angle(int a, int b, int c);
  Angle(const Angle &);
  ~Angle(){}
  void setType(int i){d_type=i;}
  Angle &operator=(const Angle &);
  int operator[](int i)const{return d_a[i];}
  int &operator[](int i){return d_a[i];}
  int type()const{return d_type;}
};

int operator<(const Angle &a, const Angle &b);

}
#endif
