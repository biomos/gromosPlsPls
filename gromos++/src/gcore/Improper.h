// gcore_Improper.h

#ifndef INCLUDED_GCORE_IMPROPER
#define INCLUDED_GCORE_IMPROPER

namespace gcore{

class Improper{
  int d_a[4];
  int d_type;
  // not implemented
  Improper();
 public:
  Improper(int a, int b, int c, int d);
  Improper(const Improper &);
  ~Improper(){}
  void setType(int i){d_type=i;}
  Improper &operator=(const Improper &);
  int operator[](int i)const{return d_a[i];}
  int &operator[](int i){return d_a[i];}
  int type()const{return d_type;}
};

int operator<(const Improper &a, const Improper &b);

}
#endif
