// gcore_Bond.h

#ifndef INCLUDED_GCORE_BOND
#define INCLUDED_GCORE_BOND

namespace gcore{

class Bond{
  int d_a[2];
  int d_type;
  // not implemented
  Bond();
 public:
  Bond(int a, int b);
  Bond(const Bond &);
  ~Bond(){}
  Bond &operator=(const Bond &b);
  void setType(int i){d_type=i;}
  int operator[](int i)const{return d_a[i];}
  int &operator[](int i){return d_a[i];}
  int type()const{return d_type;}
};

int operator<(const Bond &a, const Bond &b);
}
#endif
