// gcore_Constraint.h

#ifndef INCLUDED_GCORE_CONSTRAINT
#define INCLUDED_GCORE_CONSTRAINT

namespace gcore{

class Constraint{
  int d_a[2];
  double d_dist;
  // not implemented
  Constraint();
 public:
  Constraint(int a, int b);
  Constraint(const Constraint &);
  ~Constraint(){}
  Constraint &operator=(const Constraint &b);
  void setDist(double a){d_dist=a;}
  int operator[](int i)const{return d_a[i];}
  int &operator[](int i){return d_a[i];}
  double dist()const{return d_dist;}
};

int operator<(const Constraint &a, const Constraint &b);
}
#endif
