// gcore_Dihedral.h

#ifndef INCLUDED_GCORE_DIHEDRAL
#define INCLUDED_GCORE_DIHEDRAL

namespace gcore{

class Dihedral{
  int d_a[4];
  int d_type;
  // not implemented
  Dihedral();
 public:
  Dihedral(int a, int b, int c, int d);
  Dihedral(const Dihedral &);
  ~Dihedral(){}
  Dihedral &operator=(const Dihedral &);
  void setType(int i){d_type = i;}
  int operator[](int i)const{return d_a[i];}
  int &operator[](int i){return d_a[i];}
  int type()const{return d_type;}
};

int operator<(const Dihedral &a, const Dihedral &b);

}
#endif
