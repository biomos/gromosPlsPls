// gcore_Box.h

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif

namespace gcore{
  class Box{
    double d_dim[3];
    double d_angle;
    // not implemented
    Box &operator=(const Box&);
  public:
    // constructurs
    Box(double x=0, double y=0, double z=0, double a=90){
      d_dim[0]=x; d_dim[1]=y; d_dim[2]=z;
      d_angle=a;
    }
    Box(const Box&b):
      d_angle(b.d_angle){
      d_dim[0]=b.d_dim[0];
      d_dim[1]=b.d_dim[1];
      d_dim[2]=b.d_dim[2];
    }
    // accessors
    double &operator[](int i);
    double operator[](int i)const;
    double angle()const;
  };

  inline double &Box::operator[](int i){
    assert (i<3);
    return d_dim[i];
  }
  inline double Box::operator[](int i)const{
    assert (i<3);
    return d_dim[i];
  }
  inline double Box::angle()const{
    return d_angle;
  }

} /*namespace*/
