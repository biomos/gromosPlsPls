// gcore_Box.h

#ifndef INCLUDED_CASSERT
#include <cassert>
#define INCLUDED_CASSERT
#endif
#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_GMATH_VEC
#include "../gmath/Vec.h"
#define INCLUDED_GMATH_VEC
#endif

namespace gmath
{
  class Vec;
}

namespace gcore{
  /**
   * Class Box
   * Purpose: contains the box size and angle
   *
   * Description:
   * The Box class contains the box dimensions and angle
   *
   * @class Box
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::System
   */
  class Box{
    std::vector<gmath::Vec> d_dim;
  public:
    // constructurs
    /**
     * Box constructor
     * @param x, y, z box dimensions
     */
    Box(double x=0, double y=0, double z=0){
      d_dim.resize(3);
      d_dim[0][0]=x; d_dim[1][1]=y; d_dim[2][2]=z;
    }

    Box(gmath::Vec K, gmath::Vec L, gmath::Vec M)
    {
      d_dim.push_back(K);
      d_dim.push_back(L);
      d_dim.push_back(M);
    }
    

    /**
     * Box copy constructor
     * @param b Box to be copied
     */
    Box(const Box&b):
      d_dim(b.d_dim){}

    // accessors
    /**
     * Accessor, returns the box size in the i-th dimension
     */
    double &operator[](int i);
    /**
     * Accessor, returns the box size in the i-th dimension as a const
     */
    double operator[](int i)const;
    /**
     * assignment operator
     */
    Box &operator=(const Box&);

    /**
     * Accessor, return the K-vector of a generalized box
     */
    gmath::Vec &K();
    /**
     * Accessor, return the K-vector of a generalized box
     */
    gmath::Vec K()const;
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec &L();
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec L()const;
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec &M();
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec M()const;
    
  };

  inline gcore::Box &Box::operator=(gcore::Box const &b)
  {
    if(this!=&b){
      d_dim=b.d_dim;
    }
    return *this;
  }
  
  inline double &Box::operator[](int i){
    assert (i<3);
    return d_dim[i][i];
  }
  inline double Box::operator[](int i)const{
    assert (i<3);
    return d_dim[i][i];
  }

  inline gmath::Vec &Box::K()
  {
    return d_dim[0];
  }
  inline gmath::Vec Box::K()const
  {
    return d_dim[0];
  }
  inline gmath::Vec &Box::L()
  {
    return d_dim[1];
  }
  inline gmath::Vec Box::L()const
  {
    return d_dim[1];
  }
  inline gmath::Vec &Box::M()
  {
    return d_dim[2];
  }
  inline gmath::Vec Box::M()const
  {
    return d_dim[2];
  }
} /*namespace*/
