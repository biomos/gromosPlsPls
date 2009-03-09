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

  public:
    enum boxshape_enum { vacuum=0, rectangular=1, triclinic=2, truncoct=3 };
    enum boxformat_enum { box96, triclinicbox, genbox};

  private:
    std::vector<gmath::Vec> d_dim;
    double d_K_L_M;
    std::vector<gmath::Vec> d_cross_K_L_M;
    boxshape_enum  d_ntb;
    boxformat_enum d_boxformat;

  public:
    // constructurs
    /**
     * Box constructor
     * @param x, y, z box dimensions
     */
    Box(double x=0, double y=0, double z=0):
      d_dim(3),d_K_L_M(0),d_cross_K_L_M(3), d_ntb(vacuum), d_boxformat(genbox){
      d_dim[0][0]=x; d_dim[1][1]=y; d_dim[2][2]=z;
    }

    Box(gmath::Vec K, gmath::Vec L, gmath::Vec M):
      d_dim(3),d_K_L_M(0),d_cross_K_L_M(3), d_ntb(vacuum), d_boxformat(genbox)
    {
      d_dim[0]=K;
      d_dim[1]=L;
      d_dim[2]=M;
    }
    

    /**
     * Box copy constructor
     * @param b Box to be copied
     */
    Box(const Box&b):
      d_dim(b.d_dim), d_K_L_M(b.d_K_L_M), 
      d_cross_K_L_M(b.d_cross_K_L_M), d_ntb(b.d_ntb),
      d_boxformat(b.d_boxformat)
    {}
    /**
     * Update volume and cross product for the triclinic box
     */
    void update_triclinic();
    
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
    gmath::Vec const & K()const;
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec &L();
    /**
     * Accessor, return the L-vector of a generalized box
     */
    gmath::Vec const & L()const;
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec &M();
    /**
     * Accessor, return the M-vector of a generalized box
     */
    gmath::Vec const & M()const;
    /**
     * Accessor, return the value of ntb
     */
    void setNtb(boxshape_enum b);
    /**
     * Accessor, return the vlaue of ntb
     */
    boxshape_enum ntb()const;
    /**
     * Accessor, return the format of the box
     */
    boxformat_enum &boxformat();
    /**
     * Accessor, return the format of the box
     */
    boxformat_enum boxformat()const;
    
    /**
     * Accessor, return the (triclinic) volume of a generalized box
     */
    double K_L_M();
    /**
     * Accessor, return the (triclinic) volume of a generalized box
     */
    double K_L_M()const;
        /**
     * Accessor, return the weird cross product of a generalized box
     */
    std::vector<gmath::Vec> & cross_K_L_M();
    /**
     * Accessor, return the weird cross product of a generalized box
     */
    std::vector<gmath::Vec> const & cross_K_L_M()const;
    
  };

  inline gcore::Box &Box::operator=(gcore::Box const &b)
  {
    if(this!=&b){
      d_dim=b.d_dim;
      d_K_L_M = b.d_K_L_M;
      d_cross_K_L_M = b.d_cross_K_L_M;
      d_ntb = b.d_ntb;
      d_boxformat = b.d_boxformat;
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
  inline gmath::Vec const & Box::K()const
  {
    return d_dim[0];
  }
  inline gmath::Vec &Box::L()
  {
    return d_dim[1];
  }
  inline gmath::Vec const & Box::L()const
  {
    return d_dim[1];
  }
  inline gmath::Vec &Box::M()
  {
    return d_dim[2];
  }
  inline gmath::Vec const & Box::M()const
  {
    return d_dim[2];
  }
  inline Box::boxshape_enum Box::ntb()const
  {
    return d_ntb;
  }
  inline Box::boxformat_enum & Box::boxformat()
  {
    return d_boxformat;
  }
  inline Box::boxformat_enum Box::boxformat()const
  {
    return d_boxformat;
  }
  inline double Box::K_L_M()
  {
    return d_K_L_M;
  }
  inline double Box::K_L_M()const
  {
    return d_K_L_M;
  }
  inline std::vector<gmath::Vec> & Box::cross_K_L_M()
  {
    return d_cross_K_L_M;
  }
  inline std::vector<gmath::Vec> const & Box::cross_K_L_M()const
  {
    return d_cross_K_L_M;
  }

    
} /*namespace*/
