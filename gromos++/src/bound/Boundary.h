#ifndef INCLUDED_BOUND_BOUNDARY
#define INCLUDED_BOUND_BOUNDARY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gmath{
  class Vec;
}

namespace gcore{
  class System;
  class Box;
}

using gmath::Vec;

namespace bound{

  class Boundary_i;
  /**
   * Class Boundary
   * Class that defines some basic function for the periodic boundary 
   * conditions
   * 
   * It there more documentation needed?
   * @class Boundary
   * @author R. Buergi, M.K. Kastenholz
   * @ingroup bound
   * @todo finish documentation
   */
  class Boundary {

    Boundary_i *d_this;

    // not implemented
    Boundary& operator=(const Boundary& rhs);
    Boundary(const Boundary& b);
    Boundary();

  public:
    // Constructor
    Boundary (gcore::System *sys);
    virtual ~Boundary();

    // Methods

    // sets reference for molecule i to v.
    void setReference(int i, const gmath::Vec &v);

    /** Given the reference position r1, we give r2 back so that 
	r2 is the nearest image to r1. Used to reconnect molecules.
	Note that solvent molecules do never need to be reconnected
    */
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
				    const  gmath::Vec &r2, 
				    const gcore::Box &box) const = 0;

    // gathers the whole System...
   virtual void gathergr(){};
   virtual void gather(){};
    // gathers solute and solvent with respect to the cog of mol(0)
   virtual void coggather(){};

    // reference vector (set to pos(0) of mol(i)) of each molecule upon 
    // creation of object boundary
    const gmath::Vec &reference(int i)const;
    gcore::System &sys();
    
    char type();
    void setType(char t);
    
    typedef void (Boundary::*MemPtr)();
  };

}
#endif


