// bound_Vacuum.h

#ifndef INCLUDED_BOUND_VACUUM
#define INCLUDED_BOUND_VACUUM


#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

#ifndef INCLUDED_GMATH_VEC
#include "../gmath/Vec.h"
#endif

namespace bound{

  class Vacuum: public Boundary {

    // not implemented
    Vacuum();
    Vacuum(const Vacuum &);
    Vacuum &operator=(const Vacuum &);

  public:
    // ------  CONSTRUCTORS
    // default constructor
    Vacuum(gcore::System *sys): Boundary(sys) {}
    ~Vacuum(){}
    
    // ------  METHODS
    
    gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const
      { return r2;}
   
    void gathergr(){}; 
    void gather(){}
    void coggather(){}
  };

}

#endif
