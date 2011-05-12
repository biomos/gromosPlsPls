// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRICLINIC
#define INCLUDED_BOUND_TRICLINIC

#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace bound{
  /**
   * Class Triclinic
   * Class that defines periodic boundary conditions for a triclinic 
   * box
   *
   * @class Triclinic
   * @author R. Buergi
   * @ingroup bound
   */
  class Triclinic: public Boundary {
  public:
    // Constructor
    Triclinic(gcore::System *sys): Boundary(sys){setType('c');}
    virtual ~Triclinic(){}
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
  };
    
}

#endif
