// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRICLINIC
#define INCLUDED_BOUND_TRICLINIC

#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

namespace gmath{
class Vec;
}

using gmath::Vec;


namespace bound{
  /**
   * Class Triclinic
   * Class that defines periodic boundary conditions for a triclinic 
   * box
   *
   * @class Triclinic
   * @author R. Buergi
   * @ingroup bound
   * @todo finish documentation
   */
  class Triclinic: public Boundary {
    
    // not implemented
    Triclinic(const Triclinic &);
    Triclinic &operator=(const Triclinic &);
    Triclinic();
  public:
    // Constructor
    Triclinic(gcore::System *sys): Boundary(sys){setType('c');}
    ~Triclinic(){}
    
    gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
   
    void nogather(); 
    void gather();
    void gathergr();
    void coggather();
    void gengather();
  };
    
}

#endif
