// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRUNCOCT
#define INCLUDED_BOUND_TRUNCOCT

#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

namespace gmath{
class Vec;
}

using gmath::Vec;


namespace bound{
  /**
   * Class TruncOct
   * Class that defines periodic boundary conditions for a truncated 
   * octahedral box
   *
   * @class TruncOct
   * @author R. Buergi
   * @ingroup bound
   * @todo finish documentation
   */
  class TruncOct: public Boundary {
    
    // not implemented
    TruncOct(const TruncOct &);
    TruncOct &operator=(const TruncOct &);
    TruncOct();
  public:
    // Constructor
    TruncOct(gcore::System *sys): Boundary(sys){}
    ~TruncOct(){}
    
    gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
    
    void gather();
    void gathergr();
    void coggather();
  };
    
}

#endif
