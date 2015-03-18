// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRUNCOCT
#define INCLUDED_BOUND_TRUNCOCT

#include "Boundary.h"
#include "../gromos/Exception.h"

namespace bound{
  /**
   * Class TruncOct
   * Class that defines periodic boundary conditions for a truncated 
   * octahedral box
   *
   * @class TruncOct
   * @author R. Buergi
   * @ingroup bound
   */
  class TruncOct: public Boundary {
  public:
    // Constructor
    TruncOct(gcore::System *sys): Boundary(sys){setType('t');}
    virtual ~TruncOct(){}
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
  };
    
}

#endif
