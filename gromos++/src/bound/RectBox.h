// bound_RectBox.h

#ifndef INCLUDED_BOUND_RECTBOX
#define INCLUDED_BOUND_RECTBOX

#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace bound {
  /**
   * Class RectBox
   * Defines the periodic boundary conditions for a rectangular box
   *
   * @class RectBox
   * @author M.K. Kastenholz
   * @ingroup bound
   */
  class RectBox: public Boundary {
  public:
    // Constructor
    RectBox(gcore::System *sys): Boundary(sys){setType('r');}
    virtual ~RectBox(){}
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
  };    
}

#endif
