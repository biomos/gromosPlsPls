// bound_RectBox.h

#ifndef INCLUDED_BOUND_RECTBOX
#define INCLUDED_BOUND_RECTBOX

#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif


namespace gmath{
class Vec;
}

using gmath::Vec;

namespace bound{
  /**
   * Class RectBox
   * Defines the periodic boundary conditions for a rectangular box
   *
   * @class RectBox
   * @author M.K. Kastenholz
   * @ingroup bound
   * @todo finish documentation
   */
  class RectBox: public Boundary {
    
    // not implemented
    RectBox(const RectBox &);
    RectBox &operator=(const RectBox &);
    RectBox();
  public:
    // Constructor
    RectBox(gcore::System *sys): Boundary(sys){setType('r');}
    ~RectBox(){}
    
    gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const;
   
    virtual void nogather(); 
    virtual void gathergr();
    virtual void gathermgr();
    virtual void gather();
    virtual void coggather();
    virtual void seqgather();
    virtual void crsgather();
    virtual void gengather();
  };
    
}

#endif
