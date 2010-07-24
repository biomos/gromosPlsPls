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
   
    void nogather(); 
    void gathergr();
    void gathermgr();
    void gather();
    void coggather();
    void seqgather();
    void crsgather();
    void gengather();
    void bondgather();
    void refgather();
    // DW : add new methods
    // gather based on a general list
    void gatherlist();
    // gather in term of time
    void gathertime();
    // gather based on a reference structure
    void gatherref();
    // gather the first frame based on an atom list, then the rest in term of time
    void gatherltime();
    // gather the first frame based on a reference structure, then the rest in term of time
    void gatherrtime();
    // gather based on bond connection
    void gatherbond();
  };
    
}

#endif
