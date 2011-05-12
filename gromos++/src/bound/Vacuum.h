// bound_Vacuum.h

#ifndef INCLUDED_BOUND_VACUUM
#define INCLUDED_BOUND_VACUUM


#ifndef INCLUDED_BOUND_BOUNDARY
#include "Boundary.h"
#endif

#ifndef INCLUDED_GMATH_VEC
#include "../gmath/Vec.h"
#endif

namespace bound {

  /**
   * Class Vacuum
   * defines the periodic boundary conditions for a vacuum. Which means 
   * that there are no periodic boundary conditions
   *
   * @class Vacuum
   * @author R. Buergi
   * @ingroup bound
   */
  class Vacuum : public Boundary {
  public:
    Vacuum(gcore::System *sys) : Boundary(sys) {
      setType('v');
    }

    virtual ~Vacuum() {
    }

    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
            const gmath::Vec &r2,
            const gcore::Box &box) const {
      return r2;
    }

    // overwrite gathering methods as they do not make sense for vacuum
    virtual void nogather() {
    }

    virtual void gathergr() {
    }

    virtual void gather() {
    }

    virtual void coggather() {
    }

    virtual void seqgather() {
    }

    virtual void crsgather() {
    }

    virtual void gengather() {
    }

    virtual void bondgather() {
    }

    virtual void refgather() {
    }

    virtual void gatherlist() {
    }

    virtual void gathertime() {
    }

    virtual void gatherref() {
    }

    virtual void gatherltime() {
    }

    virtual void gatherrtime() {
    }

    virtual void gatherbond() {
    }
  };

}

#endif
