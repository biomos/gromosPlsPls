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
  /**
   * Class Vacuum
   * defines the periodic boundary conditions for a vacuum. Which means 
   * that there are no periodic boundary conditions
   *
   * @class Vacuum
   * @author R. Buergi
   * @ingroup bound
   * @todo finish documentation
   */
  class Vacuum: public Boundary {

    // not implemented
    Vacuum();
    Vacuum(const Vacuum &);
    Vacuum &operator=(const Vacuum &);

  public:
    // ------  CONSTRUCTORS
    // default constructor
    Vacuum(gcore::System *sys): Boundary(sys) {setType('v');}
    ~Vacuum(){}
    
    // ------  METHODS
    
    gmath::Vec nearestImage(const gmath::Vec &r1,
			    const  gmath::Vec &r2, 
			    const gcore::Box &box) const
      { return r2;}
 
    void nogather(){};  
    void gathergr(){}; 
    void gather(){};
    void coggather(){};
    void seqgather(){};
    void crsgather(){};
    void gengather(){};
    void bondgather(){};
  };

}

#endif
