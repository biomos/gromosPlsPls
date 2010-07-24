// bound_TruncOct.h

#ifndef INCLUDED_BOUND_TRICLINIC
#define INCLUDED_BOUND_TRICLINIC

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
   * Class Triclinic
   * Class that defines periodic boundary conditions for a triclinic 
   * box
   *
   * @class Triclinic
   * @author R. Buergi
   * @ingroup bound
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
    void seqgather();
    void crsgather();
    void gengather();
    void bondgather();
    void refgather();
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
