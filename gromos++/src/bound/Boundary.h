#ifndef INCLUDED_BOUND_BOUNDARY
#define INCLUDED_BOUND_BOUNDARY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gmath{
  class Vec;
}

namespace gcore{
  class System;
  class Box;
}

using gmath::Vec;

namespace bound{

  class Boundary_i;
  /**
   * Class Boundary
   * Class that defines some basic function for the periodic boundary 
   * conditions
   * 
   * It there more documentation needed?
   * @class Boundary
   * @author R. Buergi, M.K. Kastenholz
   * @ingroup bound
   * @todo finish documentation
   */
  class Boundary {

    Boundary_i *d_this;

    // not implemented
    Boundary& operator=(const Boundary& rhs);
    Boundary(const Boundary& b);
    Boundary();

  public:
    /** Constructor */
    Boundary (gcore::System *sys);
    /** Destructor */
    virtual ~Boundary();

    // Methods

    /**
     * sets reference position for molecule i to v.
     */
    void setReference(int i, const gmath::Vec &v);
    /**
     * sets reference position for all molecules to the first
     * atom of the corresponding molecule in sys
     */
    void setReference(gcore::System const & sys);
    
    /**
     * Given the reference position r1, we give r2 back so that 
     * r2 is the nearest image to r1. Used to reconnect molecules.
     * Note that solvent molecules do never need to be reconnected
     */
    virtual gmath::Vec nearestImage(const gmath::Vec &r1,
				    const  gmath::Vec &r2, 
				    const gcore::Box &box) const = 0;
    
    /**
     * determines whether r is in box.
     */
    //    bool isInBox(const gmath::Vec &r, const gcore::Box &box) const;
   
    /**
     * No gathering
     *
     */
    virtual void nogather(){}; 
    /**
     * gathers the whole System in gromos style (per first molecule).
     */
    virtual void gathergr(){};
    /**
     * gathers the whole system in modified gromos style (per first molecule),
     * but shifts the molecule inside the box if the centre of geometry is outside.
     */
    virtual void gathermgr(){}
    /**
     * gathers the whole system in gromos++ style (per molecule).
     */
    virtual void gather(){};
    /**
     *  gathers solute and solvent with respect to the cog of mol(0)
     */
    virtual void coggather(){};
    /**
     * gathering of e.g. amyloid crystals
     */
    virtual void crsgather(){}; 
    /**
     * same but then gathering of cogs w.r.t overall cog
     */
    virtual void seqgather(){};
    /**
     * attempt for a generalized gathering method (A. Choutko / D. Geerke / A.-P. Kunz)
     */ 
    virtual void gengather(){};
    /**
     * gather by using bonds rather than sequential atoms
     */
    virtual void bondgather(){};
    /**
     * gather by using a reference frame and previous frame
     */
    virtual void refgather(){};
    
    /**
     * reference vector (set to pos(0) of mol(i)) of each molecule upon 
     * creation of object boundary.
     * if the system does not have any coordinates yet, they are initialized
     * with 0.0.
     */
    const gmath::Vec &reference(int i)const;
    /**
     * system accessor.
     */
    gcore::System &sys();
    /**
     * reference system accessor.
     */
    gcore::System &refSys();
    /**
     * the boundary type.
     */
    char type();
    /**
     * set the boundary type.
     */
    void setType(char t);
    void setReferenceFrame(std::string file);
    /**
     * member pointer to gather function
     */
    typedef void (Boundary::*MemPtr)();
  };

}
#endif


