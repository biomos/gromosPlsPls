#ifndef INCLUDED_BOUND_BOUNDARY
#define INCLUDED_BOUND_BOUNDARY

#include <string>
#include <vector>

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
   * conditions. It provides a virtual function for the nearest image
   * calculation and virtual functions to the specific gathering methods.
   * Usually it is constructed by the @ref args::BoundaryParser class.
   * 
   * It there more documentation needed?
   * @class Boundary
   * @author R. Buergi, M.K. Kastenholz
   * @ingroup bound
   */
  class Boundary {

    Boundary_i *d_this;
    std::vector<int > d_refmol;
    bool firstframe;

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
     * No gathering
     */
    virtual void nogather(); 
    /**
     * gather based on a general list. the atom pair should be in the sequence: A B,
     * where A is an atom of the molecule to be gathered, and B is an atom of the
     * reference molecule.
     */
    virtual void gatherlist();
    /**
     * gather in term of time
     */
    virtual void gathertime();
    /**
     * gather based on a reference structure
     */
    virtual void gatherref();
    /**
     * gather the first frame based on an atom list, then the rest in term of time
     */
    virtual void gatherltime();
    /**
     * gather the first frame based on a reference, then the rest in term of time
     */
    virtual void gatherrtime();
    /**
     * gather based on bond connection
     */
    virtual void gatherbond();

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
    /**
     * set the path of to a reference frame that is used in the
     * reference gathering method
     */
    void setReferenceFrame(std::string file);
    /**
     * set the reference system to sys
     */
    void setReferenceSystem(gcore::System system);
    /**
     * add molecules to be gathered by reference
     */
    void addRefMol(int molnum);
    /**
     * member pointer to gather function
     */
    typedef void (Boundary::*MemPtr)();
    
    // old/deprecated gathering methods
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gathers the whole System in gromos style (per first molecule).
     */
    virtual void gathergr();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gathers the whole system in modified gromos style (per first molecule),
     * but shifts the molecule inside the box if the centre of geometry is outside.
     */
    virtual void gathermgr();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gathers the whole system in gromos++ style (per molecule).
     */
    virtual void gather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gathers solute and solvent with respect to the cog of mol(0)
     */
    virtual void coggather();
    /**
     * in the first frame gather all or specified molecules with respect to reference (if given) or first frame 
     * after rotational fit of reference 
     * then gather with respect to time
     * unselected molecules are always gathered with respect to cog of selected molecules
     */
    virtual void gfitgather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gathering of e.g. amyloid crystals
     */
    virtual void crsgather(); 
    /**
     * @deprecated old gathering methods which should not be used anymore
     * same but then gathering of cogs w.r.t overall cog
     */
    virtual void seqgather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * attempt for a generalized gathering method (A. Choutko / D. Geerke / A.-P. Kunz)
     */ 
    virtual void gengather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gather by using bonds rather than sequential atoms
     */
    virtual void bondgather();
    /**
     * @deprecated old gathering methods which should not be used anymore
     * gather by using a reference frame and previous frame
     */
    virtual void refgather();
  };

}
#endif


