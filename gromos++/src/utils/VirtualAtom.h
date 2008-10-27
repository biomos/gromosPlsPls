// utils_VirtualAtom.h

#ifndef INCLUDED_UTILS_VIRTUALATOM
#define INCLUDED_UTILS_VIRTUALATOM

#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif
#ifndef INCLUDED_GMATH_VEC
#include "../gmath/Vec.h"
#endif
#include <iostream>
using namespace std;


namespace gcore{
  class System;
}

namespace utils{

  class VirtualAtom_i;
  class AtomSpecifier;
  
  /**
   * @class VirtualAtom
   * @author R. Buergi and M.A. Kastenholz and M. Christen
   * @ingroup utils
   *
   * Class VirtualAtom
   * This class contains information about virtual atoms
   *
   * From this information one should be able to calculate the coordinates
   * for any hydrogen that is not actually there
   *
   * @section VirtualAtom Virtual Atom 
   * The Following virtual atom types are currently supported
   * - 0 : explicit atom
   * - 1 : aliphatic CH1 group
   * - 2 : aromatic CH1 group
   * - 3 : non-stereospecific aliphatic CH2 group (pseudo atom)
   * - 4 : stereospecific aliphatic CH2 group
   * - 5 : single CH3 group (pseudo atom)
   * - 6 : non-stereospecific CH3 groups (isopropyl; pseudo atom)
   * - 7 : aromatic flipping ring (pseudo atom)
   * - 8 : non-stereospecific NH2 group (pseudo atom)
   * - 9 : non-stereospecific (CH3)3
   * - com : centre of mass
   * - cog : centre of geometry
   *
   */
  class VirtualAtom{
    VirtualAtom_i *d_this;

    // not implemented
    VirtualAtom();
    VirtualAtom &operator=(const VirtualAtom&);
  
  public:

    /**
     * @enum Virtual atom types
     */
    enum virtual_type {
      normal = 0,
      CH1 = 1,
      aromatic = 2,
      CH2 = 3,
      stereo_CH2 = 4,
      stereo_CH3 = 5,
      CH3 = 6,
      ring = 7,
      NH2 = 8,
      CH33 = 9,
      COM = 100,
      COG = 101
    };
    
    /**
     * Constructor
     * create a virtual atom site.
     */
    
    VirtualAtom(gcore::System &sys, int mol, int atom, 
		virtual_type type, std::vector<int> const &config,
		double dish = 0.1, double disc = 0.153,
		int orientation=0); 
     

    /**
     * Constructor
     * create from atom specifier
     */
    VirtualAtom(gcore::System &sys,
		AtomSpecifier const &spec,
		virtual_type type,
		double dish = 0.1, double disc = 0.153,
		int orientation = 0);
    
    /**
     * Constructor
     * create from molecule and atom number based on the covalent neighbours
     */
    VirtualAtom(gcore::System &sys, int mol, int atom, 
		virtual_type type,
                double dish = 0.1, double disc = 0.153,
                int orientation=0);

    /**
     * copy constructor
     */
    VirtualAtom(const VirtualAtom&);

    /**
     * Destructor
     */
    ~VirtualAtom();

    ////////////////////////////////////////////////////////////
    // Methods

    /**
     * calculates the virtual atom position
     */
    gmath::Vec pos()const;

    /**
     * sets carbon-hydrogen distance (0.1 by default)
     */
    void setDish(double dish);

    /**
     * sets carbon-carbon distance (0.153 by default)
     */
    void setDisc(double disc);

    /**
     * set the system.
     */
    void setSystem(gcore::System &sys);

    ////////////////////////////////////////////////////////////
    // Accessors
    /**
     * get the type of the VirtualAtom
     */
    virtual_type type()const;

    /**
     * get the configuration of the virtual atom.
     * in other words the molecule and atom number
     * of the sites defining the virtual atom.
     */
    AtomSpecifier & conf();

    /**
     * orientation for type 4 CH1
     */
    int orientation()const;
    
    /**
     * returns AtomSpecifier format like string of contents.
     */
    std::string toString()const;
    
    /**
     * @struct Exception
     * exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what) 
	: gromos::Exception("VirtualAtom", what){}
    };
    
  }; // VirtualAtom

} // util

#endif
