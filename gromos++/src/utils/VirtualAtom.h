// utils_VirtualAtom.h

#ifndef INCLUDED_UTILS_VIRTUALATOM
#define INCLUDED_UTILS_VIRTUALATOM

#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
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
   * -  0 : explicit/real   atom
   * -  1 : aliphatic CH1 group
   * -  2 : aromatic CH1 group
   * -  3 : non-stereospecific aliphatic CH2 group (pseudo atom)
   * -  4 : stereospecific aliphatic CH2 group
   * -  5 : single CH3 group (pseudo atom)
   * -  6 : non-stereospecific CH3 groups (isopropyl; pseudo atom)
   * -  7 : non-stereospecific CH3 groups (tert-butyl; pseudo atom)
   * - -1 : centre of geometry
   * - -2 : centre of mass
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
      normal = 0,     // explicit/real atom
      CH1 = 1,        // aliphatic CH1 group
      aromatic = 2,   // aromatic CH1 group
      CH2 = 3,        // non-stereospecific aliphatic CH2 group (pseudo atom)
      stereo_CH2 = 4, // stereospecific aliphatic CH2 group
      CH31 = 5,        // single CH31 group (pseudo atom)
      CH32 = 6,       // non-stereospecific CH3 groups (isopropyl; pseudo atom)
      CH33 = 7,       // non-stereospecific CH3 groups (tert-butyl; pseudo atom)
      COG = -1,       // centre of geometry
      COM = -2        // centre of mass
    };
    
    /**
     * Constructor
     * create a virtual atom site.
     * @param config a vector containing gromos numbers
     */
    
    VirtualAtom(gcore::System &sys, virtual_type type, std::vector<int> const &config,
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
     * Constructor
     * create from molecule and atom number based on the covalent neighbours
     * for virtual_types -1 and -2 (COG and COM). This is neede vor NOE cal-
     * culations (e.g. former type NH2 -> COG).
     */
    VirtualAtom(std::string s, gcore::System &sys, int mol, int atom, 
		virtual_type type, int subtype,
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
     * get the configuration of the virtual atom.
     * in other words the molecule and atom number
     * of the sites defining the virtual atom. Const version
     */
    const AtomSpecifier & conf() const;

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
