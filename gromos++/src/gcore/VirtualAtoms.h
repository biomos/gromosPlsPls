// gcore_VirtualAtoms.h
/**
 * Class VirtualAtoms
 */

#ifndef INCLUDED_GCORE_VIRTUALATOMS
#define INCLUDED_GCORE_VIRTUALATOMS

#include <vector>
#include <cassert>

namespace gmath{
class Vec;
}

using gmath::Vec;
namespace gcore{
class GromosForceField;
}
namespace utils{
class VirtualAtom;
class AtomSpecifier;
}

using utils::VirtualAtom;
/**
 * Class VirtualAtoms
 * Purpose: Contains virtual atoms to be part of a system
 *
 * Description:
 * The gromos++ VirtualAtoms is a specific entity. It is defined as a set of 
 * virtual atoms of which the coordinates can be computed based on the 
 * coordinates of the real atoms in the system. Because a virtual atom can be 
 * defined over multiple molecules, it is not considered part of the Molecule
 *
 * @class VirtualAtoms
 * @ingroup gcore
 * @author C. Oostenbrink
 * @sa gcore::System
 * @sa utils::VirtualAtom
 */
namespace gcore{

  class System;
  class VirtualAtoms{
    std::vector<utils::VirtualAtom> d_vas;
    std::vector<int> d_iac;
    std::vector<double> d_charge; 

  public:
    /**
     * VirtualAtoms constructor
     */
    VirtualAtoms();
    VirtualAtoms(utils::AtomSpecifier as, gcore::GromosForceField &gff);
    VirtualAtoms(const VirtualAtoms &);

    /**
     * VirtualAtoms deconstructor
     */
    ~VirtualAtoms();
    
    //Accessors
    /**
     * Accessor, returns the number of atoms in the Molecule
     */
    unsigned int numVirtualAtoms()const;
    /** 
     * Accessor, returns the i-th virtual atom
     */
    utils::VirtualAtom atom(int i);
    /** 
     * Accessor, returns the i-th virtual atom as a const
     */
    const utils::VirtualAtom atom(int i)const;
    /**
     * Accessor, returns the charge of the i-th virtual atom
     */
    double charge(int i);
    /**
     * Accessor, returns the charge of the i-th virtual atom as a const
     */
    const double charge(int i)const;
    /**
     * Accessor, returns the iac of the i-th virtual atom
     */
    int iac(int i); 
    /**
     * Accessor, returns the iac of the i-th virtual atom as a const
     */
    const int iac(int i)const; 
    
    /**
     * Method to set the IAC of th i-th virtual atom
     */
    void setIac(int i, int iac);
    /**
     * Method to set the charge of the i-th virtual atom
     */
    void setCharge(int i, double charge);
    /**
     * Method to reset the system to which the virtual atoms refer
     */
    void setSystem(gcore::System &sys);
    /** 
     * Method to add a virtual atom based on an atom specifier
     */
    void addVirtualAtom(utils::AtomSpecifier as, gcore::GromosForceField &gff, int iac=-1, double charge =0.0);
    /** 
     * Method to add a virtual atom based on a list of atoms
     */
    void addVirtualAtom(gcore::System &sys, std::vector<int> conf, int type, double dish=0.1, double disc=0.153, int iac=-1, double charge =0.0);

  }; /* class VirtualAtoms */

  inline int VirtualAtoms::iac(int i){
    assert(i < d_iac.size());
    return d_iac[i];
  }
  inline const int VirtualAtoms::iac(int i)const{
    assert(i < d_iac.size());
    return d_iac[i];
  }
  inline double VirtualAtoms::charge(int i){
    assert(i < d_charge.size());
    return d_charge[i];
  }
  inline const double VirtualAtoms::charge(int i)const{
    assert(i < d_charge.size());
    return d_charge[i];
  }
  inline void VirtualAtoms::setIac(int i, int iac){
    assert(i < d_iac.size());
    d_iac[i] = iac;
    return;
  }
  inline void VirtualAtoms::setCharge(int i, double charge){
    assert(i < d_charge.size());
    d_charge[i] = charge;
    return;
  }
  
} /* Namespace */ 
#endif

