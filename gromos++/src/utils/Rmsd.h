// utils_Rmsd.h

#ifndef INCLUDED_UTILS_RMSD
#define INCLUDED_UTILS_RMSD


namespace gcore{
  class System;
}

namespace fit{
  class Reference;
  
}

namespace utils{
   class Property;
   class PropertyContainer;
  /**
   * Class Rmsd
   * A class to calculate the rmsd compared to a reference
   *
   * The atom positional root mean square deviation of the system 
   * is calculated with respect to a reference
   *
   * @class Rmsd
   * @author R. Buergi
   * @ingroup utils
   * @sa fit::Reference
   */
  class Rmsd{
    const fit::Reference *d_ref;
    const utils::PropertyContainer *d_prop_sys;
    const utils::PropertyContainer *d_prop_ref;
    // not implemented
    Rmsd();
    Rmsd(const Rmsd &);
    Rmsd &operator=(const Rmsd&);
  public:
    /**
     * Rmsd constructor
     */
    Rmsd(const fit::Reference *);    
    /**
     * Rmsd deconstructor
     */
    ~Rmsd(){}
    /** 
     * Method that actually calculate the positional rmsd between the reference
     * and the system
     * @param & System to calculate the rmsd for
     * @return the rmsd value
     */
    double rmsd(const gcore::System &);
    /** 
     * Method that actually calculate a property rmsd of a PropertyContainer
     * between the reference and the system
     * @param & System to calculate the rmsd for
     * @return the rmsd value
     */
    double rmsdproperty(const gcore::System &);

    /** 
     * Method to add a PropertyContainer for the rmsd calculation
     * @param * PropertyContainer of the reference system
     * @param * PropertyContainer of the system
     */
    void addproperty(const utils::PropertyContainer *, const utils::PropertyContainer *);

    typedef double (Rmsd::*MemPtr)(const gcore::System &);

  };
}

#endif
