/** 
 * @file   temperature.h
 * Implements a class, which can calculate the temperatur for a 
 * set of atoms
 */

#ifndef TEMPERATURE_H
#define	TEMPERATURE_H

namespace gcore {
    class System;
}

namespace utils {
  
  class AtomSpecifier;

  /**
   * @class Temperature
   * 
   * Calculated the temperatur for a system for a given set of atoms
   * @param as
   * @return 
   */
  class Temperature {
  public:
    //Temperature ();
    Temperature (const AtomSpecifier &as, double dof);
    double temperature(const gcore::System &sys);
  private:
    AtomSpecifier m_as;
    double dof;

  };
}
#endif	/* TEMPERATURE_H */

