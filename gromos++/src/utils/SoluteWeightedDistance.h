/** 
 * @file   SoluteWeightedDistance.h
 */

#ifndef SOLUTEAVERAGEDISTANCE_H
#define	SOLUTEAVERAGEDISTANCE_H

#include <fstream>
#include "AtomSpecifier.h"

namespace args {
  class Arguments;
}

namespace gcore {
  class System;
}

namespace bound {
  class Boundary;
}

namespace utils {
  
  /**
   * @struct SWD_Param 
   * 
   * Holds parameters (force constant and cut off) 
   * for the solute wigthed distance.
   */
  struct SWD_Param {
    SWD_Param(double f, double c) :
    forceConstant(f),
    cutoff(c)
    {}
    double forceConstant;
    double cutoff;
  };
  
  /**
   * @class SoluteWeightedDistance
   * 
   * Calculates the solute weighted distance for coarse grained and
   * fine grained solvents from the solute as well as their energies.
   */
  class SoluteWeightedDistance {
  public:
    SoluteWeightedDistance(gcore::System &sys, args::Arguments &args);
    //~SoluteWeightedDistance();
    
    /**
     * Calculate the distances
     */
    void calculate();
    /**
     * Return the title
     */
    std::string title();
    /**
     * Print distances to os.
     */
    void distances(std::ostream &os) const;
    /**
     * Print energies to os.
     */
    void energies(std::ostream &os) const;
    /**
     * Do we have the energy parameters?
     */
    bool withEnergy() const;
  protected:
    AtomSpecifier _solute;
    AtomSpecifier _fgSolvent;
    AtomSpecifier _cgSolvent;
    
    typedef std::vector<SWD_Param> Params;
    Params _params;
    bool _withEnergy;
    
    typedef std::vector<double> Distances;
    Distances _fgDistances;
    Distances _cgDistances;
    
    gcore::System &_sys;
    bound::Boundary *_pbc;
  };
  
  std::ostream &operator<<(std::ostream &os, SoluteWeightedDistance const & sad);
}

#endif	/* SOLUTEAVERAGEDISTANCE_H */

