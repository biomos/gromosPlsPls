/** 
 * @file:   Hbond_calc.h
 * Author: siggj
 *
 * Created on October 24, 2012, 3:31 PM
 */

#ifndef HBOND_H
#define	HBOND_H

#include "AtomSpecifier.h"


namespace gcore {
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace gmath {
  class Vec;
}

namespace args {
  class Arguments;
}

namespace bound {
  class Boundary;
}

namespace utils {
  class Pairlist;

  /**
   * Class HB_calc
   * purpose: calculate (native) intra-, intermolecular, solute-solvent and solvent-solvent 
   * hydrogen bonds over the trajectory.
   *
   * Description:
   * The HB_calc class calculates (native) intra-, intermolecular, solute-solvent and solvent-solvent
   * hydrogen bonds over the trajectory. It prints out two timeseries (total hydrogen bonds per frame
   * and the occurrence of a specific hydrogen bond at a given time) and statistics. It makes use of
   * the AtomSpecifier if specific donor or acceptor atoms want to be specified or it can read in a
   * file containing the masses of donor and acceptor atoms and thereby make the selection.
   * @class HB_calc
   * @author J. Sigg
   * @ingroup utils
   */
  class HB_calc {
  protected:
    double max_distance2;
    double min_angle;
    gcore::System *sys;
    args::Arguments *args;
    bound::Boundary *pbc;
    utils::AtomSpecifier donors, bound, acceptors;
    std::vector<double> mass_hydrogens, mass_acceptors;
    double time;
    int frames;
    int num_A_donors, num_A_acceptors;
    /**
     * Method that stores the system and all the arguments for further use.
     * @param sys
     * @param args
     */
    void setval(gcore::System& sys, args::Arguments& args);
    /**
     * Method that read in the masses from the massfile, in which the donors
     * and the acceptors are stored.
     * @param filename
     */
    void readinmasses(std::string filename);
    /**
     * Method to read in all atoms chosen by the input file.
     */
    void determineAtoms();
    /**
     * Method that cross-reference all atoms from the input file with the massfile.
     */
    void determineAtomsbymass();
    /**
     * Method to read in a frame.
     */
    void readframe();
    /**
     * Method to initialize the calculation.
     */
    void init();
    /**
     * Method that gives back true if the atom i and j are not neighbours.
     * @param i
     * @param j
     * @return 
     */
    bool neighbour(int i, int j);
    /**
     * Method that gives back true if the atom i and j have a distance smaller
     * than the maximal distance given.
     * @param dist
     * @param tmpA
     * @return 
     */
    bool distances(double &dist, gmath::Vec &tmpA);
    /**
     * Method that gives back true if the atom i and j have a angle bigger
     * than the minimal angle given.
     * @param i
     * @param angles
     * @param bound_i
     * @param tmpA
     * @return 
     */
    bool angle(int i, double &angles, gmath::Vec &bound_i, gmath::Vec &tmpA);
  public:

    /**
     * Constructor that stores the maximal distance and the minimal angle
     * @param max_distance2
     * @param min_angle
     */
    HB_calc(double max_distance2, double min_angle) :
    max_distance2(max_distance2), min_angle(min_angle) {
    }
    /**
     * Method that stores the time of the trajectories to the H-bonds.
     * @param times
     */
    void settime(double times);

  }; //end class HB_calc
}

#endif	/* HBOND_H */

