// utils_Hbondcalc.h

// Class that contains a sequential list of specific atoms

#ifndef INCLUDED_UTILS_HBONDCALC
#define INCLUDED_UTILS_HBONDCALC

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_FSTREAM
#include <fstream>
#define INCLUDED_FSTREAM
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif
#ifndef INCLUDED_HBOND
#include "Hbond.h"
#endif
#ifndef INCLUDED_HBOND3C
#include "Hbond3c.h"
#endif
#ifndef INCLUDED_UTILS_ATOMSPECIFIER
#include "AtomSpecifier.h"
#endif
#ifndef INCLUDED_ARGS_ARGUMENTS
#include "../args/Arguments.h"
#endif



namespace gcore{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace gmath
{
  class Vec;
}

namespace args
{
  class Arguments;
}

namespace bound
{
  class Boundary;
}

namespace utils
{

  /**
   * Class Hbondcalc
   * purpose: calculate (native) intra-, intermolecular, solute-solvent and solvent-solvent 
   * hydrogen bonds over the trajectory.
   *
   * Description:
   * The Hbondcalc class calculates (native) intra-, intermolecular, solute-solvent and solvent-solvent
   * hydrogen bonds over the trajectory. It prints out two timeseries (total hydrogen bonds per frame
   * and the occurance of a specific hydrogen bond at a given time) and statistics. It makes use of
   * the AtomSpecifier if specific donor or acceptor atoms want to be specified or it can read in a
   * file containing the masses of donor and acceptor atoms and thereby make the selection.
   * @class Hbondcalc
   * @author M. Kastenholz
   * @ingroup utils
   * @sa utils::AtomSpecifier
   */
  class Hbondcalc{
    std::vector<double> d_mass_hydrogens, d_mass_acceptors;
    std::map<int, Hbond> d_hbonds;
    std::map<int, Hbond3c> d_hbonds3c;
    std::vector<double> tstime, tstime3c;
    std::vector<int> tsnum, tsnum3c;
    args::Arguments *d_args;
    gcore::System *d_sys;
    utils::AtomSpecifier d_donors, d_bound, d_acceptors;
    utils::Hbond d_hbond;
    bound::Boundary *d_pbc;
    std::ofstream timeseriesHB, timeseriesHBtot;
    std::ofstream timeseriesHB3c, timeseriesHB3ctot;
    double d_maxdist2, d_minangle, d_dt, d_time;
    double d_maxdist3c2, d_minangle3c, d_minanglesum3c, d_maxdihedral3c;
    int d_frames, d_numHB, d_numHB3c, d_nummol;
    int d_num_A_donors, d_num_A_acceptors;
    bool d_do3c;

  public: 
    // Constructor
    /**
     * Hbondcalc Constructor
     * @param sys The Hbondcalc needs to know about the system. It 
     *            does not know about any atoms yet.
     * @param args all arguments are passed into Hbondcalc. 
     */
    Hbondcalc(gcore::System &sys, args::Arguments &args);

    /**
     * Hbondcalc Deconstructor
     */
    virtual ~Hbondcalc(){timeseriesHBtot.close(); timeseriesHB.close();}
    
    
    /**
     * Method to set the maximum hydrogen bond distance.
     */    
    void setmaxdist(double i);
    /**
     * Method to set the minimum angle for a hydrogen bond.
     */
    void setminangle(double i);
    /**
     * Method to set the maximum hydrogen-acceptor distance in three centered
     * hydrogen bonds
     */
    void setmaxdist3c(double i);
    /**
     * Method to set the minimum donor-hydrogen-acceptor angle in three centered
     * hydrogen bonds
     */
    void setminangle3c(double i);
     /**
     * Method to set the minimum hydrogen-angle-sum  in three centered
     * hydrogen bonds
     */
    void setminanglesum3c(double i);
     /**
     * Method to set the maximum dihedral angle value in three centered
     * hydrogen bonds
     */
    void setmaxdihedral3c(double i);
     /**
     * Method to set the time and timestep of the trajectory frames.
     */
     void settime(double i, double j);
     /**
     * Method to read in a file specifying the donor and acceptor masses.
     */
     void readinmasses(std::string fi);
     /**
     * Method to determine the atoms specified using the AtomSpecifier.
     */
     void determineAtoms();
     /**
     * Method to determine the atoms according the read in donor/acceptor mass file.
     */
     void determineAtomsbymass();
     /**
      * Method to initialize some numbers
      */
     void init();
     /**
      * Method to do the calculation of two centered hydrogen bonds
      */
     void calc();
    /**
     * Method to do the calculation of two centered hydrogen bonds AND 
     * three centered hydrogen bonds
     */
    void calc3c();
    /**
     * Method to calculate native h-bonds
     */
    void calc_native();
    /**
     * Method to clear the statistics
     */
    void clear();
    /**
     * Method that writes the total number of hydrogen bonds per frame
     * to a file.
     */
    void writets();
   
     /**
     * Method to print the statistics for (native) intra- and intermolecular hydrogen bonds.
     */
     void printstatistics();
    /**
     * Method to print the statistics for three center hydrogen bonds.
     */
    void printstatistics3c();
    
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says Hbondcalc, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): 
	gromos::Exception("Hbondcalc", what){}
    };
  protected:
    

    /**
     * Method that opens the timeseries files.
     * 
     */
    void opents(std::string fi1, std::string fi2);
    /**
     * Method that opens the timeseries files for three centered hydrogen bonds
     */
    void opents3c(std::string fi1, std::string fi2);
    
    /**
     * Method that reads a frame from either a reference coordinat file
     * or the first frame of the first trajectory file.
     */    
    void readframe();
    /**
     * Method that calculates a single hydrogen bond
     */
    void calculate_single(int i, int j);
    /**
     * Method that calculates a single three centered hydrogen bond provided
     * that the distances were calculated already
     */
    void calculate_single3c(int i, int j, int k, double d2_1, double d2_2);
    

  }; //end class Hbondcalc

}
#endif
