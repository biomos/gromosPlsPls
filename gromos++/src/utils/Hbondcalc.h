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
    std::vector<int> d_mol, d_atom, num, tsnum;
    std::vector<double> d_mass_hydrogens, d_mass_acceptors;
    std::vector<Hbond> d_hbonds;
    std::vector<Hbond*> d_hbsolv;
    std::vector<double> dist, ang, tstime;
    args::Arguments *d_args;
    gcore::System *d_sys;
    utils::AtomSpecifier d_donors, d_donors_bound_to, d_acceptors, 
                         d_donors_solv, d_donors_bound_to_solv, d_acceptors_solv;
    utils::Hbond d_hbond;
    bound::Boundary *d_pbc;
    std::ofstream timeseriesHB, timeseriesHBtot;
    double d_maxdist, d_minangle, d_dt, d_time;
    int d_frames, d_numHB, d_nummol;
    bool d_omit_self_species;
   
   
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
     * Method to initialize the calculation for intramolecular hydrogen bonds.
     */
     void calcHintra_init();
     /**
     * Method to initialize the calculation for intermolecular hydrogen bonds.
     */
     void calcHinter_init();
     /**
     * Method to initialize the calculation for native intramolecular hydrogen bonds.
     */
     void calcHintra_native_init();
     /**
     * Method to initialize the calculation for native intermolecular hydrogen bonds.
     */
     void calcHinter_native_init();
     /**
     * Method to initialize the calculation for solute-solvent hydrogen bonds.
     */
     void calcHsolusolv_init();
     /**
     * Method to initialize the calculation for solvent-solvent hydrogen bonds.
     */
     void calcHsolvsolv_init();
     /**
     * Method to calculate intramolecular hydrogen bonds over one frame.
     */
     virtual void calcHintra();
     /**
     * Method to calculate intermolecular hydrogen bonds over one frame.
     */
     virtual void calcHinter();
     /**
     * Method to calculate solute-solvent hydrogen bonds over one frame.
     */
     virtual void calcHsolusolv();
     /**
     * Method to calculate solvent-solvent hydrogen bonds over one frame.
     */
     virtual void calcHsolvsolv();
     /**
     * Method to print the statistics for (native) intra- and intermolecular hydrogen bonds.
     */
     void printstatistics();
     /**
     * Method to print the statistics for solute-solvent hydrogen bonds.
     */
     void printstatistics_solusolv();

     typedef void (Hbondcalc::*MemPtr)();
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
     * Method that writes the total number of hydrogen bonds per frame
     * to a file.
     */
    void writets();
    /**
     * Method that reads a frame from either a reference coordinat file
     * or the first frame of the first trajectory file.
     */    
    void readframe();


  }; //end class Hbondcalc

}
#endif
