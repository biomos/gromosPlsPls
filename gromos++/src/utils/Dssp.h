

#ifndef INCLUDED_UTILS_DSSP
#define INCLUDED_UTILS_DSSP

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
   * Class Dssp
   * purpose: (calculate (intramolecular) backbone-backbone 
   * hydrogen bonds, and) define secondary structure of protein (DSSP)
   * over the trajectory.
   *
   * Description:
   * Dssp within gromos++ defines secondary structure for one single (protein) solute 
   * molecule, according to the DSSP rules defined by W. Kabsch and C. Sander 
   * (Biopolymers, 22, pp2577-2637 (1983)). Within these rules it may occur that one 
   * residue is defined as two different secondary-structure elements. In order to 
   * avoid duplicates in the output, the following (ad hoc) priority rules are applied 
   * here: Beta Sheet/Bridge > 5-helix > 4-helix > 3-helix > H-bonded turn > Bend. 
   * As a consequence, there may be, for instance, helices that are shorter than 
   * their minimal length. Thes rules are easily changed in Dssp::filter_SecStruct().
   * 
   * @class Dssp
   * @author U B
   * @ingroup utils
   * @sa utils::AtomSpecifier
   */
  class Dssp{
    std::vector<int> d_mol, d_atom, num, tsnum;
    std::vector<double> dist, ang, tstime;
    std::vector<int> acc_res, don_res;
    std::vector<int> helix3, helix4, helix5, Helix, helix;
    std::vector<int> bridge, extended;
    std::vector<int> Turn, turn, Bend, Beta;
    args::Arguments *d_args;
    gcore::System *d_sys;
    utils::AtomSpecifier d_H, d_N, d_O, d_C, d_CA;
    bound::Boundary *d_pbc; 
    ofstream timeseriesTurn, timeseries3Helix, timeseries4Helix, timeseries5Helix;
    ofstream timeseriesBBridge, timeseriesBStrand, timeseriesBend;
    double  d_dt, d_time;
    int d_nummol;
    bool d_omit_self_species;
   
   
  public: 
    // Constructor
    /**
     * Dssp Constructor
     * @param sys The Dssp needs to know about the system. It 
     *            does not know about any atoms yet.
     * @param args all arguments are passed into Dssp. 
     */
    Dssp(gcore::System &sys, args::Arguments &args);

    /**
     * Dssp Deconstructor
     */
    virtual ~Dssp(){
      timeseriesTurn.close();
      timeseries3Helix.close();
      timeseries4Helix.close();
      timeseries5Helix.close();
      timeseriesBBridge.close();
      timeseriesBStrand.close();
      timeseriesBend.close();
	}

     /**
     * Method to set the time and timestep of the trajectory frames.
     */
     void settime(double i, double j);
     /**
     * Method to determine the atoms specified using the AtomSpecifier.
     */
     void determineAtoms();
     /**
     * Method to initialize the calculation for intramolecular hydrogen bonds.
     */
     void calcHintra_init();
     /**
     * Method to calculate intramolecular hydrogen bonds over one frame.
     */
     virtual void calcHb_Kabsch_Sander();
     /**
      * Method to define 3-, 4-, and 5-helices over one frame.
      */
     void calc_Helices();
     /**
      * Method to define beta structures over one frame.
      */
     void calc_Betas();
     /**
      * Method to define bends over one frame.
      */
     void calc_Bends();
     /**
      * Method to filter secondary-structure output.
      */
     void filter_SecStruct();
     /**
      * Method to write secondary-structure output to files.
      */
     void writeToFiles(int n);     

    typedef void (Dssp::*MemPtr)();
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */ 
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says Dssp, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const string &what): gromos::Exception("Dssp", what){}
    };
  protected:
    /**
     * Method that opens the timeseries files.
     * 
     */
    void opents(string fi1, string fi2, string fi3, string fi4, 
		string fi5, string fi6, string fi7);
    /**
     * Method that reads a frame from either a reference coordinate file
     * or the first frame of the first trajectory file.
     */    
    void readframe();
  }; //end class Dssp
}
