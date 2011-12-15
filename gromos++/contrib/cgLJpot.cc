/**
 * @file cgLJpot
 * used to develop a coarse-grained potential from an all-atom simulation
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor cgLJpot
 * @section cgLJpot calculates coarse-grained LJ potentials form fine grained simulations
 * @author @ref ae
 * @date 07.11.2011
 *
 * Program cgLJpot calculates the Lennard-Jones potential energy function for a coarse-grained system based on a fine-grained (all-atom) simulation trajectory
 * (V_fg2cg). The resulting potential energy is not a 12/6-Lennard-Jones potential. However, the program also calculates the coarse-grained 12/6-Lennard-Jones
 * potential energy function as an approximation to V_fg2cg, both having a minimum at the same energy and r value, i.e. the sigma- and epsilon-value for the
 * approximation to V_fg2cg is read from the minimum of V-fg2cg. Therefor, in the case the beads do not include atoms of different molecules, the calculated
 * 12/6-Lennard-Jones potential is expected to reproduce the density and heat of vaporization of the fine-grained simulation in a coarse-grained simulation.
 * However, practice showed that this is normally not the case and parameterization is still necessary, but the Lennar-Jones parameters are reasonable to start
 * the parameterization.
 * 
 * NOTE: the current program version does only work for solute
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@method</td><td>&lt;method to goarse grain: atomic or molecular&gt; </td></tr>
 * <tr><td> [\@dist</td><td>&lt;min max ngrid&gt;] </td></tr>
 * <tr><td> \@beads</td><td>&lt;number of atoms per bead (atomic)&gt; or </td></tr>
 * <tr><td>        </td><td>&lt;sequence of bead size within one molecule (molecular)&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt; </td></tr>
 * <tr><td> \@trc</td><td>&lt;simulation trajectory or coordinate file&gt;</td></tr>
 * </table>
 *
 * 
 * Example:
 * @verbatim
   rdf
     @topo   ex.top
 * @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <ctime>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomPair.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Distribution.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/System.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Physics.h"

namespace cgLJpot {

  /**
   * A class to store a pair of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJ {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;

  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     */
    IJ(int i = -1, int j = -1);
    /**
     * copy constructor 
     * @param ij the IJ class to be copied
     */
    IJ(const IJ &ij);
    void setValues(int i, int j);
    int i() const;
    int j() const;
  };
  
  /**
   * A class to store a quartet of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJK {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;
    /**
     * the third integer number
     */
    int K;


  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     * @param k third integer number; standard: -1
     */
    IJK(int i = -1, int j = -1, int k = -1);
    /**
     * copy constructor 
     * @param ijkl the IJKL class to be copied
     */
    IJK(const IJK &ijk);
    void setValues(int i, int j, int k);
    int i() const;
    int j() const;
    int k() const;
  };
  
  /**
   * A class to store a quartet of integers, e.g. the two IAC numbers of two particles interacting via a LJ potential
   */
  class IJKL {
  private:
    /**
     * the first integer number 
     */
    int I;
    /**
     * the second integer number
     */
    int J;
    /**
     * the third integer number
     */
    int K;
    /**
     * the fourth integer number
     */
    int L;

  public:
    /**
     * constructor
     * @param i first integer number; standard: -1
     * @param j second integer number; standard: -1
     * @param k third integer number; standard: -1
     * @param l fourth integer number; standard: -1
     */
    IJKL(int i = -1, int j = -1, int k = -1, int l = -1);
    /**
     * copy constructor 
     * @param ijkl the IJKL class to be copied
     */
    IJKL(const IJKL &ijkl);
    void setValues(int i, int j, int k, int l);
    int i() const;
    int j() const;
    int k() const;
    int l() const;
  };

  /**
   * just a < operator to make the class more complete... 
   * this is needed to iterate over a map using @ref IJs as the key value
   * @param ij1 first @ref IJ to be compared with...
   * @param ij2 ... the second @ref IJ
   */
  bool operator<(const IJ &ij1, const IJ &ij2);

  /**
   * stores the Lennard-Jones-potential of multiple pairs of particles on a grid */
  class LJpot {
  private:
    /**
     * the summed up LJ potential energies at each grid
     */
    vector<double> lj;
    /**
     * number of summed up potential energies at each grid
     */
    vector<int> count;
    /**
     * the number of grids, used to size the vectors @ref lj and @ref count
     */
    double dgrid;
    /**
     * the minimum distance between two particles, pairs of particles with a smaller distance are not considered
     */
    double min;
    /**
     * the maximum distance between two particles, pairs of particles with a longer distance are not considered
     */
    double max;
    /**
     * the Lennard-Jones C12 parameter
     */
    double C12;
    /**
     * the Lennard-Jones C6 parameter
     */
    double C6;

  public:
    /**
     * constructor
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     * @param c12 Lennard-Jones C12 parameter
     * @param c6 Lennard-Jones C6 parameter
     */
    LJpot(double min_, double max_, int grid, double c12, double c6);
    /**
     * constructor
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     */
    LJpot(double min_ = 0.0, double max_ = 2.0, int grid = 200);
    /**
     * copy constructor
     * @param ljp a Lennard-Jones potential class of type @ref LJ
     */
    LJpot(const LJpot &ljp);
    /**
     * initializer
     * @param min_ minimum interatomic distance
     * @param max_ maximum interatomic distance
     * @param grid number of grid points for interatomic distance
     */
    void init(double min_ = 0.1, double max_ = 2.0, int grid = 200);
    /**
     * adds the Lennard-Jones energy of two particles to the stored potential energy
     * @param pos inter-particle distance
     * @param val Lennard-Jones energy
     */
    void add(double pos, double val);
    /**
     * unifies two @ref LJpot classes into one
     * @param ljp the @LJpot class which will be unified with the current one
     */
    void unify(const LJpot &ljp);
    /**
     * accessor: returns the inter-atomic distance r corresponding to grid point i
     * @param i grid point number
     * @return interatomic distance r
     */
    double r(unsigned int i);
    /**
     * accessor: returns the inter-atomic distance r corresponding to grid point i
     * @param i grid point number
     * @return interatomic distance r
     */
    double r(unsigned int i) const;
    /**
     * accessor: returns the (averaged) Lennard-Jones potential corresponding to grid point i, i.e. the summed up Lennard-Jones potential energy divided by
     * the number if summed up terms at grid point i
     * @param i grid point number
     * @return averaged Lennard-Jones potential energy
     */
    double pot(unsigned int i);
    /**
     * accessor: returns the minimum inter-particle distance to be considered
     * @return minimum inter-particle distance to be considered
     */
    double get_min();
    /**
     * accessor: returns the maximum inter-particle distance to be considered
     * @return maximum inter-particle distance to be considered
     */
    double get_max();
    /**
     * accessor: returns the number of grid points
     * @return number of grid points
     */
    int get_grid();
    /**
     * accessor: searches the inter-particle distance r_max for which the Lennard-Jones potential is maximal and returns the inter-particle distance r > r_max
     * for which the Lennard-Jones potential is minimal
     * @return inter-particle distance r for which the Lennard-Jones energy is minimal
     */
    double rmin();
    /** 
     * accessor: searches the inter-particle distance r_max for which the Lennard-Jones potential is maximal and returns the minimum Lennard-Jones potential
     * energu V_LJ(r) with r > r_max
     * @return inter-particle distance r for which the Lennard-Jones energy is minimal
     */
    double potmin();
  };

  /**
   * stores multiple atoms to act as a bead as well as the different interaction Lennard-Jones potential energies to other beads
   */
  class bead {
  private:
    /**
     * a reference to the gromos force field
     */
    gcore::GromosForceField *gff;
    /**
     * the atoms which ar building this bead
     */
    utils::AtomSpecifier atoms;
    /**
     * the center of geometry or center of mass of all the atoms of the bead -> center of the bead
     */
    gmath::Vec centre;
    /**
     * in case all atoms of this bead belong to the same molecule this number indicates the molecule number
     */
    int memberOfMol;
    /**
     * a sequential number to tag the beads, either within a molecule or within the whole system
     */
    int beadnum;
    /**
     * stores the total Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> totLJ;
    /**
     * stores the total inter-molecular Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J,
     * respectively
     */
    map<IJ, LJpot> totinterLJ;
    /**
     * stores the total intra-molecular Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> totintraLJ;
    /**
     * stores the inter-molecular 12-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra12LJ;
    /**
     * stores the inter-molecular 13-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra13LJ;
    /**
     * stores the inter-molecular 14-neighboring Lennard-Jones potentials for the different IJ combinations; IJ indicates e.g. two IAC number or beads with size I and J, respectively
     */
    map<IJ, LJpot> intra14LJ;

  public:
    /**
     * constructor
     * @param sys the system
     * @param groff the GROMOS force field
     * @param mom the molecule the bead is a member of
     * @param bnum the sequential bead number
     * @param ij the possibly different IJ combinations (to other beads)
     * @param min the minimum distance to another bead for which the Lennard-Jones potential energy is calculated and stored
     * @param max the maximum distance to another bead for which the Lennard-Jones potential energy is calculated and stored 
     * @param grid the number of grid points of the Lennard-Jones potentials
     */
    bead(gcore::System &sys, gcore::GromosForceField &groff, int mom, int bnum, set<IJ> &ij, double min = 0.0, double max = 2.0, double grid = 200);
    /**
     * copy constructor
     * @param b another bead to be copied
     */
    bead(bead const &b);
    /**
     * destructor
     */
    ~bead() {
    };
    /**
     * adds an atom to the current bead
     * @param m molecule number of the atom to be added
     * @param a atom number of the atom to be added
     * @return total number of atoms within the bead
     */
    int addAtom(int m, int a);
    /**
     * returns the number of atoms within that bead
     */
    int size();
    /**
     * calculates, defines and returns the center of geometry of the bead
     * @param pbc periodic boundary, see @ref bound::Boundary
     * @param sys the system
     */
    gmath::Vec cog(bound::Boundary *pbc, gcore::System &sys);
    /**
     * calculates, defines and returns the center of mass
     * @param pbc periodic boundary, see @ref bound::Boundary
     * @param sys the system
     */
    gmath::Vec com(bound::Boundary *pbc, gcore::System &sys);
    /**
     * accessor to the (center) position of the bead (com or cog)
     * */
    gmath::Vec pos();
    /**
     * add a value to the total LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtot(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtotinter(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total intramolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJtotintra(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 12-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra12(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 13-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra13(const IJ &ij, double r, const double &lj);
    /**
     * add a value to the total 14-intermolecular LJ potential energy
     * @param ij the combination @ref IJ
     * @param r the inter-particle distance
     * @param lj the Lennard-Jones energy
     */
    void addLJintra14(const IJ &ij, double r, const double &lj);
    /** calculate the LJ interaction to another bead
     * @param b a @ref bead
     * @param pbc
     * @param sys
     * @return 
     */
    double calcLJ(bead &b, bound::Boundary *pbc, gcore::System &sys);
    /**
     * accessor: returns the total LJpot energy
     * @return total LJpot energy
     */
    map<IJ, LJpot> get_totLJ();
    /**
     * accessor: returns the total intermolecular LJpot energy
     * @return total intermolecular LJpot energy
     */
    map<IJ, LJpot> get_totinterLJ();
    /**
     * accessor: returns the total intramolecular LJpot energy
     * @return total intramolecular LJpot energy
     */
    map<IJ, LJpot> get_totintraLJ();
    /**
     * accessor: returns the 12-intermolecular LJpot energy
     * @return 12-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra12LJ();
    /**
     * accessor: returns the 13-intermolecular LJpot energy
     * @return 13-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra13LJ();
    /**
     * accessor: returns the 14-intermolecular LJpot energy
     * @return 14-intermolecular LJpot energy
     */
    map<IJ, LJpot> get_intra14LJ();
    /**
     * accessor: returns the molecule number the bead/atoms of the bead belong to
     */
    int mol();
    /**
     * accessor: returns the bead number within the molecule (sequential bead number)
     */
    int beadNum();
  };

  /**
   * prints the different potentials in columns where the first column is the inter-particle distance r
   * @param os the output stream
   * @param totLJpot the total Lennard-Jones potential
   * @param totLJinter the total inter-particle Lennard-Jones potential
   * @param totLJintra the total intra-particle Lennard-Jones potential
   * @param totLJ12 the total 12-Lennard-Jones potential
   * @param totLJ13 the total 13-Lennard-Jones potential
   * @param totLJ14 the total 14-Lennard-Jones potential
   */
  void printPot(ostream &os, std::map<cgLJpot::IJ, LJpot> &totLJpot,
          std::map<IJ, LJpot> &totLJinter,
          std::map<IJ, LJpot> &totLJintra,
          std::map<IJ, LJpot> &totLJ12,
          std::map<IJ, LJpot> &totLJ13,
          std::map<IJ, LJpot> &totLJ14);

  /**
   * prints some header information to remember what the program was analyzing
   */
  void printTitleBlock(std::vector<int> &beadsizes, utils::AtomSpecifier &allatoms);
  /**
   * calculates the sigma and epsilon values of a Lennard-Jones potential @ref LJpot
   */
  void calcEpsSigma(std::map<IJ, double> &epsilons, std::map<IJ, double> &sigmas, std::map<IJ, LJpot> &LJ);
  /**
   * calculates the Lennard-Jones C12 and C6 parameters based on the corresponding epsilon- and sigma values
   */
  void calcC12C6(std::map<IJ, double> &C12s, std::map<IJ, double> &C6s, std::map<IJ, double> &epsilons, std::map<IJ, double> &sigmas);
  /**
   * prints the various Lennard-Jones parameters
   */
  void printLennardJonesParamters(map<IJ, double> &epsilons_tot,
          map<IJ, double> &epsilons_totinter,
          map<IJ, double> &epsilons_totintra,
          map<IJ, double> &epsilons_intra12,
          map<IJ, double> &epsilons_intra13,
          map<IJ, double> &epsilons_intra14,
          map<IJ, double> &sigmas_tot,
          map<IJ, double> &sigmas_totinter,
          map<IJ, double> &sigmas_totintra,
          map<IJ, double> &sigmas_intra12,
          map<IJ, double> &sigmas_intra13,
          map<IJ, double> &sigmas_intra14,
          map<IJ, double> &C12_tot,
          map<IJ, double> &C12_totinter,
          map<IJ, double> &C12_totintra,
          map<IJ, double> &C12_intra12,
          map<IJ, double> &C12_intra13,
          map<IJ, double> &C12_intra14,
          map<IJ, double> &C6_tot,
          map<IJ, double> &C6_totinter,
          map<IJ, double> &C6_totintra,
          map<IJ, double> &C6_intra12,
          map<IJ, double> &C6_intra13,
          map<IJ, double> &C6_intra14);
  
}

using namespace args;
using namespace bound;
using namespace cgLJpot;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;
using namespace utils;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "method" << "beads" << "pbc" << "trc" << "output" << "dist";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@method        <method to goarse grain: atomic or molecular>\n";
  usage += "\t[@dist         <min max ngrid>]\n";
  usage += "\t@beads         <number of atoms per bead (atomic)> or\n";
  usage += "\t               <sequence of bead size within one molecule (molecular)>\n";
  usage += "\t[@pbc          <boundary type (read from GENBOX block if not specified)> [<gather method>]]\n";
  usage += "\t[@output       fg2cg   <file name>\n";
  usage += "\t               cg      <file name>\n";
  usage += "\t               bbdist  <file name>]\n";
  usage += "\t@trc           <simulation trajectory or coordinate file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read the method:
    // - atomic: the program does not care about the moleculs and where they start and end
    // - the program cares about the molecules (no beads over more than one molecule
    args.check("method", 1);
    string method = args["method"];
    if (method != "atomic" && method != "molecular") {
      stringstream msg;
      msg << "method \"" << method << "\" (@method) not implemented, chose \"atomic\" or \"molecular\"";
      throw gromos::Exception(argv[0], msg.str());
    }
    
    string fname_LJpot_FG2CG = "LJpot_FG2CG.dat";
    string fname_LJpot_CG = "LJpot_CG.dat";
    string fname_beadbead_dist = "bead-bead_dist.dat";
    if (args.count("output") > 0) {
      int fg2cg = 0;
      int cg=0;
      int bb=0;
      if(args.check("output", 6) == 0) {
        Arguments::const_iterator it = args.lower_bound("output");
        for(; it != args.upper_bound("output"); ++it) {
          if(it->second == "fg2cg") {
            fg2cg = 1;
            it++;
            fname_LJpot_FG2CG = it->second;
          } else if (it->second == "cg") {
            cg = 1;
            it++;
            fname_LJpot_CG = it->second;
          } else if(it->second == "bbdist") {
            bb = 1;
            it++;
            fname_beadbead_dist = it->second;
          } else {
            stringstream ss;
            ss << it->second << " unknown, check arguments for @output";
            throw gromos::Exception(argv[0], ss.str());
          }
        }
        if(fg2cg + cg + bb != 3) {
          throw gromos::Exception(argv[0], "check arguments for @output");
        }
      }
    }

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff = it.forceField();
    
    // read the bead size
    args.check("beads", 1);
    vector<int> beadsizes;
    {
      Arguments::const_iterator start = args.lower_bound("beads");
      Arguments::const_iterator stop = args.upper_bound("beads");
      Arguments::const_iterator it;
      for (it = start; it != stop; ++it) {
        stringstream ss;
        ss << it->second;
        int b;
        ss >> b;
        if (ss.fail() || ss.bad() || !ss.eof()) {
          stringstream msg;
          msg << "cannot use " << it->second << " as a bead size";
          throw gromos::Exception(argv[0], msg.str());
        }
        beadsizes.push_back(b);
      }
    }
    if (method == "atomic" && beadsizes.size() != 1) {
      throw gromos::Exception(argv[0], "method \"atomic\" (@atom) does not allow for more than one bead size");
    }

    // read the distribution parameters
    double distmin = 0.0;
    double distmax = 2.0;
    int distgrid = 200;
    if (args.count("dist") >= 0) {
      if (args.count("dist") == 3) {
        stringstream ss;
        Arguments::const_iterator it = args.lower_bound("dist");
        ss << it->second << endl;
        ++it;
        ss << it->second << endl;
        ++it;
        ss << it->second;
        ss >> distmin >> distmax >> distgrid;
        if (ss.fail() || ss.bad() || !ss.eof()) {
          stringstream msg;
          msg << "cannot convert the arguments of @dist to min, max and number of grid points";
          throw gromos::Exception(argv[0], msg.str());
        }
      } else {
        stringstream msg;
        msg << "cannot convert the arguments of @dist to min, max and number of grid points";
        throw gromos::Exception(argv[0], msg.str());
      }
      if (distmin >= distmax) {
        throw gromos::Exception(argv[0], "distmin >= distmax but should be distmin < distmax");
      }
    }

    // do some checks depending on the coarse-grain method
    int numMol = sys.numMolecules();
    int numAtMol = sys.mol(0).numAtoms();
    int checkMol = -1; // -1 if all molecules have the same number of atoms, otherwise
    // it gets is assigned to the molecule number of the first molecule
    // having a different atom number
    AtomSpecifier allAtoms(sys);
    for (int m = 0; m < numMol; ++m) {
      int numAt = sys.mol(m).numAtoms();
      if (numAt != numAtMol && checkMol == -1) {
        checkMol = m;
      }
      for (int a = 0; a < numAt; ++a) {
        allAtoms.addAtom(m, a);
      }
    }
    // checks for method = atomic:
    // is the total number of atoms "dividable" (without rest) by the bead size?
    int numBeads = 0;
    if (method == "atomic") {
      if (allAtoms.size() % beadsizes[0] != 0) {
        stringstream msg;
        msg << "the total number of atoms (" << allAtoms.size() << ") cannot be divided"
                " by the bead size (" << beadsizes[0] << "): " << allAtoms.size() % beadsizes[0]
                << " atoms left";
        throw gromos::Exception(argv[0], msg.str());
      }
      numBeads = allAtoms.size() / beadsizes[0];
    }
    // checks for method = molecular
    if (method == "molecular") {
      if (checkMol != -1) {
        stringstream msg;
        msg << "molecule " << checkMol + 1 << " has a different number of atoms"
                " than the previous molecules: method \"molecular\" (@method) "
                "does not work therefore";
        throw gromos::Exception(argv[0], msg.str());
      }
      int numAtMolBeads = 0;
      for (unsigned int b = 0; b < beadsizes.size(); ++b) {
        numAtMolBeads += beadsizes[b];
      }
      if (numAtMol != numAtMolBeads) {
        stringstream msg;
        msg << "the number of atoms per molecule is " << numAtMol << ", but the"
                " different bead sizes sum up to " << numAtMolBeads;
        throw gromos::Exception(argv[0], msg.str());
      }
      numBeads = allAtoms.size() / numAtMolBeads * beadsizes.size();
    }

    // get all possible combinations of beads (with respect to its size
    set<IJ> IJs;
    set<IJK> IJKs;
    set<IJKL> IJKLs;
    for (unsigned int bs1 = 0; bs1 < beadsizes.size(); ++bs1) {
      for (unsigned int bs2 = bs1; bs2 < beadsizes.size(); ++bs2) {
        IJ ij(beadsizes[bs1], beadsizes[bs2]);
        if (IJs.find(ij) == IJs.end()) {
          IJs.insert(ij);
        }
        for (unsigned int bs3 = 0; bs3 < beadsizes.size(); ++bs3) {
          IJK ijk(beadsizes[bs1], beadsizes[bs3], beadsizes[bs2]);
          if (IJKs.find(ijk) == IJKs.end()) {
            IJKs.insert(ijk);
            for (unsigned int bs4 = 0; bs4 < beadsizes.size(); ++bs4) {
              IJKL ijkl(beadsizes[bs3], beadsizes[bs1], beadsizes[bs2], beadsizes[bs4]);
              if (IJKLs.find(ijkl) == IJKLs.end()) {
                IJKLs.insert(ijkl);
              }
            }
          }
        }
      }
    }

    // construct the beads
    vector<bead> beads;
    {
      int a = 0; // the current position in allAtoms
      while (a < allAtoms.size()) {
        int bnum = 0;
        for (int bs = 0; bs < (int) beadsizes.size(); ++bs, ++bnum) {
          bead B(sys, gff, allAtoms.mol(a), bnum, IJs, distmin, distmax, distgrid);
          for (int i = 0; i < beadsizes[bs]; ++i) {
            int mol = allAtoms.mol(a);
            int at = allAtoms.atom(a);
            B.addAtom(mol, at);
            ++a;
          }
          beads.push_back(B);
        }
      }
    }
    
    printTitleBlock(beadsizes, allAtoms);

    // a map of distributions (for each IJ one) to remember the intramolecular
    // neighboring bead-bead distance (min and max automatically set based on the
    // first configuration of the trajectories
    map<IJ, Distribution> beadbeadDist;
    map<IJK, Distribution> angleDist;
    map<IJKL, Distribution> dihedralDist;
    double bondlength_min = sys.box().K().abs() + sys.box().L().abs() + sys.box().M().abs();
    double bondlength_max = 0.0;
    {
      Arguments::const_iterator trcfirs = args.lower_bound("trc");
      InG96 ic;
      ic.open(trcfirs->second.c_str());
      ic.select("SOLUTE");
      ic >> sys;
      ic.close();
      bound::Boundary *pbc;
      if (args.count("pbc") > 0) { // read from arguments
        pbc = args::BoundaryParser::boundary(sys, args);
      } else { // read from GENBOX block
        if (args.count("pbc") == 0) {
          cerr << "WARNING: @pbc given with no argument(s), reading boundary type "
                  "from GENBOX block of the trajectory/coordinate file\n";
        }
        pbc = args::BoundaryParser::boundary(sys);
      }
      // loop over the beads and get the min/max distance
      for (int b1 = 0; b1 < (int) beads.size() - 1; ++b1) {
        beads[b1].com(pbc, sys);
        for (int b2 = b1 + 1; b2 < (int) beads.size(); ++b2) {
          if (beads[b1].mol() == beads[b2].mol() &&
                  abs(beads[b1].beadNum() - beads[b2].beadNum()) == 1) {
            beads[b2].com(pbc, sys);
            double r = (beads[b1].pos() - pbc->nearestImage(beads[b1].pos(), beads[b2].pos(), sys.box())).abs();
            if (r < bondlength_min) {
              bondlength_min = r;
            }
            if (r > bondlength_max) {
              bondlength_max = r;
            }
          } else {
            continue;
          }
        }
      }
      double min = bondlength_min - (1.2 * bondlength_max - bondlength_max);
      bondlength_max *= 1.2;
      bondlength_min = min > 0 ? min * 0.8 : 0.0;
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        beadbeadDist.insert(pair<IJ, Distribution > (*it, Distribution(bondlength_min, bondlength_max, 200)));
      }
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
        angleDist.insert(pair<IJK, Distribution> (*it, Distribution(0, 180, 180)));
      }
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        dihedralDist.insert(pair<IJKL, Distribution> (*it, Distribution(0, 180, 180)));
      }
    }

    // loop over the different trajectory files
    if (args.count("trc") < 1) {
      throw gromos::Exception(argv[0], "no coordinate or trajectory file specified (@trc)");
    }
    Arguments::const_iterator trcfirs = args.lower_bound("trc");
    Arguments::const_iterator trclast = args.upper_bound("trc");
    for (args::Arguments::const_iterator trc = trcfirs;
            trc != trclast; ++trc) {

      // the input coordinates
      InG96 ic;

      // the boundary
      bound::Boundary *pbc;

      // read boundary type, either from @pbc or GENBOX block
      if (args.count("pbc") > 0) { // read from arguments
        pbc = args::BoundaryParser::boundary(sys, args);
      } else { // read from GENBOX block
        if (args.count("pbc") == 0) {
          cerr << "WARNING: @pbc given with no argument(s), reading boundary type "
                  "from GENBOX block of the trajectory/coordinate file\n";
        }
        ic.open(trc->second.c_str());
        ic.select("SOLUTE");
        ic >> sys;
        pbc = args::BoundaryParser::boundary(sys);
        ic.close();
      }

      // loop over the configurations of the trajectory file
      ic.open(trc->second.c_str());
      ic.select("SOLUTE");
      while (!ic.eof()) {
        ic >> sys;

        // check if the box length is as least as big as twice the cut-off radius
        double L = sys.box().K().abs2();
        if (sys.box().L().abs2() < L) {
          L = sys.box().L().abs2();
        }
        if (sys.box().M().abs2() < L) {
          L = sys.box().M().abs2();
        }
        if (distmax > (sqrt(L) / 2.0)) {
          throw gromos::Exception(argv[0], "maximal @dist value bigger than "
                  "1/2 of the box length");
        }

        // calculate all centres of the beads (com)
        for (int b = 0; b < (int) beads.size(); ++b) {
          beads[b].com(pbc, sys);
        }

        // double loop over the beads
        int b2;
        double lj;
        double r2;
        double r;
        IJ ij;
#ifdef OMP
#pragma omp parallel for private(b2, ij, lj, r2, r) schedule(dynamic)
#endif
        for (int b1 = 0; b1 < (int) beads.size() - 1; ++b1) {
          for (b2 = b1 + 1; b2 < (int) beads.size(); ++b2) {
            r2 = (beads[b1].pos() - pbc->nearestImage(beads[b1].pos(),
                    beads[b2].pos(), sys.box())).abs2();
            // if the two beads are within the range of the distribution range,
            // calculate the LJ potential energy
            if (distmin * distmin <= r2 && distmax * distmax > r2) {
              if (method == "molecular") {
                ij.setValues(beads[b1].size(), beads[b2].size());
                // calculate the bead-bead LJ potential energy
                lj = beads[b1].calcLJ(beads[b2], pbc, sys);
                r = std::sqrt(r2);
                // and add the result to the two beads (corresponding potentials)
#ifdef OMP
#pragma omp critical
#endif               
                {
                  beads[b1].addLJtot(ij, r, lj);
                  beads[b2].addLJtot(ij, r, lj);
                  if (beads[b1].mol() != beads[b2].mol()) { // intermolecular LJ potential
                    beads[b1].addLJtotinter(ij, r, lj);
                    beads[b2].addLJtotinter(ij, r, lj);
                  } else { // intramolecular LJ potential energy
                    beads[b1].addLJtotintra(ij, r, lj);
                    beads[b2].addLJtotintra(ij, r, lj);
                    // 12-LJ pot (neighboring beads)
                    if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 1) {
                      beads[b1].addLJintra12(ij, r, lj);
                      beads[b2].addLJintra12(ij, r, lj);
                      // add the distance to the distribution
                      beadbeadDist[ij].add(r);
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 2) {
                      beads[b1].addLJintra13(ij, r, lj);
                      beads[b2].addLJintra13(ij, r, lj);
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 3) {
                      beads[b1].addLJintra14(ij, r, lj);
                      beads[b2].addLJintra14(ij, r, lj);
                    }
                  }
                }
              } else {
                stringstream msg;
                msg << "method " << method << " not implemented";
                throw gromos::Exception(argv[0], msg.str());
              }
            }
          }
        }
        
        // bond, angle and dihedral distribution (actually, bonds are done above
        // and are not done here any more...)
        if (method == "molecular") {
          // do the angles
          for (int b = 0; b < (int) beads.size() - 2; ++b) {
            if((beads[b].beadNum() < beads[b + 1].beadNum()) &&
                    (beads[b].beadNum() < beads[b + 2].beadNum())) {
              Vec p2 = beads[b + 1].pos();
              Vec p1 = pbc->nearestImage(p2, beads[b].pos(), sys.box());
              Vec p3 = pbc->nearestImage(p2, beads[b + 2].pos(), sys.box());
              Vec v1 = p1 - p2;
              Vec v2 = p3 - p2;
              double a = acos(v1.dot(v2) / (v1.abs() * v2.abs())) / physConst.get_pi() * 180;
              IJK ijk(beads[b].size(), beads[b+1].size(), beads[b+2].size());
              angleDist[ijk].add(a);
            }
            // do the dihedrals
            if(b < int(beads.size() - 3) && (beads[b].beadNum() < beads[b + 1].beadNum()) && 
                    (beads[b].beadNum() < beads[b + 2].beadNum()) && 
                    (beads[b].beadNum() < beads[b + 3].beadNum())) {
              Vec tmpA = beads[b].pos() - pbc->nearestImage(beads[b].pos(), beads[b+1].pos(), sys.box());
              Vec tmpB = beads[b+3].pos() - pbc->nearestImage(beads[b+3].pos(), beads[b+2].pos(), sys.box());
              Vec tmpC = beads[b+2].pos() - pbc->nearestImage(beads[b+2].pos(), beads[b+1].pos(), sys.box());
              Vec p1 = tmpA.cross(tmpC);
              Vec p2 = tmpB.cross(tmpC);
              double cosphi = ((p1.dot(p2)) / (p1.abs() * p2.abs()));
              if (cosphi > 1.0) cosphi = 1.0;
              if (cosphi <-1.0) cosphi = -1.0;
              double d = acos(cosphi) * 180 / physConst.get_pi();
              IJKL ijkl(beads[b].size(), beads[b+1].size(), beads[b+2].size(), beads[b+3].size());
              dihedralDist[ijkl].add(d);
            }
          }
        } else {
          stringstream msg;
          msg << "method " << method << " not implemented";
          throw gromos::Exception(argv[0], msg.str());
        }
        
      } // end of loop over the configurations of the trajectory file
      ic.close();

    } // end of loop over the different trajectory files

    // unify the calculated potential of all beads
    map<IJ, LJpot> totLJ;
    map<IJ, LJpot> totinterLJ;
    map<IJ, LJpot> totintraLJ;
    map<IJ, LJpot> intra12LJ;
    map<IJ, LJpot> intra13LJ;
    map<IJ, LJpot> intra14LJ;
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      totLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totinterLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      totintraLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra12LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra13LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
      intra14LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid)));
    }
    for (unsigned int b = 0; b < beads.size(); ++b) {
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        totLJ[*it].unify(beads[b].get_totLJ()[*it]);
        totinterLJ[*it].unify(beads[b].get_totinterLJ()[*it]);
        totintraLJ[*it].unify(beads[b].get_totintraLJ()[*it]);
        intra12LJ[*it].unify(beads[b].get_intra12LJ()[*it]);
        intra13LJ[*it].unify(beads[b].get_intra13LJ()[*it]);
        intra14LJ[*it].unify(beads[b].get_intra14LJ()[*it]);
      }
    }

    // print the different potentials
    {
      ofstream fout(fname_LJpot_FG2CG.c_str());
      printPot(fout, totLJ, totinterLJ, totintraLJ, intra12LJ, intra13LJ, intra14LJ);
      fout.close();
    }

    
    // calculate and print the resulting LJ pot for the CG system
    map<IJ, double> sigmas_tot;
    map<IJ, double> sigmas_totinter;
    map<IJ, double> sigmas_totintra;
    map<IJ, double> sigmas_intra12;
    map<IJ, double> sigmas_intra13;
    map<IJ, double> sigmas_intra14;
    map<IJ, double> epsilons_tot;
    map<IJ, double> epsilons_totinter;
    map<IJ, double> epsilons_totintra;
    map<IJ, double> epsilons_intra12;
    map<IJ, double> epsilons_intra13;
    map<IJ, double> epsilons_intra14;
    map<IJ, double> C12_tot;
    map<IJ, double> C12_totinter;
    map<IJ, double> C12_totintra;
    map<IJ, double> C12_intra12;
    map<IJ, double> C12_intra13;
    map<IJ, double> C12_intra14;
    map<IJ, double> C6_tot;
    map<IJ, double> C6_totinter;
    map<IJ, double> C6_totintra;
    map<IJ, double> C6_intra12;
    map<IJ, double> C6_intra13;
    map<IJ, double> C6_intra14;
    calcEpsSigma(epsilons_tot, sigmas_tot, totLJ);
    calcEpsSigma(epsilons_totinter, sigmas_totinter, totinterLJ);
    calcEpsSigma(epsilons_totintra, sigmas_totintra, totintraLJ);
    calcEpsSigma(epsilons_intra12, sigmas_intra12, intra12LJ);
    calcEpsSigma(epsilons_intra13, sigmas_intra13, intra13LJ);
    calcEpsSigma(epsilons_intra14, sigmas_intra14, intra14LJ);
    calcC12C6(C12_tot, C6_tot, epsilons_tot, sigmas_tot);
    calcC12C6(C12_totinter, C6_totinter, epsilons_totinter, sigmas_totinter);
    calcC12C6(C12_totintra, C6_totintra, epsilons_totintra, sigmas_totintra);
    calcC12C6(C12_intra12, C6_intra12, epsilons_intra12, sigmas_intra12);
    calcC12C6(C12_intra13, C6_intra13, epsilons_intra13, sigmas_intra13);
    calcC12C6(C12_intra14, C6_intra14, epsilons_intra14, sigmas_intra14);
    printLennardJonesParamters(epsilons_tot, epsilons_totinter, epsilons_totintra, epsilons_intra12, epsilons_intra13, epsilons_intra14,
            sigmas_tot, sigmas_totinter, sigmas_totintra, sigmas_intra12, sigmas_intra13, sigmas_intra14,
            C12_tot, C12_totinter, C12_totintra, C12_intra12, C12_intra13, C12_intra14,C6_tot,
            C6_totinter, C6_totintra, C6_intra12, C6_intra13, C6_intra14);

    // calculate the LJ potentials using the C12 and C6 values above...
    map<IJ, LJpot> ftotLJ;
    map<IJ, LJpot> ftotinterLJ;
    map<IJ, LJpot> ftotintraLJ;
    map<IJ, LJpot> fintra12LJ;
    map<IJ, LJpot> fintra13LJ;
    map<IJ, LJpot> fintra14LJ;
    for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
      double c6 = 4 * epsilons_tot[*it] * pow(sigmas_tot[*it], 6);
      double c12 = 4 * epsilons_tot[*it] * pow(sigmas_tot[*it], 12);
      ftotLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totinter[*it] * pow(sigmas_totinter[*it], 6);
      c12 = 4 * epsilons_totinter[*it] * pow(sigmas_totinter[*it], 12);
      ftotinterLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_totintra[*it] * pow(sigmas_totintra[*it], 6);
      c12 = 4 * epsilons_totintra[*it] * pow(sigmas_totintra[*it], 12);
      ftotintraLJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra12[*it] * pow(sigmas_intra12[*it], 6);
      c12 = 4 * epsilons_intra12[*it] * pow(sigmas_intra12[*it], 12);
      fintra12LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra13[*it] * pow(sigmas_intra13[*it], 6);
      c12 = 4 * epsilons_intra13[*it] * pow(sigmas_intra13[*it], 12);
      fintra13LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      c6 = 4 * epsilons_intra14[*it] * pow(sigmas_intra14[*it], 6);
      c12 = 4 * epsilons_intra14[*it] * pow(sigmas_intra14[*it], 12);
      fintra14LJ.insert(pair<IJ, LJpot > (*it, LJpot(distmin, distmax, distgrid, c12, c6)));
    }
    {
      ofstream fout(fname_LJpot_CG.c_str());
      printPot(fout, ftotLJ, ftotinterLJ, ftotintraLJ, fintra12LJ, fintra13LJ, fintra14LJ);
      fout.close();
    }

    // normalise and print the distribution
    {
      ofstream fout(fname_beadbead_dist.c_str());
      double dgrid = (bondlength_max - bondlength_min) / 200;
      fout.precision(9);
      fout << "#" << setw(19) << "r / nm";
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j();
        fout << scientific << setw(20) << ss.str();
      }
      fout << endl;
      for (int i = 0; i < 200; ++i) {
        double r = (i + 0.5) * dgrid + bondlength_min;
        fout << scientific << setw(20) << r;
        for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
          if (beadbeadDist[*it].nVal() > 0) {
            fout << scientific << setw(20) << (double) beadbeadDist[*it][i] / beadbeadDist[*it].nVal() * 100;
          } else {
            fout << scientific << setw(20) << beadbeadDist[*it][i];
          }
        }
        fout << endl << "#" << endl;
      }
      fout.close();
    }
    
    // print the angle distribution
    {
      ofstream fout("angle.dist");
      double dgrid = 1.0;
      fout.precision(9);
      fout << "#" << setw(19) << "angle / degree";
      for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j() << "-" << it->k();
        fout << scientific << setw(20) << ss.str();
      }
      fout << endl;
      for (int i = 0; i < 180; ++i) {
        double a = (i + 0.5) * dgrid;
        fout << scientific << setw(20) << a;
        for (set<IJK>::const_iterator it = IJKs.begin(); it != IJKs.end(); ++it) {
          if (angleDist[*it].nVal() > 0) {
            fout << scientific << setw(20) << (double) angleDist[*it][i] / angleDist[*it].nVal() * 100;
          } else {
            fout << scientific << setw(20) << (double) angleDist[*it][i];
          }
        }
        fout << endl << "#" << endl;
      }
      fout.close();
    }
    
    // print the dihedral distribution
    {
      ofstream fout("dihedral.dist");
      double dgrid = 1.0;
      fout.precision(9);
      fout << "#" << setw(19) << "dihedal / degree";
      for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j() << "-" << it->k() << "-" << it->l();
        fout << scientific << setw(20) << ss.str();
      }
      fout << endl;
      for (int i = 0; i < 180; ++i) {
        double a = (i + 0.5) * dgrid;
        fout << scientific << setw(20) << a;
        for (set<IJKL>::const_iterator it = IJKLs.begin(); it != IJKLs.end(); ++it) {
          if (dihedralDist[*it].nVal() > 0) {
            fout << scientific << setw(20) <<  (double) dihedralDist[*it][i] / dihedralDist[*it].nVal() * 100;
          } else {
            fout << scientific << setw(20) <<  (double) dihedralDist[*it][i];
          }
        }
        fout << endl << "#" << endl;
      }
      fout.close();
    }

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

namespace cgLJpot {
  
  IJ::IJ(int i, int j) {
    if (j < i) {
      int t = i;
      i = j;
      j = t;
    }
    I = i;
    J = j;
  }

  IJ::IJ(const IJ &ij) {
    I = ij.I;
    J = ij.J;
  }

  int IJ::i() const {
    return I;
  }

  int IJ::j() const {
    return J;
  }

  void IJ::setValues(int i, int j) {
    if (j < i) {
      int t = i;
      i = j;
      j = t;
    }
    I = i;
    J = j;
  }
  
  IJK::IJK(int i, int j, int k) {
    if (k < i) {
      int t = i;
      i = k;
      k = t;
    }
    I = i;
    J = j;
    K = k;
  }

  IJK::IJK(const IJK &ijk) {
    I = ijk.I;
    J = ijk.J;
    K = ijk.K;
  }

  int IJK::i() const {
    return I;
  }

  int IJK::j() const {
    return J;
  }
  
  int IJK::k() const {
    return K;
  }

  void IJK::setValues(int i, int j, int k) {
    if (k < i) {
      int t = i;
      i = k;
      k = t;
    }
    I = i;
    J = j;
    K = k;
  }
  
  IJKL::IJKL(int i, int j, int k, int l) {
    if (k < j) {
      int t = j;
      j = k;
      k = t;
      t = i;
      i = l;
      l = t;
    }
    I = i;
    J = j;
    K = k;
    L = l;
  }

  IJKL::IJKL(const IJKL &ijkl) {
    I = ijkl.I;
    J = ijkl.J;
    K = ijkl.K;
    L = ijkl.L;
  }

  int IJKL::i() const {
    return I;
  }

  int IJKL::j() const {
    return J;
  }
  
  int IJKL::k() const {
    return K;
  }
  
  int IJKL::l() const {
    return L;
  }

  void IJKL::setValues(int i, int j, int k, int l) {
    if (k < j) {
      int t = j;
      j = k;
      k = t;
      t = i;
      i = l;
      l = t;
    }
    I = i;
    J = j;
    K = k;
    L = l;
  }

  LJpot::LJpot(double min_, double max_, int grid) {
    lj.resize(grid);
    count.resize(grid);
    min = min_;
    max = max_;
    dgrid = (max - min) / grid;
    C12 = -1.0;
    C6 = -1.0;
    for (int i = 0; i < grid; ++i) {
      lj[i] = 0.0;
      count[i] = 0;
    }
  }

  LJpot::LJpot(double min_, double max_, int grid, double c12, double c6) {
    lj.resize(grid);
    count.resize(grid);
    min = min_;
    max = max_;
    dgrid = (max - min) / grid;
    C12 = c12;
    C6 = c6;
    for (int i = 0; i < grid; ++i) {
      double R = r(i);
      double R3 = R * R * R;
      double R6 = R3 * R3;
      lj[i] = (C12 / R6 - C6) / R6;
      count[i] = 1;
    }
  }

  LJpot::LJpot(const LJpot &ljp) {
    lj = ljp.lj;
    count = ljp.count;
    min = ljp.min;
    max = ljp.max;
    dgrid = ljp.dgrid;
    C12 = ljp.C12;
    C6 = ljp.C6;
  }

  void LJpot::add(double pos, double val) {
    if (min <= pos && pos < max) {
      int r = floor((pos - min) / dgrid);
      count[r] += 1;
      lj[r] += val;
    }
  }

  void LJpot::unify(const LJpot & ljp) {
    assert(lj.size() == ljp.lj.size());
    for (unsigned int i = 0; i < lj.size(); ++i) {
      lj[i] += ljp.lj[i];
      count[i] += ljp.count[i];
    }
  }

  double LJpot::r(unsigned int i) {
    assert(i < lj.size());
    return (i + 0.5) * dgrid + min;
  }

  double LJpot::r(unsigned int i) const {
    assert(i < lj.size());
    return (i + 0.5) * dgrid + min;
  }

  double LJpot::pot(unsigned int i) {
    assert(i < lj.size());
    return count[i] == 0 ? 0.0 : lj[i] / count[i];
  }

  double LJpot::get_min() {
    return min;
  }

  double LJpot::get_max() {
    return max;
  }

  int LJpot::get_grid() {
    return lj.size();
  }

  double LJpot::rmin() {
    // find the grid with maximal pot first
    double i_max = 0;
    double pot_max = pot(0);
    for (int i = 1; i < get_grid(); ++i) {
      if (pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double r_min = r(i_max);
    double pot_min = pot(i_max);
    for (int i = i_max + 1; i < get_grid(); ++i) {
      if (pot(i) < pot_min) {
        r_min = r(i);
        pot_min = pot(i);
      }
    }
    return r_min;
  }

  double LJpot::potmin() {
    // find the grid with maximal pot first
    double i_max = 0;
    double pot_max = pot(0);
    for (int i = 1; i < get_grid(); ++i) {
      if (pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double r_min = r(i_max);
    double pot_min = pot(i_max);
    for (int i = i_max + 1; i < get_grid(); ++i) {
      if (pot(i) < pot_min) {
        r_min = r(i);
        pot_min = pot(i);
      }
    }
    return pot_min;
  }

  bead::bead(gcore::System& sys, GromosForceField &groff, int mom, int bnum, set<IJ> &ij, double min, double max, double grid) {
    gff = &groff;
    atoms.setSystem(sys);
    memberOfMol = mom;
    beadnum = bnum;
    set<IJ>::const_iterator it = ij.begin();
    for (; it != ij.end(); ++it) {
      totLJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totinterLJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      totintraLJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra12LJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra13LJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
      intra14LJ.insert(pair<IJ, LJpot > (*it, LJpot(min, max, grid)));
    }
  }

  bead::bead(bead const &b) {
    gff = b.gff;
    atoms = b.atoms;
    centre = b.centre;
    memberOfMol = b.memberOfMol;
    beadnum = b.beadnum;
    totLJ = b.totLJ;
    totinterLJ = b.totinterLJ;
    totintraLJ = b.totintraLJ;
    intra12LJ = b.intra12LJ;
    intra13LJ = b.intra13LJ;
    intra14LJ = b.intra14LJ;
  }

  int bead::addAtom(int m, int a) {
    return atoms.addAtom(m, a);
  }

  int bead::size() {
    return atoms.size();
  }

  Vec bead::cog(Boundary *pbc, System &sys) {
    assert(atoms.size() > 0);
    Vec cog(atoms.pos(0));
    for (int i = 1; i < atoms.size(); ++i) {
      cog += pbc->nearestImage(atoms.pos(0), atoms.pos(i), sys.box());
    }
    cog /= atoms.size();
    centre = cog;
    return cog;
  }

  Vec bead::com(Boundary *pbc, System &sys) {
    assert(atoms.size() > 0);
    double sumMass = atoms.mass(0);
    Vec com(atoms.pos(0) * sumMass);
    for (int i = 1; i < atoms.size(); ++i) {
      double mass = atoms.mass(i);
      sumMass += mass;
      com += pbc->nearestImage(atoms.pos(0), atoms.pos(i), sys.box()) * mass;
    }
    com = com / (sumMass);
    centre = com;
    return com;
  }

  Vec bead::pos() {
    return centre;
  }

  map<IJ, LJpot> bead::get_totLJ() {
    return totLJ;
  }

  map<IJ, LJpot> bead::get_totinterLJ() {
    return totinterLJ;
  }

  map<IJ, LJpot> bead::get_totintraLJ() {
    return totintraLJ;
  }

  map<IJ, LJpot> bead::get_intra12LJ() {
    return intra12LJ;
  }

  map<IJ, LJpot> bead::get_intra13LJ() {
    return intra13LJ;
  }

  map<IJ, LJpot> bead::get_intra14LJ() {
    return intra14LJ;
  }

  int bead::mol() {
    return memberOfMol;
  }

  int bead::beadNum() {
    return beadnum;
  }

  double bead::calcLJ(bead &b, Boundary *pbc, System &sys) {
    double LJsum = 0.0;
    for (int i1 = 0; i1 < size(); ++i1) {
      for (int i2 = 0; i2 < b.size(); ++i2) {
        // the stuff for atom 1 is also defined in here since
        // the atoms at1 and a2 might be interchanged sometimes which has to be
        // undone before continuing with the nest at2
        int m1 = atoms.mol(i1);
        int a1 = atoms.atom(i1);
        int gnum1 = atoms.gromosAtom(i1);
        AtomTopology at1 = sys.mol(m1).topology().atom(a1);
        int m2 = b.atoms.mol(i2);
        int a2 = b.atoms.atom(i2);
        int gnum2 = b.atoms.gromosAtom(i2);
        AtomTopology at2 = sys.mol(m2).topology().atom(a2);
        // to memorize if the inter-bead atoms are excluded from each other
        bool excluded = false;
        bool excluded14 = false;
        // check if and which (normal, special 1,4; only happens if the atoms are
        // located in the same molecule) LJ interaction is calculated
        if (m1 == m2) {
          // make sure at1 is the atom first listed in the topology
          if (gnum2 < gnum1) {
            int mtmp = m1;
            m1 = m2;
            m2 = mtmp;
            int atmp = a1;
            a1 = a2;
            a2 = atmp;
            int gnumtmp = gnum1;
            gnum1 = gnum2;
            gnum2 = gnumtmp;
            AtomTopology attmp(at1);
            at1 = at2;
            at2 = attmp;
          }
          // check if and which (normal, special 1,4) LJ interaction is calculated
          for (int e = 0; e < at1.exclusion().size(); ++e) {
            if (at1.exclusion().atom(e) == a2) {
              excluded = true;
              break;
            }
          }
          for (int e = 0; e < at1.exclusion14().size(); ++e) {
            if (at1.exclusion14().atom(e) == a2) {
              excluded14 = true;
              break;
            }
          }
        }
        // do the appropriate LJ potential energy calculation
        double c12;
        double c6;
        int iac1 = at1.iac();
        int iac2 = at2.iac();
        if (excluded) {
          continue;
        } else if (excluded14) {
          c12 = gff->ljType(AtomPair(iac1, iac2)).cs12();
          c6 = gff->ljType(AtomPair(iac1, iac2)).cs6();
        } else {
          c12 = gff->ljType(AtomPair(iac1, iac2)).c12();
          c6 = gff->ljType(AtomPair(iac1, iac2)).c6();
        }
        double r2a = (atoms.pos(i1) -
                (pbc->nearestImage(atoms.pos(i1), b.atoms.pos(i2), sys.box()))).abs2();
        double r6a = r2a * r2a * r2a;
        LJsum += (c12 / r6a - c6) / r6a;
      }
    }
    return LJsum;
  }

  void bead::addLJtot(const IJ &ij, double r, const double &lj) {
    totLJ[ij].add(r, lj);
  }

  void bead::addLJtotinter(const IJ &ij, double r, const double &lj) {
    totinterLJ[ij].add(r, lj);
  }

  void bead::addLJtotintra(const IJ &ij, double r, const double &lj) {
    totintraLJ[ij].add(r, lj);
  }

  void bead::addLJintra12(const IJ &ij, double r, const double &lj) {
    intra12LJ[ij].add(r, lj);
  }

  void bead::addLJintra13(const IJ &ij, double r, const double &lj) {
    intra13LJ[ij].add(r, lj);
  }

  void bead::addLJintra14(const IJ &ij, double r, const double &lj) {
    intra14LJ[ij].add(r, lj);
  }

  bool operator<(const IJ &ij1, const IJ &ij2) {
    if (ij1.i() < ij2.i() || (ij1.i() == ij2.i() && ij1.j() < ij2.j())) {
      return true;
    }
    return false;
  }
  
  bool operator<(const IJK &ijk1, const IJK &ijk2) {
    if (ijk1.i() < ijk2.i() || (ijk1.i() == ijk2.i() && ijk1.k() < ijk2.k())) {
      return true;
    }
    return false;
  }
  
  bool operator<(const IJKL &ijkl1, const IJKL &ijkl2) {
    if (ijkl1.j() < ijkl2.j() || 
            (ijkl1.j() == ijkl2.j() && ijkl1.k() < ijkl2.k()) || 
            (ijkl1.j() == ijkl2.j() && ijkl1.k() == ijkl2.k() && ijkl1.i() < ijkl2.i()) ||
            (ijkl1.j() == ijkl2.j() && ijkl1.k() == ijkl2.k() && ijkl1.i() == ijkl2.i() && ijkl1.l() < ijkl2.l())) {
      return true;
    }
    return false;
  }

  void printPot(ostream &os, map<IJ, LJpot> &totLJpot, map<IJ, LJpot> &totLJinter,
          map<IJ, LJpot> &totLJintra, map<IJ, LJpot> &intraLJ12, map<IJ, LJpot> &intraLJ13,
          map<IJ, LJpot> &intraLJ14) {
    vector<string> header;
    vector<string> potentialnames;
    potentialnames.push_back("totLJpot");
    potentialnames.push_back("totLJinter");
    potentialnames.push_back("totLJintra");
    potentialnames.push_back("totLJ12");
    potentialnames.push_back("totLJ13");
    potentialnames.push_back("totLJ14");
    header.push_back("radius / nm");
    for (unsigned int i = 0; i < potentialnames.size(); ++i) {
      for (map<IJ, LJpot>::const_iterator it = totLJpot.begin(); it != totLJpot.end(); ++it) {
        stringstream ss;
        ss << potentialnames[i] << "_" << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    // print the header
    os << "#";
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        os << setw(19) << header[i];
      } else {
        os << setw(20) << header[i];
      }
    }
    os << endl;
    os.precision(9);
    map<IJ, LJpot>::const_iterator iter = totLJpot.begin();
    for (int i = 0; i < totLJpot.begin()->second.get_grid(); ++i) {
      os << scientific << setw(20) << iter->second.r(i);
      for (map<IJ, LJpot>::const_iterator it = totLJpot.begin(); it != totLJpot.end(); ++it) {
        os << scientific << setw(20) << totLJpot[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = totLJinter.begin(); it != totLJinter.end(); ++it) {
        os << scientific << setw(20) << totLJinter[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = totLJintra.begin(); it != totLJintra.end(); ++it) {
        os << scientific << setw(20) << totLJintra[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ12.begin(); it != intraLJ12.end(); ++it) {
        os << scientific << setw(20) << intraLJ12[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ13.begin(); it != intraLJ13.end(); ++it) {
        os << scientific << setw(20) << intraLJ13[it->first].pot(i);
      }
      for (map<IJ, LJpot>::const_iterator it = intraLJ14.begin(); it != intraLJ14.end(); ++it) {
        os << scientific << setw(20) << intraLJ14[it->first].pot(i);
      }
      os << endl;
    }
    
    /*
    os << endl << endl;
    for (map<IJ, LJpot>::iterator it = totLJpot.begin(); it != totLJpot.end(); it++) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    for (map<IJ, LJpot>::iterator it = totLJinter.begin(); it != totLJinter.end(); ++it) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    for (map<IJ, LJpot>::iterator it = totLJintra.begin(); it != totLJintra.end(); ++it) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    for (map<IJ, LJpot>::iterator it = intraLJ12.begin(); it != intraLJ12.end(); ++it) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    for (map<IJ, LJpot>::iterator it = intraLJ13.begin(); it != intraLJ13.end(); ++it) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    for (map<IJ, LJpot>::iterator it = intraLJ14.begin(); it != intraLJ14.end(); ++it) {
      os << scientific << setw(20) << it->second.rmin() << setw(20) << it->second.potmin() << endl << endl;
    }
    */

  }

  void printTitleBlock(vector<int> &beadsizes, AtomSpecifier &allAtoms) {
    cout << "TITLE\n";
    cout << "   number of beads per molecule: " << beadsizes.size() << endl;
    cout << "   number of atoms per bead: ";
    for (unsigned int b = 0; b < beadsizes.size(); ++b) {
      cout << beadsizes[b] << " ";
    }
    cout << endl;
    cout << "   the molecule: |";
    int at = 0;
    for (unsigned int b = 0; b < beadsizes.size(); ++b) {
      for (int a = 0; a < beadsizes[b]; ++a) {
        cout << allAtoms.name(at);
        if (a < beadsizes[b] - 1) {
          cout << "--";
        } else if (b < beadsizes.size() - 1) {
          cout << "-|-";
        }
        at++;
      }
    }
    cout << "|\n";
    time_t rawtime;
    time(&rawtime);
    cout << "   Timestamp: " << ctime(&rawtime) << "END\n";
  }

  void calcEpsSigma(map<IJ, double> &epsilons, map<IJ, double> &sigmas, map<IJ, LJpot> &LJ) {
    if (sigmas.size() != 0) {
      sigmas.clear();
    }
    if (epsilons.size() != 0) {
      epsilons.clear();
    }
    for (map<IJ, LJpot>::const_iterator it = LJ.begin(); it != LJ.end(); ++it) {
      double sigma = LJ[it->first].rmin() / pow(2.0, 1.0 / 6.0);
      double epsilon = -LJ[it->first].potmin();
      sigmas.insert(pair<IJ, double> (it->first, sigma));
      epsilons.insert(pair<IJ, double>(it->first, epsilon));
    }
  }
  
  void calcC12C6(map<IJ, double> &C12s, map<IJ, double> &C6s, map<IJ, double> &epsilons, map<IJ, double> &sigmas) {
    for(map<IJ, double>::const_iterator it =epsilons.begin(); it != epsilons.end(); ++it) {
      C6s.insert(pair<IJ, double>(it->first, 4 * epsilons[it->first] * pow(sigmas[it->first], 6)));
      C12s.insert(pair<IJ, double>(it->first, 4 * epsilons[it->first] * pow(sigmas[it->first], 12)));
    }
  }
  
  void printLennardJonesParamters(map<IJ, double> &epsilons_tot,
          map<IJ, double> &epsilons_totinter,
          map<IJ, double> &epsilons_totintra,
          map<IJ, double> &epsilons_intra12,
          map<IJ, double> &epsilons_intra13,
          map<IJ, double> &epsilons_intra14,
          map<IJ, double> &sigmas_tot,
          map<IJ, double> &sigmas_totinter,
          map<IJ, double> &sigmas_totintra,
          map<IJ, double> &sigmas_intra12,
          map<IJ, double> &sigmas_intra13,
          map<IJ, double> &sigmas_intra14,
          map<IJ, double> &C12_tot,
          map<IJ, double> &C12_totinter,
          map<IJ, double> &C12_totintra,
          map<IJ, double> &C12_intra12,
          map<IJ, double> &C12_intra13,
          map<IJ, double> &C12_intra14,
          map<IJ, double> &C6_tot,
          map<IJ, double> &C6_totinter,
          map<IJ, double> &C6_totintra,
          map<IJ, double> &C6_intra12,
          map<IJ, double> &C6_intra13,
          map<IJ, double> &C6_intra14) {
    vector<string> header;
    vector<string> names;
    names.push_back("eps_tot_");
    names.push_back("eps_totinter_");
    names.push_back("eps_totintra_");
    names.push_back("eps_intra12_");
    names.push_back("eps_intra13_");
    names.push_back("eps_intra14_");
    cout << "LENNARD-JONES\n";
    cout << "# var_xy_i-j:\n";
    cout << "#   var ... eps:      epsilon\n";
    cout << "#       ... sig:      sigma\n";
    cout << "#   C12 ... C12 LJ-potential parameter\n";
    cout << "#   C6  ... C6 LJ-potential parameter\n";
    cout << "#   xy  ... tot:      total LJ-potential energy\n";
    cout << "#       ... totinter: total inter-molecular LJ potential\n";
    cout << "#       ... totintra: total intra-molecular LJ potential\n";
    cout << "#       ... intra12:  total 12-intra-molecular LJ potential\n";
    cout << "#       ... intra13:  total 13-intra-molecular LJ potential\n";
    cout << "#       ... intra14:  total 14-intra-molecular LJ potential\n";
    cout << "#   i   ... interaction from bead of size i to ...\n";
    cout << "#   j   ... bead of size j\n";
    cout << "#\n";
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::const_iterator it = epsilons_tot.begin(); it != epsilons_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    cout.precision(9);
    cout.setf(ios::scientific);
    for (map<IJ, double>::const_iterator it = epsilons_tot.begin(); it != epsilons_tot.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_totinter.begin(); it != epsilons_totinter.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_totintra.begin(); it != epsilons_totintra.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra12.begin(); it != epsilons_intra12.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra13.begin(); it != epsilons_intra13.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = epsilons_intra14.begin(); it != epsilons_intra14.end(); ++it) {
      cout << setw(20) << it->second;
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("sig_tot_");
    names.push_back("sig_totinter_");
    names.push_back("sig_totintra_");
    names.push_back("sig_intra12_");
    names.push_back("sig_inter13_");
    names.push_back("sig_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = sigmas_tot.begin(); it != sigmas_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = sigmas_tot.begin(); it != sigmas_tot.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_totinter.begin(); it != sigmas_totinter.end(); ++it) {
      cout << setw(20) << it->second;
    }
   for (map<IJ, double>::const_iterator it = sigmas_totintra.begin(); it != sigmas_totintra.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra12.begin(); it != sigmas_intra12.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra13.begin(); it != sigmas_intra13.end(); ++it) {
      cout << setw(20) << it->second;
    }
    for (map<IJ, double>::const_iterator it = sigmas_intra14.begin(); it != sigmas_intra14.end(); ++it) {
      cout << setw(20) << it->second;
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("C12_tot_");
    names.push_back("C12_totinter_");
    names.push_back("C12_totintra_");
    names.push_back("C12_intra12_");
    names.push_back("C12_inter13_");
    names.push_back("C12_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
      cout << setw(20) << C12_tot[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_totinter.begin(); it != C12_totinter.end(); ++it) {
      cout << setw(20) << C12_totinter[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_totintra.begin(); it != C12_totintra.end(); ++it) {
      cout << setw(20) << C12_totintra[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra12.begin(); it != C12_intra12.end(); ++it) {
      cout << setw(20) << C12_intra12[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra13.begin(); it != C12_intra13.end(); ++it) {
      cout << setw(20) << C12_intra13[it->first];
    }
    for (map<IJ, double>::const_iterator it = C12_intra14.begin(); it != C12_intra14.end(); ++it) {
      cout << setw(20) << C12_intra14[it->first];
    }
    cout << "\n#\n";
    header.clear();
    names.clear();
    names.push_back("C6_tot_");
    names.push_back("C6_totinter_");
    names.push_back("C6_totintra_");
    names.push_back("C6_intra12_");
    names.push_back("C6_inter13_");
    names.push_back("C6_inter14_");
    for (unsigned int i = 0; i < names.size(); ++i) {
      for (map<IJ, double>::iterator it = C12_tot.begin(); it != C12_tot.end(); ++it) {
        stringstream ss;
        ss << names[i] << it->first.i() << "-" << it->first.j();
        header.push_back(ss.str());
      }
    }
    for (unsigned int i = 0; i < header.size(); ++i) {
      if (i == 0) {
        cout << "#" << setw(19) << header[i];
      } else {
        cout << setw(20) << header[i];
      }
    }
    cout << endl;
    for (map<IJ, double>::const_iterator it = C6_tot.begin(); it != C6_tot.end(); ++it) {
      cout << setw(20) << C6_tot[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_totinter.begin(); it != C6_totinter.end(); ++it) {
      cout << setw(20) << C6_totinter[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_totintra.begin(); it != C6_totintra.end(); ++it) {
      cout << setw(20) << C6_totintra[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra12.begin(); it != C6_intra12.end(); ++it) {
      cout << setw(20) << C6_intra12[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra13.begin(); it != C6_intra13.end(); ++it) {
      cout << setw(20) << C6_intra13[it->first];
    }
    for (map<IJ, double>::const_iterator it = C6_intra14.begin(); it != C6_intra14.end(); ++it) {
      cout << setw(20) << C6_intra14[it->first];
    }
    cout << "\nEND\n";
  }

}

