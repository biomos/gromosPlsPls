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
 * @date 29.06.2011
 *
 * Program cgLJpot ...
 * 
 * NOTE: the current program version does only work for solute atoms/molecules.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@method</td><td>&lt;method to goarse grain: atomic or molecular&gt; </td></tr>
 * <tr><td> \@dist</td><td>&lt;min max ngrid&gt; </td></tr>
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

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

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

namespace cgLJpot {

  class IJ {
  private:
    int I;
    int J;

  public:
    IJ(int i = -1, int j = -1);
    IJ(const IJ &ij);
    void setValues(int i, int j);
    int i() const;
    int j() const;
  };
  
  bool operator<(const IJ &ij1, const IJ &ij2);

  class LJpot {
  private:
    vector<double> lj;
    vector<int> count;
    double dgrid;
    double min;
    double max;
    double C12;
    double C6;

  public:
    LJpot(double min_, double max_, int grid, double c12, double c6);
    LJpot(double min_ = 0.0, double max_ = 2.0, int grid = 200);
    LJpot(const LJpot &ljp);
    void init(double min_ = 0.0, double max_ = 2.0, int grid = 200);
    void add(double pos, double val);
    void print(ostream &os);
    void unify(const LJpot &ljp);
    double r(unsigned int i);
    double r(unsigned int i) const;
    double pot(unsigned int i);
    double get_min();
    double get_max();
    int get_grid();
    int size();
    // returns the radius with the minimal LJ potential energy
    double rmin();
    double potmin();
  };

  class bead {
  private:
    gcore::GromosForceField *gff;
    utils::AtomSpecifier atoms;
    gmath::Vec centre;
    // molecule number the bead belongs to
    int memberOfMol;
    // bead number within the molecule, starting at 0
    int beadnum;
    map<IJ, LJpot> totLJ;
    map<IJ, LJpot> totinterLJ;
    map<IJ, LJpot> totintraLJ;
    map<IJ, LJpot> intra12LJ;
    map<IJ, LJpot> intra13LJ;
    map<IJ, LJpot> intra14LJ;

  public:
    // Constructor, needs to know about the system
    bead(gcore::System &sys, gcore::GromosForceField &groff, int mom, int bnum, set<IJ> &ij, double min = 0.0, double max = 2.0, double grid = 200);
    bead(bead const &b);
    // Destructor
    ~bead() {};
    // add an atom to the bead:
    // m: moleucle number
    // a: atom number
    int addAtom(int m, int a);
    // return the bead size
    int size();
    // calculates and returns the centre of geometry; setc the centre to cog
    gmath::Vec cog(bound::Boundary *pbc, gcore::System &sys);
    // calculates and returns the centre of mass; sets the centre to com
    gmath::Vec com(bound::Boundary *pbc, gcore::System &sys);
    // Accessor to the position of the bead (com or cog)
    gmath::Vec pos();
    // add a value to the total LJ potential energy
    void addLJtot(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // add a value to the total intermolecular LJ potential energy
    void addLJtotinter(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // add a value to the total intramolecular LJ potential energy
    void addLJtotintra(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // add a value to the total 12-intermolecular LJ potential energy
    void addLJintra12(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // add a value to the total 13-intermolecular LJ potential energy
    void addLJintra13(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // add a value to the total 14-intermolecular LJ potential energy
    void addLJintra14(const IJ &ij, double r, const double &lj, double min, double max, int grid);
    // calculate the LJ interaction to another bead
    double calcLJ(bead &b, bound::Boundary *pbc, gcore::System &sys);
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_totLJ();
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_totinterLJ();
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_totintraLJ();
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_intra12LJ();
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_intra13LJ();
    // Accessor: returns the intermolecular LJpot energy
    map<IJ, LJpot> get_intra14LJ();
    // Accessor: returns the molecule number the bead belongs to
    int mol();
    // Accessor: returns the bead number within the molecule
    int beadNum();
  };
  
  LJpot leastSquareLJ(LJpot &ljp);
  
  void printPot(ostream &os, std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJpot,
          std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJinter,
          std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJintra,
          std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJ12,
          std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJ13,
          std::map<cgLJpot::IJ, cgLJpot::LJpot> &totLJ14);

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
  knowns << "topo" << "method" << "beads" << "pbc" << "trc" << "dist" << "verbose" << "outfit" << "outbonddist" << "hvap";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@method        <method to goarse grain: atomic or molecular>\n";
  usage += "\t[@dist         <min max ngrid>]\n";
  usage += "\t@beads         <number of atoms per bead (atomic)> or\n";
  usage += "\t               <sequence of bead size within one molecule (molecular)>\n";
  usage += "\t[@pbc          <boundary type (read from GENBOX block if not specified)> [<gather method>]]\n";
  usage += "\t[@hvap         <experimental heat of vaporisation / kJ mol^(-1)>]\n";
  usage += "\t[@outfit       <output file name for fitted LJ potentials>\n";
  usage += "\t[@outbonddist  <output file for bead-bead bond distributions\n";
  usage += "\t@trc           <simulation trajectory or coordinate file>\n";
  usage += "\t@verbose       gives more verbose output";

  try {
    Arguments args(argc, argv, knowns, usage);
    
    // is there a heat of vaporisation to be used to adept epsilon?
    double hvap;
    if (args.count("hvap") >= 0) {
      if (args.count("outfit") >= 0) {
        throw gromos::Exception(argv[0], "@outfit cannot be used at the same time as @hvap");
      }
      if (args.count("hvap") != 1) {
        throw gromos::Exception(argv[0], "there must be exactly one argument for @hvap");
      } else {
        hvap = atof(args["hvap"].c_str());
      }
    }
    
    // more than only standard output?
    bool printfit;
    string fitout;
    if(args.count("outfit") == -1) {
      printfit = false;
    } else if(args.count("outfit") == 1) {
      fitout = args["outfit"];
      printfit = true;
    } else {
      throw gromos::Exception(argv[0], "wrong number of arguments for @outfit");
    }
    bool printBeadBeadDist;
    string BeadBeadDistOut;
    if(args.count("outbonddist") == -1) {
      printBeadBeadDist = false;
    } else if(args.count("outbonddist") == 1) {
      BeadBeadDistOut = args["outbonddist"];
      printBeadBeadDist = true;
    } else {
      throw gromos::Exception(argv[0], "wrong number of arguments for @outbonddist");
    }
    
    // verbose or not?
    bool verbose = false;
    if (args.count("verbose") >= 0) {
      verbose = true;
    }

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
    for (unsigned int bs1 = 0; bs1 < beadsizes.size(); ++bs1) {
      for (unsigned int bs2 = bs1; bs2 < beadsizes.size(); ++bs2) {
        IJ ij(beadsizes[bs1], beadsizes[bs2]);
        if (IJs.find(ij) == IJs.end()) {
          IJs.insert(ij);
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
    if (method == "molecular" && verbose) {
      cout << "# Molecular setup of the beads\n";
      cout << "# ============================\n";
      cout << "#\n";
      cout << "# number of beads per molecule: " << beadsizes.size() << endl;
      cout << "# number of atoms per bead: ";
      for (unsigned int b = 0; b < beadsizes.size(); ++b) {
        cout << beadsizes[b] << " ";
      }
      cout << endl;
      cout << "# the molecule: |";
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
      cout << "|\n#\n";
    }

    // a map of distributions (for each IJ one) to remember the intramolecular
    // neighboring bead-bead distance (min and max automatically set based on the
    // first configuration of the trajectories
    map<IJ, Distribution> beadbeadDist;
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
        if(sys.box().L().abs2() < L) {
          L = sys.box().L().abs2();
        }
        if(sys.box().M().abs2() < L) {
          L = sys.box().M().abs2();
        }
        if(distmax > (sqrt(L) / 2.0)) {
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
                if(lj == 0.0) {
                  continue;
                }
                r = std::sqrt(r2);
                // and add the result to the two beads (corresponding potentials)
#ifdef OMP
#pragma omp critical
#endif               
                {
                  beads[b1].addLJtot(ij, r, lj, distmin, distmax, distgrid);
                  beads[b2].addLJtot(ij, r, lj, distmin, distmax, distgrid);
                  if (beads[b1].mol() != beads[b2].mol()) { // intermolecular LJ potential
                    beads[b1].addLJtotinter(ij, r, lj, distmin, distmax, distgrid);
                    beads[b2].addLJtotinter(ij, r, lj, distmin, distmax, distgrid);
                  } else { // intramolecular LJ potential energy
                    beads[b1].addLJtotintra(ij, r, lj, distmin, distmax, distgrid);
                    beads[b2].addLJtotintra(ij, r, lj, distmin, distmax, distgrid);
                    // 12-LJ pot (neighboring beads)
                    if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 1) {
                      beads[b1].addLJintra12(ij, r, lj, distmin, distmax, distgrid);
                      beads[b2].addLJintra12(ij, r, lj, distmin, distmax, distgrid);
                      // add the distance to the distribution
                      beadbeadDist[ij].add(r);
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 2) {
                      beads[b1].addLJintra13(ij, r, lj, distmin, distmax, distgrid);
                      beads[b2].addLJintra13(ij, r, lj, distmin, distmax, distgrid);
                    } else if (abs(beads[b1].beadNum() - beads[b2].beadNum()) == 3) {
                      beads[b1].addLJintra14(ij, r, lj, distmin, distmax, distgrid);
                      beads[b2].addLJintra14(ij, r, lj, distmin, distmax, distgrid);
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
      totLJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
      totinterLJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
      totintraLJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
      intra12LJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
      intra13LJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
      intra14LJ.insert(pair<IJ, LJpot> (*it, LJpot(distmin, distmax, distgrid)));
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
    printPot(cout, totLJ, totinterLJ, totintraLJ, intra12LJ, intra13LJ, intra14LJ);
    
    // calculate and print the resulting LJ pot based on the heat of vaporisation, if requested
    if(args.count("hvap") > 0) {
      map<IJ, double> sigmas_tot;
      map<IJ, double> sigmas_totinter;
      map<IJ, double> sigmas_totintra;
      map<IJ, double> sigmas_intra12;
      map<IJ, double> sigmas_intra13;
      map<IJ, double> sigmas_intra14;
      map<IJ, double> epsilon_tot;
      map<IJ, double> epsilon_totinter;
      map<IJ, double> epsilon_totintra;
      map<IJ, double> epsilon_intra12;
      map<IJ, double> epsilon_intra13;
      map<IJ, double> epsilon_intra14;
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        double sigma = totLJ[*it].rmin() / pow(2.0, 1.0/6.0);
        double epsilon = -totLJ[*it].potmin();
        sigmas_tot.insert(pair<IJ, double> (*it, sigma));
        epsilon_tot.insert(pair<IJ, double>(*it, epsilon));
        sigma = totinterLJ[*it].rmin() / pow(2.0, 1.0/6.0);
        epsilon = -totinterLJ[*it].potmin();
        sigmas_totinter.insert(pair<IJ, double> (*it, sigma));
        epsilon_totinter.insert(pair<IJ, double>(*it, epsilon));
        sigma = totintraLJ[*it].rmin() / pow(2.0, 1.0/6.0);
        epsilon  = -totinterLJ[*it].potmin();
        sigmas_totintra.insert(pair<IJ, double> (*it, sigma));
        epsilon_totintra.insert(pair<IJ, double> (*it, epsilon));
        sigma = intra12LJ[*it].rmin() / pow(2.0, 1.0/6.0);
        epsilon = -intra12LJ[*it].potmin();
        sigmas_intra12.insert(pair<IJ, double> (*it, sigma));
        epsilon_intra12.insert(pair<IJ, double> (*it, epsilon));
        sigma = intra13LJ[*it].rmin() / pow(2.0, 1.0/6.0);
        epsilon = -intra13LJ[*it].potmin();
        sigmas_intra13.insert(pair<IJ, double> (*it, sigma));
        epsilon_intra13.insert(pair<IJ, double> (*it, epsilon));
        sigma = intra14LJ[*it].rmin() / pow(2.0, 1.0/6.0);
        epsilon = -intra14LJ[*it].potmin();
        sigmas_intra14.insert(pair<IJ, double> (*it, sigma));
        epsilon_intra14.insert(pair<IJ, double>(*it, epsilon));
      }
      vector<string> header;
      vector<string> names;
      names.push_back("eps_tot_");
      names.push_back("eps_totinter_");
      names.push_back("eps_totintra_");
      names.push_back("eps_intra12_");
      names.push_back("eps_inter13_");
      names.push_back("eps_inter14_");
      for (unsigned int i = 0; i < names.size(); ++i) {
        for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
          stringstream ss;
          ss << names[i] << it->i() << "-" << it->j();
          header.push_back(ss.str());
        }
      }
      for(unsigned int i = 0; i < header.size(); ++i) {
        if (i == 0) {
          cout << "#" << setw(19) << header[i];
        } else {
          cout << setw(20) << header[i];
        }
      }
      cout << endl;
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_tot[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_totinter[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_totintra[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_intra12[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_intra13[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << epsilon_intra14[*it];
      }
      cout << endl <<"#" << endl;
      header.clear();
      names.clear();
      names.push_back("sig_tot_");
      names.push_back("sig_totinter_");
      names.push_back("sig_totintra_");
      names.push_back("sig_intra12_");
      names.push_back("sig_inter13_");
      names.push_back("sig_inter14_");
      for (unsigned int i = 0; i < names.size(); ++i) {
        for (set<IJ>::iterator it = IJs.begin(); it != IJs.end(); ++it) {
          stringstream ss;
          ss << names[i] << it->i() << "-" << it->j();
          header.push_back(ss.str());
        }
      }
      for(unsigned int i = 0; i < header.size(); ++i) {
        if (i == 0) {
          cout << "#" << setw(19) << header[i];
        } else {
          cout << setw(20) << header[i];
        }
      }
      cout << endl;
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_tot[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_totinter[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_totintra[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_intra12[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_intra13[*it];
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << sigmas_intra14[*it];
      }
      cout << endl << "#" << endl;
      header.clear();
      names.clear();
      names.push_back("C12_tot_");
      names.push_back("C12_totinter_");
      names.push_back("C12_totintra_");
      names.push_back("C12_intra12_");
      names.push_back("C12_inter13_");
      names.push_back("C12_inter14_");
      for (unsigned int i = 0; i < names.size(); ++i) {
        for (set<IJ>::iterator it = IJs.begin(); it != IJs.end(); ++it) {
          stringstream ss;
          ss << names[i] << it->i() << "-" << it->j();
          header.push_back(ss.str());
        }
      }
      for(unsigned int i = 0; i < header.size(); ++i) {
        if (i == 0) {
          cout << "#" << setw(19) << header[i];
        } else {
          cout << setw(20) << header[i];
        }
      }
      cout << endl;
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_tot[*it] * pow(sigmas_tot[*it], 12);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_totinter[*it] * pow(sigmas_totinter[*it], 12);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_totintra[*it] * pow(sigmas_totintra[*it], 12);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra12[*it] * pow(sigmas_intra12[*it], 12);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra13[*it] * pow(sigmas_intra13[*it], 12);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra14[*it] * pow(sigmas_intra14[*it], 12);
      }
      cout << endl << "#" << endl;
      header.clear();
      names.clear();
      names.push_back("C6_tot_");
      names.push_back("C6_totinter_");
      names.push_back("C6_totintra_");
      names.push_back("C6_intra12_");
      names.push_back("C6_inter13_");
      names.push_back("C6_inter14_");
      for (unsigned int i = 0; i < names.size(); ++i) {
        for (set<IJ>::iterator it = IJs.begin(); it != IJs.end(); ++it) {
          stringstream ss;
          ss << names[i] << it->i() << "-" << it->j();
          header.push_back(ss.str());
        }
      }
      for(unsigned int i = 0; i < header.size(); ++i) {
        if (i == 0) {
          cout << "#" << setw(19) << header[i];
        } else {
          cout << setw(20) << header[i];
        }
      }
      cout << endl;
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_tot[*it] * pow(sigmas_tot[*it], 6);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_totinter[*it] * pow(sigmas_totinter[*it], 6);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_totintra[*it] * pow(sigmas_totintra[*it], 6);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra12[*it] * pow(sigmas_intra12[*it], 6);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra13[*it] * pow(sigmas_intra13[*it], 6);
      }
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        cout << setw(20) << 4 * epsilon_intra14[*it] * pow(sigmas_intra14[*it], 6);
      }
      cout << endl;
      
      // calculate the LJ potentials using the C12 and C6 values above...
      map<IJ, LJpot> ftotLJ;
      map<IJ, LJpot> ftotinterLJ;
      map<IJ, LJpot> ftotintraLJ;
      map<IJ, LJpot> fintra12LJ;
      map<IJ, LJpot> fintra13LJ;
      map<IJ, LJpot> fintra14LJ;
      for(set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        double c6 = 4 * epsilon_tot[*it] * pow(sigmas_tot[*it], 6);
        double c12 = 4 * epsilon_tot[*it] * pow(sigmas_tot[*it], 12);
        ftotLJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
        c6 = 4 * epsilon_totinter[*it] * pow(sigmas_totinter[*it], 6);
        c12 = 4 * epsilon_totinter[*it] * pow(sigmas_totinter[*it], 12);
        ftotinterLJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
        c6 = 4 * epsilon_totintra[*it] * pow(sigmas_totintra[*it], 6);
        c12 = 4 * epsilon_totintra[*it] * pow(sigmas_totintra[*it], 12);
        ftotintraLJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
        c6 = 4 * epsilon_intra12[*it] * pow(sigmas_intra12[*it], 6);
        c12 = 4 * epsilon_intra12[*it] * pow(sigmas_intra12[*it], 12);
        fintra12LJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
        c6 = 4 * epsilon_intra13[*it] * pow(sigmas_intra13[*it], 6);
        c12 = 4 * epsilon_intra13[*it] * pow(sigmas_intra13[*it], 12);
        fintra13LJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
        c6 = 4 * epsilon_intra14[*it] * pow(sigmas_intra14[*it], 6);
        c12 = 4 * epsilon_intra14[*it] * pow(sigmas_intra14[*it], 12);
        fintra14LJ.insert(pair<IJ, LJpot>(*it, LJpot(distmin, distmax, distgrid, c12, c6)));
      }
      ofstream fit("test.out");
      printPot(fit, ftotLJ, ftotinterLJ, ftotintraLJ, fintra12LJ, fintra13LJ, fintra14LJ);
      fit.close();
    }
    
    // fitted potentials requested?
    if (printfit) {
      // the least-square fitted potentials
      map<IJ, LJpot> ftotLJ;
      map<IJ, LJpot> ftotinterLJ;
      map<IJ, LJpot> ftotintraLJ;
      map<IJ, LJpot> fintra12LJ;
      map<IJ, LJpot> fintra13LJ;
      map<IJ, LJpot> fintra14LJ;
      LJpot fitpot;
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        fitpot = leastSquareLJ(totLJ[*it]);
        ftotLJ.insert(pair<IJ, LJpot > (*it, fitpot));
        fitpot = leastSquareLJ(totinterLJ[*it]);
        ftotinterLJ.insert(pair<IJ, LJpot > (*it, fitpot));
        fitpot = leastSquareLJ(totintraLJ[*it]);
        ftotintraLJ.insert(pair<IJ, LJpot > (*it, fitpot));
        fitpot = leastSquareLJ(intra12LJ[*it]);
        fintra12LJ.insert(pair<IJ, LJpot > (*it, fitpot));
        fitpot = leastSquareLJ(intra13LJ[*it]);
        fintra13LJ.insert(pair<IJ, LJpot > (*it, fitpot));
        fitpot = leastSquareLJ(intra14LJ[*it]);
        fintra14LJ.insert(pair<IJ, LJpot > (*it, fitpot));
      }
      ofstream fit(fitout.c_str());
      printPot(fit, ftotLJ, ftotinterLJ, ftotintraLJ, fintra12LJ, fintra13LJ, fintra14LJ);
      fit.close();
    }

    // print the distribution, if requested
    if (printBeadBeadDist) {
      ofstream dout(BeadBeadDistOut.c_str());
      double dgrid = (bondlength_max - bondlength_min) / 200;
      dout.precision(9);
      dout << "#" << setw(19) << "r / nm";
      for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
        stringstream ss;
        ss << it->i() << "-" << it->j();
        dout << scientific << setw(20) << ss.str();
      }
      dout << endl;
      for (int i = 0; i < 200; ++i) {
        double r = (i + 0.5) * dgrid + bondlength_min;
        dout << scientific << setw(20) << r;
        for (set<IJ>::const_iterator it = IJs.begin(); it != IJs.end(); ++it) {
          dout << scientific << setw(20) << beadbeadDist[*it][i];
        }
        dout << endl << "#" << endl;
      }
      dout.close();
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
    max =ljp.max;
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
  
  void LJpot::print(ostream& os) {
    os.precision(9);
    if(C12 > 0.0 && C6 > 0.0) {
      os << scientific << "# C12 = " << C12 << " , C6 = " << C6 << endl;
    }
    for (unsigned int i = 0; i < lj.size(); ++i) {
      if (count[i] == 0) {
        os << scientific << setw(20) << r(i) << setw(20) << 0.0 << endl;
      } else {
        os << scientific << setw(20) << r(i) << setw(20) << lj[i] / count[i] << endl;
      }
    }
    os << endl;
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
  
  int LJpot::size() {
    return lj.size();
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
    for(int i = 1; i < get_grid(); ++i) {
      if(pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double r_min = r(i_max);
    double pot_min = pot(i_max);
    for(int i = i_max + 1; i < get_grid(); ++i) {
      if(pot(i) < pot_min) {
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
    for(int i = 1; i < get_grid(); ++i) {
      if(pot(i) > pot_max) {
        i_max = i;
        pot_max = pot(i);
      }
    }
    // now find the minimum
    double r_min = r(i_max);
    double pot_min = pot(i_max);
    for(int i = i_max + 1; i < get_grid(); ++i) {
      if(pot(i) < pot_min) {
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
    for(; it != ij.end(); ++it) {
      totLJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
      totinterLJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
      totintraLJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
      intra12LJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
      intra13LJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
      intra14LJ.insert(pair<IJ, LJpot>(*it, LJpot(min, max, grid)));
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

  void bead::addLJtot(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    totLJ[ij].add(r, lj);
  }

  void bead::addLJtotinter(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    totinterLJ[ij].add(r, lj);
  }

  void bead::addLJtotintra(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    totintraLJ[ij].add(r, lj);
  }

  void bead::addLJintra12(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    intra12LJ[ij].add(r, lj);
  }

  void bead::addLJintra13(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    intra13LJ[ij].add(r, lj);
  }

  void bead::addLJintra14(const IJ &ij, double r, const double &lj, double min, double max, int grid) {
    intra14LJ[ij].add(r, lj);
  }
  
  bool operator<(const IJ &ij1, const IJ &ij2) {
    if(ij1.i() < ij2.i() || (ij1.i() == ij2.i() && ij1.j() < ij2.j())) {
      return true;
    }
    return false;
  }
  
  LJpot leastSquareLJ(LJpot &ljp) {
    // for the fitting we only consider the values V_LJ(r) which are not bigger
    // than -min(LJ(r))
    double LJ_min = 0.0;
    int i_LJ_min = 0;
    for(int i = 0 ;i < ljp.size(); ++i) {
      if(ljp.pot(i) < LJ_min) {
        LJ_min = ljp.pot(i);
        i_LJ_min = i;
      }
    }
    int n = ljp.size();
    bool b = false;
    // find the smallest r with LJ(r) < - min(LJ)
    for(int i = 0 ;i < ljp.size(); ++i) {
      if(ljp.pot(i) < - LJ_min && b) {
        break;
      }
      if(ljp.pot(i) > - LJ_min) {
        b = true;
      }
      n--;
    }
    // just in case it is never more positive than more negative...
    if(n == 0) {
      n = ljp.size();
    }

    gsl_matrix * X = gsl_matrix_alloc(n, 2);
    gsl_vector *y = gsl_vector_alloc(n);
    for (int l = i_LJ_min; l < ljp.size(); ++l) {
      double r3 = ljp.r(l) * ljp.r(l) * ljp.r(l);
      double r6 = r3 * r3;
      double r12 = r6 * r6;
      gsl_matrix_set(X, l - i_LJ_min, 0, 1/r12);
      gsl_matrix_set(X, l - i_LJ_min, 1, -1/r6);
      gsl_vector_set(y, l - i_LJ_min, ljp.pot(l));
    }
    gsl_vector *c = gsl_vector_alloc(2);
    gsl_matrix * cov = gsl_matrix_alloc (2, 2);
    double chisq;
    gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc (n, 2);
    // let's do the least-square fit
    gsl_multifit_linear(X, y, c, cov, &chisq, work);
    // free what is not needed any more
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_multifit_linear_free(work);
    
    // put the result in a LJpot
    LJpot lj(ljp.get_min(), ljp.get_max(), ljp.get_grid(), gsl_vector_get(c, 0), gsl_vector_get(c, 1));
    // and free the rest
    gsl_vector_free(c);
    
    return lj;
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
    
  }
  
}
