/**
 * @file hvap
 * Calculates the heat of vaporistaion from a coordinate trajectory
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor hvap
 * @section hval calculates the heat of vaporistaion from a coordinate trajectory
 * @author @ref ae
 * @date 07.07.2011
 *
 * Program hvap ...
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
#include <iomanip>
#include <iostream>
#include <vector>

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_vector.h>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gmath/Physics.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"

using namespace args;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace std;

double LJ(int iac1, int iac2, double r2, GromosForceField &gff);
double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut);
void putIntoBox(Vec &v, Box &box);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pbc" << "RD" << "trc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo     <molecular topology file>\n";
  usage += "\t[@pbc     <boundary type (read from GENBOX block if not specified)> [<gather method>]>]\n";
  usage += "\t@rf       <cut off radiuse> <epsilon> <kappa>\n";
  usage += "\t@trc      <simulation trajectory or coordinate file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read reaction fiel parameters
    double eps, kap, cut;
    if (args.count("rf") != 3) {
      //throw gromos::Exception(argv[0], "too few arguments for the reaction-field parameters (@rf)");
    } else {
      Arguments::const_iterator it = args.lower_bound("rf");
      stringstream ss;
      ss << it->second;
      ss >> cut;
      ++it;
      ss << it->second;
      ss >> eps;
      ++it;
      ss << it->second;
      ss >> kap;
    }

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff = it.forceField();

    if (args.count("trc") < 1) {
      throw gromos::Exception(argv[0], "no coordinate or trajectory file specified (@trc)");
    }

    vector<double> hvaps;

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
      ic.select("ALL");
      int numMol = sys.numMolecules();
      while (!ic.eof()) {

        double hvap = 0.0;

        ic >> sys;

        // some numnbers needed in case tere is solvent
        int totAtSolv = sys.sol(0).numAtoms();
        int numAtSolv = sys.sol(0).topology().numAtoms();
        int numSolvMol = totAtSolv / numAtSolv;

        // loop over the intermolecular atom pairs
        for (int m1 = 0; m1 < numMol; ++m1) {
          int numAt1 = sys.mol(m1).numAtoms();
          for (int a1 = 0; a1 < numAt1; ++a1) {
            int iac1 = sys.mol(m1).topology().atom(a1).iac();
            Vec pos1 = sys.mol(m1).pos(a1);
            double charge1 = sys.mol(m1).topology().atom(a1).charge();
            // the solute-solute unteractions
            for (int m2 = m1 + 1; m2 < numMol; ++m2) {
              int numAt2 = sys.mol(m2).numAtoms();
              for (int a2 = 0; a2 < numAt2; ++a2) {
                int iac2 = sys.mol(m2).topology().atom(a2).iac();
                Vec pos2 = sys.mol(m2).pos(a2);
                double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
                double charge2 = sys.mol(m2).topology().atom(a2).charge();
                hvap += LJ(iac1, iac2, r2, gff);
                hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
              }
            }
            // the solute-solvent interactions
            for (int s = 0; s < totAtSolv; ++s) {
              int iac2 = sys.sol(0).topology().atom(s % numAtSolv).iac();
              Vec pos2 = sys.sol(0).pos(s % numAtSolv);
              double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
              double charge2 = sys.sol(0).topology().atom(s % numAtSolv).charge();
              hvap += LJ(iac1, iac2, r2, gff);
              hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
            }
          }
        }
        // and the solvent-solvent interactions
        for (int s = 0; s < totAtSolv; s += numAtSolv) {
          for (int s1 = s; s1 < s + numAtSolv; ++s1) {
            for (int s2 = s + numAtSolv; s2 < totAtSolv; ++s2) {
              int iac1 = sys.sol(0).topology().atom(s1 % numAtSolv).iac();
              int iac2 = sys.sol(0).topology().atom(s2 % numAtSolv).iac();
              Vec pos1 = sys.sol(0).pos(s1 % numAtSolv);
              Vec pos2 = sys.sol(0).pos(s2 % numAtSolv);
              double r2 = (pos1 - pbc->nearestImage(pos1, pos2, sys.box())).abs2();
              double charge1 = sys.sol(0).topology().atom(s1 % numAtSolv).charge();
              double charge2 = sys.sol(0).topology().atom(s2 % numAtSolv).charge();
              hvap += LJ(iac1, iac2, r2, gff);
              hvap += Coulomb(charge1, charge2, r2, eps, kap, cut);
            }
          }
        }

        // divide by the number of molecules (solute and solvent)
        hvap /= (numMol + numSolvMol);
        hvaps.push_back(hvap);

      } // end of loop over configuration of the current trajectory file
    }

    double hvap = 0.0;
    for (unsigned int i = 0; i < hvaps.size(); ++i) {
      hvap += hvaps[i];
    }
    hvap /= hvaps.size();

    cout << hvap << endl;

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

double LJ(int iac1, int iac2, double r2, GromosForceField &gff) {
  double r6 = r2 * r2 * r2;
  AtomPair ap(iac1, iac2);
  double c12 = gff.ljType(ap).c12();
  double c6 = gff.ljType(ap).c6();
  return (c6 - c12 / r6) / r6;
}

double Coulomb(double c1, double c2, double r2, double eps, double kap, double cut) {
  return -(c1 * c2 * physConst.get_four_pi_eps_i() / r2);
}
