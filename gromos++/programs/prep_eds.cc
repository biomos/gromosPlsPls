/**
 * @file prep_eds.cc
 * Generate molecular topology and perturbation topology files for EDS
 */

/**
 * @page programs Program Documentation
 *
 * @anchor prep_eds
 * @section prep_eds Generate molecular topology and perturbation topology files for EDS
 * @author @ref sr
 * @date 19-1-10
 *
 * The molecular topology file for EDS is generated from N 'normal' topologies,
 * where N is the number of states. A state is in EDS a molecule, e.g. a ligand.
 * In the EDS-topology, all states or molecules, respectively, are combined and
 * excluded from another. The resulting molecular topology file is written out
 * to a file called com_eds.top.
 *
 * In the EDS perturbation topology, a molecule is 'visible' in one state and
 * in all other states it consists of dummy atoms. For this the MPERTATOM block
 * is used. The resulting perturbation topology file is written out to a file
 * called pert_eds.ptp.
 *
 * The argument \@inG96 converts GROMOS96 topologies to the current formats. On
 * the other hand \@outG96 converts topologys in the current format to the
 * GROMOS96 format.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology files&gt; </td></tr>
 * <tr><td> \@numstat</td><td>&lt;number of EDS states&gt; </td></tr>
 * <tr><td> \@param</td><td>&lt;index number of molecular topology file to take parameters from&gt; </td></tr>
 * <tr><td> \@solv</td><td>&lt;index number of molecular topology file to take solvent from&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  prep_eds
    @topo      brbn.top clbn.top mebn.top
    @numstat   3
    @param     1
    @solv      1
 @endverbatim
 *
 * <hr>
 */


#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/OutPtTopology.h"
#include "../src/gcore/PtTopology.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "numstat" << "param" << "solv";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology files>\n";
  usage += "\t@numstat <number of EDS states>\n";
  usage += "\t@param <index number of molecular topology file to take parameters from>\n";
  usage += "\t@solv  <index number of molecular topology file to take solvent from>\n";

  try {
    Arguments args(argc, argv, knowns, usage);
    // set some values
    int numstat = 1;
    if (args.count("numstat") > 0) numstat = atoi(args["numstat"].c_str());
    if (args.count("topo") <= 0)
      throw gromos::Exception("prep_eds", "needs at least one topology\n" + usage);
    if (args.count("topo") > numstat)
      throw gromos::Exception("prep_eds", "too many topologies\n" + usage);
    if (args.count("topo") < numstat)
      throw gromos::Exception("prep_eds", "too few topologies\n" + usage);
    int parnum=1;
    if(args.count("param")>0) parnum=atoi(args["param"].c_str());
    int solnum=1;
    if(args.count("solv")>0)  solnum=atoi(args["solv"].c_str());

    System sys;

    //----------------------------
    //First, the molecular topology
    // ---------------------------
    string paramname, toponame;

    ostringstream title;
    title << "prep_eds: Combined EDS topology using:\n";

    ostringstream title_ptp;
    title_ptp << "prep_eds: EDS perturbation topology using: \n";

    int totNumAt = 0, topnum = 1;
    for (Arguments::const_iterator iter = args.lower_bound("topo"),
            to = args.upper_bound("topo"); iter != to; ++iter) {

      toponame = iter->second;

      // read topology
      InTopology it(toponame);
      // Directly add pressure and temperature groups
      for (int j = 0; j < it.system().numTemperatureGroups(); j++) {
        sys.addTemperatureGroup(it.system().temperatureGroup(j) + totNumAt);
      }

      for (int j = 0; j < it.system().numPressureGroups(); j++) {
        sys.addPressureGroup(it.system().pressureGroup(j) + totNumAt);
      }

      // Add molecules and count new number of atoms in sys
      for (int j = 0; j < it.system().numMolecules(); j++) {
        sys.addMolecule(it.system().mol(j));
        totNumAt += it.system().mol(j).numAtoms();
      } // molecules

      if (solnum == topnum)
        sys.addSolvent(it.system().sol(0));
      if (parnum == topnum)
        paramname = toponame;

      title << toponame << endl;
      title_ptp << toponame << endl;
    }

    // Add the additional exclusions to the atoms
    int start_atom = 0;
    int num_atom = start_atom;
    for (int j = 0; j < (sys.numMolecules() - 1); j++) {
      start_atom = sys.mol(j).numAtoms();
      num_atom = start_atom;
      for (int i = 0; i < sys.mol(j).numAtoms(); i++) {
        for (int l = j + 1; l < sys.numMolecules(); l++) {
          for (int k = 0; k < sys.mol(l).topology().numAtoms(); k++) {
            sys.mol(j).topology().atom(i).exclusion().insert(num_atom);
            num_atom++;
          } // atom k of mol l
        } // mol l
        num_atom = start_atom;
      } // atom i of mol j
    } // mol j

    InTopology it(paramname);
    title << "Parameters from " << parnum
            << ", solvent from " << solnum;

    ofstream topo("com_eds.top");
    OutTopology ot(topo);

    ot.setTitle(title.str());
    ot.write(sys, it.forceField());

    //-------------------------------
    // Now, the perturbation topology
    //-------------------------------

    PtTopology pt;
    pt.setSize(totNumAt, numstat);
    for (int p = 0; p < numstat; p++) {
      int atm = 0;
      std::stringstream statename;
      statename << "state" << p;
      pt.setPertName(p, statename.str());
      for (int i = 0; i < sys.numMolecules(); i++) {
        if (i == p) {
          for (int k = 0; k < sys.mol(i).numAtoms(); k++) {
            pt.setIac(atm, p, sys.mol(i).topology().atom(k).iac());
            pt.setAtomName(atm, sys.mol(i).topology().atom(k).name());
            pt.setCharge(atm, p, sys.mol(i).topology().atom(k).charge());
            pt.setAtomNum(atm, atm);
            atm++;
          } // atoms
        } else {
          for (int k = 0; k < sys.mol(i).numAtoms(); k++) {
            pt.setIac(atm, p, 21);
            pt.setAtomName(atm, sys.mol(i).topology().atom(k).name());
            pt.setCharge(atm, p, 0.0);
            pt.setAtomNum(atm, atm);
            atm++;
          } // atoms
        }
      } // molecules
    } // number of states

    ofstream ptp("pert_eds.ptp");
    OutPtTopology op(ptp);

    op.setTitle(title_ptp.str());
    op.write_multiple(pt);

    
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




