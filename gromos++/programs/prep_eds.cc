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
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/OutPtTopology.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
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
    int numstat = args.getValue<int>("numstat", false, 1);
    if (args.count("topo") <= 0)
      throw gromos::Exception("prep_eds", "needs at least one topology\n" + usage);
    if (args.count("topo") > numstat)
      throw gromos::Exception("prep_eds", "too many topologies\n" + usage);
    if (args.count("topo") < numstat)
      throw gromos::Exception("prep_eds", "too few topologies\n" + usage);
    int parnum = args.getValue<int>("param", false, 1);
    int solnum = args.getValue<int>("solv", false, 1);

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
    vector<int> size_topo(numstat,0), last_atom(numstat,0);
    int i_topo = 0;
    for (Arguments::const_iterator iter = args.lower_bound("topo"),
            to = args.upper_bound("topo"); iter != to; ++iter, ++i_topo) {

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

      size_topo[i_topo] = it.system().numMolecules();
      for (int i = 0; i < it.system().numMolecules(); ++i) {
        last_atom[i_topo] += it.system().mol(i).numAtoms();
      }
    }
    int last_mol = 0;
    for (int i = 0; i < (numstat-1); ++i) {
      last_mol += size_topo[i];
    }

    for (int i = 1; i < numstat; ++i) {
      last_atom[i] += last_atom[i-1];
    }

    // Add the additional exclusions to the atoms
    int start_atom = last_atom[0], end_atom = last_atom[numstat-1];
    int n = 0, adjust_atom = 0;
    for (int j = 0; j < last_mol; j++) {
      for (int i = 0; i < sys.mol(j).numAtoms(); i++) {
        for (int l = start_atom; l < end_atom; l++) {
          if (j == 0) {
            sys.mol(j).topology().atom(i).exclusion().insert(l);
          } else {
            sys.mol(j).topology().atom(i).exclusion().insert(l - adjust_atom);
          }
        } // mol l
      } // atom i of mol j
      adjust_atom += sys.mol(j).numAtoms();
      if ((j+1) == size_topo[n]) {
        ++n;
        start_atom = last_atom[n];
      }
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
    int topo_mol_end = 0, topo_mol_start = 0;
    for (int p = 0; p < numstat; p++) {
      int atm = 0;
      std::stringstream statename;
      statename << "state" << p;
      pt.setPertName(p, statename.str());
      topo_mol_start = topo_mol_end;
      topo_mol_end += size_topo[p];
      for (int i = 0; i < sys.numMolecules(); i++) {
        if (i >= topo_mol_start && i < topo_mol_end) {
          for (int k = 0; k < sys.mol(i).numAtoms(); k++) {
            pt.setIac(atm, p, sys.mol(i).topology().atom(k).iac());
            pt.setAtomName(atm, sys.mol(i).topology().atom(k).name());
            pt.setCharge(atm, p, sys.mol(i).topology().atom(k).charge());
            pt.setAtomNum(atm, atm);
            pt.setAlphaLJ(atm, 1.0);
            pt.setAlphaCRF(atm, 1.0);
            atm++;
          } // atoms
        } else {
          for (int k = 0; k < sys.mol(i).numAtoms(); k++) {
            pt.setIac(atm, p, 21);
            pt.setAtomName(atm, sys.mol(i).topology().atom(k).name());
            pt.setCharge(atm, p, 0.0);
            pt.setAtomNum(atm, atm);
            pt.setAlphaLJ(atm, 1.0);
            pt.setAlphaCRF(atm, 1.0);
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




