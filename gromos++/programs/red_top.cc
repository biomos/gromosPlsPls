/**
 * @file red_top.cc
 * Reduce a molecular topology to a subset of atoms
 */

/**
 * @page programs Program Documentation
 *
 * @anchor red_top
 * @section red_top Reduca a molecular topology to a subset of atoms
 * @author @ref co
 * @date 7-6-07
 *
 * For large molecular complexes, one would sometimes like to consider only a 
 * part of the many atoms, thereby reducing the computational effort required
 * by a simulation.
 * Programs @ref tstrip and @ref filter can filter an atomic coordinate file 
 * (see sections V-4.1 and V-4.2, respectively). Accordingly, program red_top 
 * can cut out parts of a molecular topology.
 *
 * The user has to list atoms of the molecular topology to be reduced. All 
 * atoms, exclusions, bonds, bond angles etc. that involve atoms that are not
 * in this list are removed. Note that to cut out all atoms within a sphere
 * around a part of the system, one could first generate a list of the 
 * corresponding atoms using the program @ref pairlist (section V-5.5). 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" in the system to keep&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  red_top
    @topo   ex.top
    @atoms  1:4-34,40
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){
  Argument_List knowns;
  knowns << "topo" << "atoms";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@atoms <atoms in the system to keep>\n";

  try{
    Arguments args(argc, argv, knowns, usage);

    InTopology it(args["topo"]);
    System sys(it.system());
    
    // Parse atom specifiers
    utils::AtomSpecifier as(sys);
    for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
          to=args.upper_bound("atoms"); iter!=to;++iter) {
      as.addSpecifier(iter->second);
    }
    vector<string> as_str = as.toString();

    // flag all atoms that are not in the list with a negative iac
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a<sys.mol(m).numAtoms(); a++){
	if(as.findAtom(m,a)==-1) sys.mol(m).topology().atom(a).setIac(-1);
      }
    }

    // create a linear topology
    gcore::LinearTopology lt(sys);

    // remove all flagged atoms
    lt.removeAtoms();
    
    // calculate the new 1,4 interactions
    lt.get14s();

    // delete those 1,4 interactions if they are in LJ exceptions
    for(set<LJException>::const_iterator it = lt.ljexceptions().begin();
            it != lt.ljexceptions().end(); it++) { 
      int at1 = (*it)[0];
      int at2 = (*it)[1];
      int n14 = lt.atoms()[at1].exclusion14().size();
      for(int i = 0; i < n14; i++) {
        if(lt.atoms()[at1].exclusion14().atom(i) == at2) {
          lt.atoms()[at1].exclusion14().erase(at2);
        }
      }
    }

    // and parse the linearized thing back into a topology (which might
    // have a zillion molecules, because bonds have been removed)
    System syo = lt.parse();
    
    // take the old solvent
    syo.addSolvent(sys.sol(0));
    
    // set the temperature and pressure groups
    {
      int a = 0;
      for (int m = 0; m < syo.numMolecules(); ++m) {
        a += syo.mol(m).numAtoms();
        syo.addPressureGroup(a);
        syo.addTemperatureGroup(a);
      }
      // compare the current number with the previous one and warn if it changed for other reasons than deleted molecuels
      if (sys.numMolecules() != sys.numPressureGroups()) {
        if (sys.numPressureGroups() != syo.numPressureGroups()) {
          cerr << "WARNING: The number of pressure groups has changed. manual check recommended.\n";
          cerr << "         Number of pressure greoups: " << sys.numPressureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                     " << syo.numPressureGroups() << " (reduced topology)\n";
        }
      }
      if (sys.numMolecules() != sys.numTemperatureGroups()) {
        if (sys.numTemperatureGroups() != syo.numTemperatureGroups()) {
          cerr << "WARNING: The number of temperature groups has changed. manual check recommended.\n";
          cerr << "         Number of temperature greoups: " << sys.numTemperatureGroups() << " (" << args["topo"] << ")" << endl
                  << "                                        " << syo.numTemperatureGroups() << " (reduced topology)\n";
        }
      }
    }
    
    // and write out the new topology
    OutTopology ot(cout);
    ostringstream os;
    os << "Reduced topology based on " << args["topo"] << endl;
    os << "using atoms ";
    for(unsigned int i=0; i< as_str.size(); i++)
      os << as_str[i] << " ";
    
    ot.setTitle(os.str());
    
    ot.write(syo, it.forceField());
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}



