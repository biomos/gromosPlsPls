/**
 * @file atominfo.cc
 * list characteristics of atoms and convert AtomSpecifier to gromos 
 * numbering and vv
 */

/**
 * @page programs Program Documentation
 *
 * @anchor atominfo
 * @section atominfo list atom characteristics
 * @author @ref co
 * @date 16. 3. 2005
 *
 * Lists names, atomtypes, mass charge etc. for individual atoms in a 
 * topology. Takes gromosnumbers or AtomSpecifier for input and
 * creates a table for all specified atoms.
 *
 * arguments:
 * - topo      <topology>
 * - gromosnum <numbers of atoms according to numbering in topology>
 * - atomspec  <atoms specified through the specifier>
 *
 * Example:
 * @verbatim
 * atominfo
 *   @topo ex.top
 *   @gromosnum 43
 *   @atomspec 1:CA
 * @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace gio;
using namespace args;

using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "gromosnum", "atomspec"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo      <topology>\n";
  usage += "\t@gromosnum <gromos atom number>\n";
  usage += "\t@atomspec  <atomspecifier>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    utils::AtomSpecifier as(sys);
    
    
    // read gromos numbers
    Arguments::const_iterator iter=args.lower_bound("gromosnum"),
      to=args.upper_bound("gromosnum");
    for(;iter!=to;++iter){
      int grom;
      istringstream is(iter->second);
      if(!(is >> grom))
	throw gromos::Exception("atominfo", 
				"failed to read an integer from input"
				+ iter->second);
      as.addGromosAtom(grom-1);
      
    }
    
    // and atomspecifiers
    iter=args.lower_bound("atomspec");
    to=args.upper_bound("atomspec");
    for(;iter!=to; ++iter){
      as.addSpecifier(iter->second);
    }


    if(as.size()){
      
      cout << setw(10) << "Atom"
	   << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	   << setw(10) << "Residue"
	   << setw(10) << "Atom"
	   << setw(12) << "Integer"
	   << setw(10) << "Charge" << endl;
      cout << setw(10) << "Specifier"
	   << setw(10) << "number"
	   << setw(10) << "number"
	   << setw(10) << "name"
	   << setw(10) << "name"
	   << setw(12) << "Atom Code"
	   << endl;
    }
    
    for(int i=0; i < as.size(); ++i){

      //determine the residue name
      cout << setw(10) << as.toString(i)
	   << setw(10) << as.gromosAtom(i)+1
	   << setw(10) << as.resnum(i)+1
	   << setw(10) << as.resname(i)
	   << setw(10) << as.name(i)
	   << setw(12) << as.iac(i)+1
	   << setw(10) << as.charge(i)
	   << endl;
    }
    
      
	
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
