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
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "atoms"};
  int nknowns = 2;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@atoms <atoms to keep>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

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

    // and parse the linearized thing back into a topology (which might
    // have a zillion molecules, because bonds have been removed)
    System syo = lt.parse();
    
    // take the old solvent
    syo.addSolvent(sys.sol(0));
    
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



