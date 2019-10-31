/**
 * @file addvirt_top.cc
 * create virtual atoms in a topology
 */

/**
 * @page programs Program Documentation
 *
 * @anchor addvirt_top
 * @section addvirt_top Add virtual atoms to a topology
 * @author @ref co
 * @date 29-10-2019
 *
 *
 * <hr>
 */


#include <cassert>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/VirtualAtom.h"
#include "../src/utils/Neighbours.h"
using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "atomspec" << "hydrogens" << "dish" << "disc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology files>\n";
  usage += "\t[@atomspec <virtual atom specified>]\n";
  usage += "\t[@hydrogens  <iac of united atoms to add H's for>]\n";
  usage += "\t[@dish      <dish> (default: 0.1)]\n";
  usage += "\t[@disc      <disc> (default: 0.153)]\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());    

    gcore::VirtualAtoms vas = sys.vas();
    int numOrigVas = vas.numVirtualAtoms();

    // read dish and disc
    double dish = args.getValue<double>("dish", false, 0.1);
    double disc = args.getValue<double>("disc", false, 0.153);
    if(vas.numVirtualAtoms()){
      if(args.count("dish") > 0 && dish!=vas.dish()){ 
        std::cerr << "WARNING: specified value of dish does not match current value in topology\n"
                << "         Specified: " << dish
                << " in topology: " << vas.dish()  
                << "\n         Taking specified value\n";
      } else {
        dish = vas.dish();
      }
      if(args.count("disc") > 0 && disc!=vas.disc()){ 
        std::cerr << "WARNING: specified value of disc does not match current value in topology\n"
                << "         Specified: " << disc
                << " in topology: " << vas.disc()  
                << "\n         Taking specified value\n";
      } else {
        disc = vas.disc();
      }
    }
    vas.setDis(dish, disc);

    // read atomspecifier
    utils::AtomSpecifier as(sys);
    for(Arguments::const_iterator iter=args.lower_bound("atomspec"), 
          to=args.upper_bound("atomspec"); iter!=to;++iter) {

      as.addSpecifier(iter->second);
    } 
    vas.addVirtualAtom(as);
    

    // any explicit hydrogens
    std::set<int> uas;
    for(Arguments::const_iterator iter=args.lower_bound("hydrogens"), 
          to=args.upper_bound("hydrogens"); iter!=to;++iter) {

      int iac;
      istringstream is(iter->second);
      if(!(is >> iac))
	throw gromos::Exception("addvirt_top", 
				"failed to read an integer from @hydrogens"
				+ iter->second);

      uas.insert(iac -1);
    }
    // loop over all atoms in the system
    int totNumAtoms = 0;
    for(int m = 0; m < sys.numMolecules(); m++){
      for(int a = 0; a < sys.mol(m).numAtoms(); a++){
        if (uas.count(sys.mol(m).topology().atom(a).iac())) {
          utils::Neighbours n(sys, m, a);
          std::vector<int> n1, n2;
          switch(n.size()){
            case 1:
              n1.push_back(totNumAtoms + a); 
              n1.push_back(totNumAtoms + n[0]);
              vas.addVirtualAtom(sys, n1, 5, dish, disc);
              break;
            case 2:
              n1.push_back(totNumAtoms + a);
              n1.push_back(totNumAtoms + n[0]);
              n1.push_back(totNumAtoms + n[1]);
              vas.addVirtualAtom(sys, n1, 4, dish, disc);
              n2.push_back(totNumAtoms + a);
              n2.push_back(totNumAtoms + n[1]);
              n2.push_back(totNumAtoms + n[0]);
              vas.addVirtualAtom(sys, n2, 4, dish, disc);
              break;
            case 3:
              n1.push_back(totNumAtoms + a);
              n1.push_back(totNumAtoms + n[0]);
              n1.push_back(totNumAtoms + n[1]);
              n1.push_back(totNumAtoms + n[2]);
              vas.addVirtualAtom(sys, n1, 1, dish, disc);
              break;
            default:
              ostringstream error;
              error << "don't know how to add a virtual atom to atom " 
                    << a+1 << " in molecule " << m+1 
                    << " (option @hydrogens)";
              throw gromos::Exception("addvirt_top", error.str());
          }
        } 
      }
      totNumAtoms += sys.mol(m).numAtoms();
    }

    // put the VirtualAtoms back in the system
    sys.addVirtualAtoms(vas);

    OutTopology ot(cout);
  
    ostringstream title;
    title << it.title()
          << "addvirt_top added " << vas.numVirtualAtoms() - numOrigVas << " virtual atoms"; 
    ot.setTitle(title.str());
    ot.write(sys,it.forceField());
  } catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




