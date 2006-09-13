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
 * Internally the gromos preparation and analysis tools determine which atoms
 * belong to one molecule based on bonds specified in the topology. These
 * programs can make use of the convenient @ref AtomSpecifier "atomspecifier"
 * to select atoms, molecular @ref PropertySpecifier "properties" etc. for 
 * efficiency reasons, both MD-engines, promd and md, number all atoms in the
 * molecular system sequentially. Program atominfo can read both 
 * @ref AtomSpecifier "atom specifiers" and sequential numbers (gromos-numbers)
 * and will list the properties of the selected atoms.
 *
 * The atomlist can be sorted, according to the following priority: solute
 * molecule &lt; virtual atom &lt; solvent molecule. All programs that make
 * use of atom specifiers can also read in a file containing the output of
 * atominfo, by specifying "file(&lt;atominfo output file&gt;)". This allows
 * the user to store complicated selections in a file for future use.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;topology&gt; </td></tr>
 * <tr><td> \@gromosnum</td><td>&lt;gromos atom number&gt; </td></tr>
 * <tr><td> \@atomspec</td><td>&lt;atomspecifier&gt; </td></tr>
 * <tr><td> [\@sort</td><td>(sort the atoms)] </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo       ex.top
     @gromosnum  43
     @atomspec   1:CA
     @sort
   @endverbatim

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
#include "../src/utils/VirtualAtom.h"

using namespace gcore;
using namespace gio;
using namespace args;

using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "gromosnum", "atomspec", "sort"};
  int nknowns = 4;

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <topology>\n";
  usage += "\t@gromosnum <gromos atom number>\n";
  usage += "\t@atomspec  <atomspecifier>\n";
  usage += "\t[@sort     (sort the atoms)]\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());


    // add one set of solvent coordinates to the solvent, so that we can 
    // access them
    sys.sol(0).setNumPos(sys.sol(0).topology().numAtoms());
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

      if (args.count("sort") != -1)
	as.sort(); 

      cout << "TITLE\n\tatominfo\n";
      cout << "\t" << as.toString()[0];
      cout << "\nEND\n";
      cout << "ATOMS\n";
      
      cout << "#"
	   << setw(12) << "Atom"
	   << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	   << setw(10) << "Residue"
	   << setw(10) << "Atom"
	   << setw(12) << "Integer"
	   << setw(10) << "Charge" << endl;
      cout << "#"
	   << setw(12) << "Specifier"
	   << setw(10) << "number"
	   << setw(10) << "number"
	   << setw(10) << "name"
	   << setw(10) << "name"
	   << setw(12) << "Atom Code"
	   << endl;
    }
    
    for(int i=0; i < as.size(); ++i){
      if(as.atom()[i]->type()==utils::spec_virtual){
	utils::AtomSpecifier conf=as.atom()[i]->conf();
	cout << "----------------------------------------"
	     << "--------------------------------\n"
	     << "virtual atom, ";
	switch(as.atom()[i]->virtualType()){
	  case utils::VirtualAtom::normal: 
	    cout << "explicit atom:\n";
	    break;
	  case utils::VirtualAtom::CH1:
	    cout << "aliphatic CH1 group, based on " << conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::aromatic:
	    cout << "aromatic CH1 group, based on "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::CH2:
	    cout << "non-stereospecific aliphatic CH2 group (pseudo atom),"
		 << " based on "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::stereo_CH2:
	    cout << "stereospecific aliphatic CH2, based on "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::stereo_CH3:
	    cout << "single CH3 groups (psuedo atom), based on "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::CH3:
	    cout << "non-stereospecific CH3 groups (isopropyl; pseudo atom), "
		 << "based on "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::ring:
	    cout << "aromatic flipping ring (pseudo atom), positioned at:\n";
	    break;
	  case utils::VirtualAtom::NH2:
	    cout << "non-stereospecific NH2 group (pseudo atom), "
		 << "based on "<< conf.size() << " atoms:\n";
	    break;
	  case utils::VirtualAtom::COM:
	    cout << "center of mass for "<< conf.size() 
		 << " atoms:\n";
	    break;
	  case utils::VirtualAtom::COG:
	    cout << "center of geometry for "<< conf.size() 
		 << " atoms:\n";
	    break;
	    
	}
	
	for(int j=0; j< conf.size(); ++j){
	  cout << setw(13) << conf.toString(j)
	       << setw(10) << conf.gromosAtom(j)+1
	       << setw(10) << conf.resnum(j)+1
	       << setw(10) << conf.resname(j)
	       << setw(10) << conf.name(j)
	       << setw(12) << conf.iac(j)+1
	       << setw(10) << conf.charge(j)
	       << endl;
	}
	cout << "----------------------------------------"
	     << "--------------------------------\n";
	
      }
      else{
	
	// print out normal atoms
	cout << setw(13) << as.toString(i)
	     << setw(10) << as.gromosAtom(i)+1
	     << setw(10) << as.resnum(i)+1
	     << setw(10) << as.resname(i)
	     << setw(10) << as.name(i)
	     << setw(12) << as.iac(i)+1
	     << setw(10) << as.charge(i)
	     << endl;
      }
    }
    
    cout << "END\n";
	
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
