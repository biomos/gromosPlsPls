#include <iostream>
#include <sstream>
#include <set>
#include <map>
#include <vector>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/BbSolute.h"
#include "../src/gcore/BuildingBlock.h"
#include "../src/utils/AtomSpecifier.h"



using namespace gcore;
using namespace gio;
using namespace args;

// we use functionality of maketop
#include "../src/utils/maketop.h"

void linearizeTopology(vector<AtomTopology> &atoms, 
		       vector<Bond> &bonds,
		       vector<Angle> &angle,
		       vector<Improper> &imp,
		       vector<Dihedral> &dih,
		       map<int, int> &resnum,
		       vector<string> &resname,
		       const System &sys, 
		       utils::AtomSpecifier &keep);

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
 
    // linearize the topology
    // flag all atoms that are to be removed with a negative iac
    vector<AtomTopology> atoms;
    vector<string> resNames;
    vector<Bond> bonds;
    vector<Angle> angles;
    vector<Improper> imps;
    vector<Dihedral> dihs;
    map<int, int> resMap;

    linearizeTopology(atoms, bonds, angles, imps, dihs, 
		      resMap, resNames, sys, as);

    // now call the function removeAtom from maketop
    removeAtoms(&atoms, &bonds, &angles, &imps, &dihs, &resMap, &resNames);
 
    // calculate the new 1,4 interactions
    get14s(&atoms, &bonds);
    
    // and parse the linearized thing back into a topology (which might
    // have a zillion molecules, because bonds have been removed)
    System syo = parseTopology(&atoms, &bonds, &angles, &imps, &dihs, 
			       &resNames, &resMap);
    
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





void linearizeTopology(vector<AtomTopology> &atom, 
		       vector<Bond> &bond,
		       vector<Angle> &angle,
		       vector<Improper> &imp,
		       vector<Dihedral> &dih,
		       map<int, int> &resmap,
		       vector<string> &resname,
		       const System &sys,
		       utils::AtomSpecifier &keep)
{
  int lastAtom=0;
  int lastRes=0;

  for(int m=0; m< sys.numMolecules(); m++){
    for(int a=0; a<sys.mol(m).numAtoms(); a++){
      atom.push_back(sys.mol(m).topology().atom(a));
      resmap[a+lastAtom]=sys.mol(m).topology().resNum(a)+lastRes;
      if(keep.findAtom(m,a) == -1) atom[a+lastAtom].setIac(-1);
    }
    for(int i=0; i<sys.mol(m).topology().numRes(); i++)
      resname.push_back(sys.mol(m).topology().resName(i));
      
    BondIterator bi(sys.mol(m).topology());
    for(; bi; ++bi){
      Bond b=bi();
      b[0]+=lastAtom; b[1]+=lastAtom;
      bond.push_back(b);
    }

    AngleIterator ai(sys.mol(m).topology());
    for(; ai; ++ai){
      Angle a=ai();
      a[0]+=lastAtom; a[1]+=lastAtom; a[2]+=lastAtom;
      angle.push_back(a);
    }
    
    DihedralIterator di(sys.mol(m).topology());
    for(; di; ++di){
      Dihedral d=di();
      d[0]+=lastAtom; d[1]+=lastAtom; d[2]+=lastAtom; d[3]+=lastAtom;
      dih.push_back(d);
    }

    ImproperIterator ii(sys.mol(m).topology());
    for(; ii; ++ii){
      Improper i=ii();
      i[0]+=lastAtom; i[1]+=lastAtom; i[2]+=lastAtom; i[3]+=lastAtom;
      imp.push_back(i);
    }

    lastRes+=sys.mol(m).topology().numRes();
    lastAtom+=sys.mol(m).numAtoms();
  }
}
