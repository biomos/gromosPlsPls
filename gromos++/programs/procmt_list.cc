// procmt_list writes out an atom specifier that selects all atoms within a
//             cutoff and a list of atoms that are in a shell within rcutp 
//             and rcutl

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace utils;

using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "coord", "cutl", "cutp", "refpos", "type"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo    <topology>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t@coord   <coordinates to base the list on>\n";
  usage += "\t@cutp    <small cutoff>\n";
  usage += "\t@cutl    <large cutoff>\n";
  usage += "\t@refpos  <atomspecifier> or <vector>\n";
  usage += "\t@type    <ATOMIC or CHARGEGROUP>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // define in and output coordinates
    InG96 ic(args["coord"]);
    ic.select("ALL");
    
    ic >> sys;
    
    // read in the cutoffs
    double cutp=0.0, cutl=0.0;
    bool do_short=false, do_long=false, do_diff=false;
    if(args.count("cutp")>0) cutp=atof(args["cutp"].c_str());
    if(args.count("cutl")>0) cutl=atof(args["cutl"].c_str());
    if(cutp!=0.0) do_short=true;
    if(cutl!=0.0) do_long =true;
    do_diff = do_short && do_long;
    if(do_diff && cutp > cutl)
      throw gromos::Exception("procmt_list", "cutp should be shorter than cutl");
    
    if(!do_short && !do_long)
      throw gromos::Exception("procmt_list", 
			      "Please specify at least one of cutp or cutl");
    
    // read in the refpos
    gmath::Vec ref;
    AtomSpecifier ref_as(sys);
    bool ref_is_atom=false;
    int numMol=sys.numMolecules();

    if(args.count("refpos")==3){
      Arguments::const_iterator iter=args.lower_bound("refpos");
      for(int i=0;iter != args.upper_bound("refpos"); ++iter, ++i)
	ref[i]=atof(iter->second.c_str());
      MoleculeTopology mt;
      AtomTopology at;
      at.setName("ref");
      at.setChargeGroup(1);
      mt.addAtom(at);
      sys.addMolecule(mt);
      sys.mol(numMol).pos(0) = ref;
      ref_as.addAtom(numMol,0);
    }
    if(args.count("refpos")==1){
      ref_as.addSpecifier(args["refpos"]);
      if(ref_as.size()!=1)
	throw gromos::Exception("procmt_list",
		"only one atom should be specified with refpos");
      ref_is_atom=true;
      
    }

    // read in the type
    std::string t="ATOMIC";
    if(args.count("type")>0){
      if(args["type"]=="CHARGEGROUP") t=args["type"];
      else if(args["type"]!="ATOMIC") throw gromos::Exception("procmt_list",
		"only ATOMIC and CHARGEGROUP are accepted as type");
    }
    
    SimplePairlist pp(sys, *pbc, cutp);
    SimplePairlist pl(sys, *pbc, cutl);
    AtomSpecifier  diff(sys);
    
    if(do_short){
      pp.setAtom(ref_as.mol(0), ref_as.atom(0));
      pp.setType(t);
      pp.calc();
      if(ref_is_atom)
	pp.addAtom(ref_as.mol(0), ref_as.atom(0));
    }
    if(do_long){
      pl.setAtom(ref_as.mol(0), ref_as.atom(0));
      pl.setType(t);
      pl.calc();
      if(ref_is_atom)
	pl.addAtom(ref_as.mol(0), ref_as.atom(0));
    }
    if(do_diff){
      for(int i=0; i<pl.size(); i++)
	if(pp.findAtom(pl.mol(i), pl.atom(i))==-1)
	  diff.addAtom(pl.mol(i), pl.atom(i));
    }

    // now produce the output
    cout << "Making atom lists based on " << args["coord"] << endl << endl;
    cout << "Reference point is " << (*ref_as.coord(0))[0] << " " 
	 <<  (*ref_as.coord(0))[1] << " " <<(*ref_as.coord(0))[2];
    if(ref_is_atom){
      cout << " (coordinates of atom ";
      if(ref_as.mol(0)<0) cout << "s"; else cout << ref_as.mol(0) + 1;
      cout << ":" << ref_as.atom(0) + 1 << ")";
    }
    cout << endl << endl;
    cout << "Using a";
    if(t=="ATOMIC") cout << "n atomic";
    else cout << " charge group based";
    cout << " cutoff criterion" << endl << endl;
    
    if(do_short){
      cout << "Within the short range ( r < " << cutp << " ) there are "
	   << pp.size() << " elements:" << endl;
      vector<string> s = pp.toString();
      for(unsigned int i=0; i< s.size(); i++){
	cout << setw(15) << s[i];
	if(i%5==4) cout << endl;
      }
      cout << endl << endl;
    }
    if(do_long){
      cout << "Within the long range ( r < " << cutl << " ) there are "
	   << pl.size() << " elements:" << endl;
      vector<string> s = pl.toString();
      for(unsigned int i=0; i< s.size(); i++){
	cout << setw(15) << s[i];
	if(i%5==4) cout << endl;
      }
      cout << endl << endl;
    }
    if(do_diff){
      cout << "Within the shell ( " << cutp << " < r < " << cutl 
	   << " ) there are " << diff.size() << " elements:" << endl;
      vector<string> s = diff.toString();
      for(unsigned int i=0; i< s.size(); i++){
	cout << setw(15) << s[i];
	if(i%5==4) cout << endl;
      }
      cout << endl << endl;
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

