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

  Argument_List knowns;
  knowns << "topo" << "pbc" << "coord" << "cutl" << "cutp" << "refpos"
         << "type" << "atominfo";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <topology>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t@coord   <coordinates to base the list on>\n";
  usage += "\t@cutp    <small cutoff>\n";
  usage += "\t@cutl    <large cutoff>\n";
  usage += "\t@refpos  <atomspecifier> or <vector>\n";
  usage += "\t@type    <ATOMIC or CHARGEGROUP>\n";
  usage += "\t@atominfo <write in atominfo style>\n";
    
  try{
    Arguments args(argc, argv, knowns, usage);

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
      sys.mol(numMol).initPos();
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
    cout << "# Making atom lists based on " << args["coord"] << endl << endl;
    cout << "# Reference point is " << (*ref_as.coord(0))[0] << " " 
	 <<  (*ref_as.coord(0))[1] << " " <<(*ref_as.coord(0))[2];
    if(ref_is_atom){
      cout << " (coordinates of atom ";
      if(ref_as.mol(0)<0) cout << "s"; else cout << ref_as.mol(0) + 1;
      cout << ":" << ref_as.atom(0) + 1 << ")";
    }
    cout << endl << endl;
    cout << "# Using a";
    if(t=="ATOMIC") cout << "n atomic";
    else cout << " charge group based";
    cout << " cutoff criterion" << endl << endl;
    
    if(do_short){
      cout << "# Within the short range ( r < " << cutp << " ) there are "
	   << pp.size() << " elements:" << endl << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tprocmt_list: short range\n";
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
	
	for(int i=0; i < pp.size(); ++i){
	  cout << setw(13) << pp.toString(i)
	       << setw(10) << pp.gromosAtom(i)+1
	       << setw(10) << pp.resnum(i)+1
	       << setw(10) << pp.resname(i)
	       << setw(10) << pp.name(i)
	       << setw(12) << pp.iac(i)+1
	       << setw(10) << pp.charge(i)
	       << endl;
	}
	cout << "END\n";
      }
      else{
	vector<string> s = pp.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
      }
      cout << endl << endl;
    }
    if(do_long){
      cout << "# Within the long range ( r < " << cutl << " ) there are "
	   << pl.size() << " elements:" << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tprocmt_list: long (and short) range\n";
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

	for(int i=0; i < pl.size(); ++i){
	  cout << setw(13) << pl.toString(i)
	       << setw(10) << pl.gromosAtom(i)+1
	       << setw(10) << pl.resnum(i)+1
	       << setw(10) << pl.resname(i)
	       << setw(10) << pl.name(i)
	       << setw(12) << pl.iac(i)+1
	       << setw(10) << pl.charge(i)
	       << endl;
	}
	cout << "END\n";
      }
      else{
	vector<string> s = pl.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
      }
      cout << endl << endl;
    }
    if(do_diff){
      cout << "# Within the shell ( " << cutp << " < r < " << cutl 
	   << " ) there are " << diff.size() << " elements:" << endl;

      if (args.count("atominfo") >=0){

	cout << "TITLE\n\tprocmt_list: long (only) range\n";
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


	for(int i=0; i < diff.size(); ++i){
	  cout << setw(13) << diff.toString(i)
	       << setw(10) << diff.gromosAtom(i)+1
	       << setw(10) << diff.resnum(i)+1
	       << setw(10) << diff.resname(i)
	       << setw(10) << diff.name(i)
	       << setw(12) << diff.iac(i)+1
	       << setw(10) << diff.charge(i)
	       << endl;
	}
	cout << "END\n";
	
      }
      else{
	vector<string> s = diff.toString();
	for(unsigned int i=0; i< s.size(); i++){
	  cout << setw(15) << s[i];
	  if(i%5==4) cout << endl;
	}
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

