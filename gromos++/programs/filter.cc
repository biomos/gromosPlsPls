// filter.cc -- remove all solvent molecules that are further than
//              a cutoff distance away from a set of atoms

#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"

using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;
using namespace std;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "traj", "cut", "atoms", "select", "reject", "pairlist", "outformat"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@pbc   <boundary type>\n";
  usage += "\t@cut   <cutoff distance>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@select   <specific atoms to keep>\n";
  usage += "\t@reject   <specific atoms not to keep>\n";
  usage += "\t@pairlist <ATOMIC or CHARGEGROUP>\n";
  usage += "\t@outformat <g96 or pdb>\n";
  usage += "\t@traj  <trajectory files>\n";
   

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    System osys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define in and output coordinates
    InG96 ic;

    // read in the atom list the has to be kept definitely
    utils::AtomSpecifier ls(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("select"),
	to=args.upper_bound("select");
      for(;iter!=to; ++iter)
	ls.addSpecifier(iter->second);
    }    
    // read in the atom list for atoms to throw away for sure
    utils::AtomSpecifier rej(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("reject"),
	to=args.upper_bound("reject");
      for(;iter!=to; ++iter)
	rej.addSpecifier(iter->second);
    }
    // check that there are no doubles in ls and rej
    for(int i=0; i<rej.size(); i++)
      if( ls.findAtom(rej.mol(i), rej.atom(i))!=-1 )
	throw gromos::Exception("filter", "select and reject show overlap");

    // read in the reference atoms for the additional distance criterion
    utils::AtomSpecifier ref(sys);
    double cut=0;
    {
      Arguments::const_iterator iter=args.lower_bound("atoms"),
	to=args.upper_bound("atoms");
      for(;iter!=to; ++iter)
	ref.addSpecifier(iter->second);
    }
    if(ref.size()){
      // then we need a cutoff
      if(args.count("cut")<=0)
	throw gromos::Exception("filter", "If you specify reference atoms, then you need a cutoff");
      cut=atof(args["cut"].c_str());
    }

    // read in the type
    std::string t="ATOMIC";
    if(args.count("pairlist")>0){
      if(args["pairlist"]=="CHARGEGROUP") t=args["pairlist"];
      else if(args["pairlist"]!="ATOMIC") throw gromos::Exception("filter",
                "only ATOMIC and CHARGEGROUP are accepted for pairlist");
    }

    bool pdb_format=false;
    
    if(args.count("outformat")>0){
      if(args["outformat"]=="pdb") pdb_format = true;
    }

    
    // loop over all trajectories
    ostream &os(cout);

    if(!pdb_format){
      os <<"TITLE" << endl;
      
      os << "Filtered trajectory. Keeping" << endl;
      if(ls.size()){
	vector<string> s=ls.toString();
	os << "* atoms ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
      }
      if(ref.size()){
	vector<string> s=ref.toString();
	os << "* atoms within " << cut << " of atoms ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
	os << "  using a";
	if(t=="ATOMIC") os << "n atomic"; else os << " chargegroup based";
	os << " cutoff" << endl;
      }
      if(rej.size()){
	vector<string> s=rej.toString();
	os << "* as long as the above atoms do not belong to ";
	for(unsigned int i=0; i<s.size(); i++) os << s[i] << " ";
	os << endl;
      }
      os << "END" << endl;
    }
    
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while(!ic.eof()){
        ic >> sys;
	(*pbc.*gathmethod)();

	utils::AtomSpecifier rls=ls;
	
	Vec center(0.0,0.0,0.0);
	// add all atoms that need to be added according to the 
	// distances to reference atoms
	for(int i=0; i<ref.size(); i++){
	  utils::SimplePairlist spl(sys, *pbc, cut);
	  spl.setAtom(ref.mol(i), ref.atom(i));
	  spl.setType(t);
	  spl.calc();
	  spl.addAtom(ref.mol(i), ref.atom(i));
	  
	  rls = rls + spl;
	  center += *ref.coord(i);
	  
	}
	if(ref.size()) center/=ref.size();
	
	// remove atoms that are to be rejected
	for(int i=0; i<rej.size(); i++){
	  rls.removeAtom(rej.mol(i), rej.atom(i));
	}
	
	rls.sort();
	if(pdb_format){
	  os.setf(ios::fixed, ios::floatfield);
	  os.setf(ios::unitbuf);
	  os.precision(3);
	  int res=0;
	  int count=0;
	  int resoff=0;
	  
	  for (int i=0; i<rls.size(); ++i){
	    
	    int maxmol=rls.mol(i);
	    if(maxmol<0) maxmol=sys.numMolecules();
	    count=rls.atom(i);
	    resoff=0;
	    for(int j=0; j< maxmol; j++) {
	      count+=sys.mol(j).numAtoms();
	      resoff+=sys.mol(j).topology().numRes();
	    }
	    
	    if(rls.mol(i)<0) res=rls.atom(i)/sys.sol(0).topology().numAtoms();
	    else res=sys.mol(rls.mol(i)).topology().resNum(rls.atom(i));
	    os << "ATOM";
	    os.setf(ios::right, ios::adjustfield);
	    os << setw(7) << count+1;
	    os.setf(ios::left, ios::adjustfield);
	    os << "  " <<setw(4) << rls.name(i).c_str();
	    if(rls.mol(i)<0) os << setw(4) << "SOLV";
	    else os << setw(4) << sys.mol(rls.mol(i)).topology().resName(res).c_str();
	    os.setf(ios::right, ios::adjustfield);
	    Vec pos=pbc->nearestImage(center, *rls.coord(i), sys.box());
	    
	    os << setw(5) << res+resoff+1 << "    "
	       << setw(8) << pos[0]*10
	       << setw(8) << pos[1]*10
	       << setw(8) << pos[2]*10
	       << "  1.00  0.00" << endl;
	  }
	  os << "TER\n";
	}
	else{
	  
	  os  << "POSITIONRED" << endl;
	  os.setf(ios::fixed, ios::floatfield);
	  os.precision(9);
	  os << "# filter selected " << rls.size() << " atoms" << endl;
	  for (int i=0;i<rls.size();++i){
	    os << setw(15) << (*rls.coord(i))[0]
	       << setw(15) << (*rls.coord(i))[1]
	       << setw(15) << (*rls.coord(i))[2] 
	       << "  # ";
	    if(rls.mol(i)<0) os << "s"; else os << rls.mol(i)+1;
	    os << ":" << rls.atom(i)+1 << endl;
	    if(!((i+1)%10))
	      os << "#" << setw(10) << i+1 << endl;
	  }
	  os << "END" << endl;
	
	
	  // and write the box block
	  os << "BOX" << endl;
	  os << setw(15) << sys.box()[0]
	     << setw(15) << sys.box()[1]
	     << setw(15) << sys.box()[2] << endl;
	  os << "END" << endl;
	}
      }
      
      
      ic.close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
