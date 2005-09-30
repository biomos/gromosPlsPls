/**
 * @file fit.cc
 */

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/FastRotationalFit.h"
//added
#include "../src/fit/PositionUtils.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace bound;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

// copied from frameout
std::string fileName(int i, std::string const & ext);

int main(int argc, char **argv){
  char *knowns[] = {"topo1", "topo2", "ref", "traj", "pbc1", "pbc2", "atoms1", "atoms2",
		    "skip", "stride", "outformat", "exclude1", "exclude2", "solvent"}; 

  const int nknowns = 14;

  string usage = argv[0];

  usage += "\n\t@topo1         <topology of system 1>\n";
  usage += "\t@topo2         <topology of system 2>\n";
  usage += "\t@pbc1          <boundary conditions> <gather type> of system 1\n";
  usage += "\t@pbc2          <boundary conditions> <gather type> of system 2\n";
  usage += "\t@atoms1        <atomspecifier: atoms to consider for fit from system 1\n";
  usage += "\t@atoms2        <atomspecifier: atoms to consider for fit from system 2\n";
  usage += "\t@ref           <reference coordinates (correspond to system 1)>\n";
  usage += "\t@traj          <trajectory files (correspond to system 2)>\n";
  usage += "\t[@skip         <skip frames at beginning>]\n";
  usage += "\t[@stride       <use only every step frame>]\n";
  usage += "\t[@outformat    <output format. either pdb, g96 (default) or vmdam>]\n"; 
  usage += "\t[@exclude1     <molecules from system 1 that should not be written to output>]\n"; 
  usage += "\t[@exclude2     <molecules from system 2 that should not be written to output>]\n";
  usage += "\t[@solvent      <1 / 2 (take from system 1 or from system 2)>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    // system 1 (reference coordinates)
    // read topology
    InTopology it(args["topo1"]);
    System sys1(it.system());
    
    // system 2 (trajectory files or coordinates)
    // read topology
    InTopology it2(args["topo2"]);
    System sys2(it2.system());

    // system 1 (reference coordinates)
    // read the fit atoms
    AtomSpecifier atoms1(sys1);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms1"),
	to=args.upper_bound("atoms1");
      for( ; iter!=to; ++iter){
	atoms1.addSpecifier(iter->second);
      }
      if(atoms1.size()==0)
	throw gromos::Exception("fit",
				"No fit atoms for system 1 given\n");
    }
   
    // system 2 (trajectory)
    // read the fit atoms
    AtomSpecifier atoms2(sys2);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms2"),
	to=args.upper_bound("atoms2");
      for( ; iter!=to; ++iter){
	atoms2.addSpecifier(iter->second);
      }
      if(atoms2.size()==0)
	throw gromos::Exception("fit",
				"No fit atoms for system 2 given\n");
    }
    if(atoms1.size()!=atoms2.size())
      throw gromos::Exception("fit",
			      "(Number of fit atoms must be the same for both systems.\n)");
    
    //molecules to be excluded
    set<int> excl1, excl2;
    {
      Arguments::const_iterator iter=args.lower_bound("exclude1"),
	to=args.upper_bound("exclude1");
      for(;iter!=to; ++iter){
	std::istringstream is(iter->second);
	int i;
	if (!(is >> i))
	  throw gromos::Exception("fit", "could not read exclude1");
	
	excl1.insert(i-1);
      }
      iter=args.lower_bound("exclude2");
      to=args.upper_bound("exclude2");

      for(;iter!=to; ++iter){
	std::istringstream is(iter->second);
	int i;
	if (!(is >> i))
	  throw gromos::Exception("fit", "could not read exclude2");
	
	excl2.insert(i-1);
      }
    }

    int solvent = 0;
    if (args.count("solvent") > 0){
      std::istringstream is(args["solvent"]);
      if (!(is >> solvent))
	throw gromos::Exception("fit", "could not read solvent");
    }
    
    // or can we say FastRotationalFit frf();?
    // vector<bool> fit_spec, rmsd_spec;
    // FastRotationalFit frf(fit_spec, rmsd_spec);
    FastRotationalFit frf;
    
    int skip=0, stride = 1;
    {
      if(args.count("skip")>0){
	istringstream is(args["skip"]);
	if (!(is >> skip))
	  throw gromos::Exception("fit", "could not read skip");
      }
      if(args.count("stride")>0){
	istringstream is(args["stride"]);
	if (!(is >> stride))
	  throw gromos::Exception("fit", "could not read stride");
      }
    }
    
    // parse boundary conditions
    Boundary *pbc1 = BoundaryParser::boundary(sys1, args, "pbc1");
    Boundary *pbc2 = BoundaryParser::boundary(sys2, args, "pbc2");
    
    // GatherParser
    Boundary::MemPtr gathmethod1 = args::GatherParser::parse(args, "pbc1");
    Boundary::MemPtr gathmethod2 = args::GatherParser::parse(args, "pbc2");

    // create the vector to store the trajectory

    // vector< vector < Vec > > traj;
    // vector< Vec > frame(atoms1.size());
 
    // read reference coordinates...

    if(args.count("ref") <= 0){
      throw gromos::Exception("fit","no reference coordinates specified");
    }
    
    {
      InG96 ic;
      ic.open(args["ref"]);
      if (solvent) ic.select("ALL");
      ic >> sys1;
      ic.close();
    }
    (*pbc1.*gathmethod1)();
    fit::PositionUtils::shiftToCog(&sys1, atoms1);
    
    InG96 ic(skip, stride);
    int framenum=0;
    
    // parse outformat
    bool single_file = false;
    OutCoordinates *oc;
    string ext = ".g96";
    if(args.count("outformat")>0){
      string format = args["outformat"];
      if(format == "pdb"){
	oc = new OutPdb();
	ext = ".pdb";}
      else if(format == "g96"){
        oc = new OutG96S();
        ext = ".g96";}
      else if (format == "vmdam"){
        oc = new Outvmdam();
        ext = ".vmd";
	single_file = true;
      }
      else{
        throw gromos::Exception
	  ("fit","output format "+format+" unknown.\n");
      }
    }
    else{
      oc = new OutG96S();
    }

    ofstream os;
 
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);
      if (solvent) ic.select("ALL");

      // loop over all frames
      while(!ic.eof()){
	ic >> sys2;
	if (ic.stride_eof()) break;
	
	(*pbc2.*gathmethod2)();

	// translational and rotational fit
	PositionUtils::shiftToCog(&sys2, atoms2);
	int err = frf.fit(atoms1, atoms2, sys2);
	
	if(err) {
	  ostringstream os;
	  os << "Error while fitting to the reference structure\n"
	     << "Error code " << err << " in frame number " << framenum+1;
	  
	  throw gromos::Exception("fit", os.str());
	}

	// we need an output system
	gcore::System sys3;
	for(int i=0; i<sys1.numMolecules(); ++i){
	  if(excl1.count(i) == 0)
	    sys3.addMolecule(sys1.mol(i));
	}
	for(int i=0; i<sys2.numMolecules(); ++i){
	  if (excl2.count(i) == 0)
	    sys3.addMolecule(sys2.mol(i));
	}

	if (solvent == 1){
	  sys3.addSolvent(sys1.sol(0));
	  for(int i=0; i < sys1.sol(0).numPos(); ++i){
	    sys3.sol(0).addPos(sys1.sol(0).pos(i));
	  }
	}
	else if (solvent == 2){
	  assert(sys2.numSolvents() == 1);
	  sys3.addSolvent(sys2.sol(0));
	  for(int i=0; i < sys2.sol(0).numPos(); ++i){
	    sys3.sol(0).addPos(sys2.sol(0).pos(i));
	  }
	}
	
	string file=fileName(framenum, ext);
	os.open(file.c_str());
	oc->open(os);

	if (solvent)
	  oc->select("ALL");
	// do we have to add anything to the header???
	oc->writeTitle(file);
	
	*oc << sys3;
	os.close();
	
	framenum++;
      } // loop over frames in trajectory

      ic.close();
    } // loop over trajectories
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}

std::string fileName(int numFrames, std::string const & ext)
{
  ostringstream out;
  string outFile="FRAME";
  if (numFrames < 10){
    out <<outFile<<"_"<<"0000"<<numFrames<<ext;
  }
  else if (numFrames < 100){
    out <<outFile<<"_"<<"000"<<numFrames<<ext;
  }
  else if (numFrames < 1000){
    out <<outFile<<"_"<<"00"<<numFrames<<ext;
  }
  else {
    out <<outFile<<"_"<<"0"<<numFrames<<ext;
  } 
  
  return out.str();
}
