// frameout.cc

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"


using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

bool writeFrame(int i, std::string const & spec, vector<int> const & fnum);
std::string fileName(int i, std::string const & ext);


int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "pbc", "spec", "frames", "outformat", 
		    "include", "ref", "atomsfit"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo       <topology>\n";
  usage += "\t@pbc        <boundary type> <gather method>\n";
  usage += "\t[@spec      <specification for writing out frames. "
    "either ALL (default), EVERY or SPEC>]\n";
  usage += "\t[@frames    <frames to be written out>]\n";
  usage += "\t[@outformat <output format. either pdb, g96 (default) or vmdam>]\n"; 
  usage += "\t[@include   <either SOLUTE (default), SOLVENT or ALL>]\n";
  usage += "\t[@ref       <reference structure to fit to>]\n";
  usage += "\t[@atomsfit  <atoms to fit to>]\n";
  usage += "\t@traj       <trajectory files>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
   
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // do we want to fit to a reference structure
    bool fit=false;
    System refSys(sys);
    Reference reffit(&refSys);
    Vec cog(0.0,0.0,0.0);
    
    if(args.count("ref")>0){
      fit=true;
      // refSys = sys;

      // Parse boundary conditions
      Boundary *pbc = BoundaryParser::boundary(refSys, args);
      Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
   
      // read reference coordinates...
      InG96 ic(args["ref"]);
      ic >> refSys;
      (*pbc.*gathmethod)();
 
      delete pbc;
    
      AtomSpecifier fitatoms(refSys);
      
      //try for fit atoms
      if(args.count("atomsfit") > 0){
	Arguments::const_iterator iter = args.lower_bound("atomsfit");
	Arguments::const_iterator to = args.upper_bound("atomsfit");
	for(;iter!=to;iter++) fitatoms.addSpecifier(iter->second);
      }
      else{
	throw gromos::Exception("frameout", 
				"If you want to fit (@ref) then give "
				"atoms to fit to (@atomsfit)");
      }
      reffit.addAtomSpecifier(fitatoms);
      cog=PositionUtils::cog(refSys, reffit);
    }
    // does this work if nothing is set?
    RotationalFit rf(&reffit);

    // parse includes
    string inc = "SOLUTE";
    if(args.count("include")>0){
      inc = args["include"];
      if(inc != "SOLUTE" && inc !="ALL" && inc!="SOLVENT")
	throw gromos::Exception("frameout",
				"include format "+inc+" unknown.\n");
    }
    
    // parse spec
    string spec = "ALL";
    vector<int> fnum;
    if(args.count("spec")>0){
      spec = args["spec"];
      if(spec!="ALL" && spec !="EVERY" && spec !="SPEC")
	throw gromos::Exception("frameout",
				"spec format "+spec+" unknown.\n");
      if(spec=="EVERY" || spec=="SPEC"){
	//smack in the framenumbers
	for(Arguments::const_iterator it=args.lower_bound("frames");
	    it != args.upper_bound("frames"); ++it){
	  int bla=atoi(it->second.c_str());
	  fnum.push_back(bla);
	}      
	if(fnum.size()==0){
	  throw gromos::Exception("frameout", 
				  "if you give EVERY or SPEC you have to use "
				  "@frames as well");
	}
	if(fnum.size()!=1 && spec=="EVERY"){
	  throw gromos::Exception("frameout",
				  "if you gice EVERY you have to give exactly"
				  " one number with @frames");
	}
      }
    }
    
    // parse outformat
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
        ext = ".vmd";}
      else 
        throw gromos::Exception("frameout","output format "+format+" unknown.\n");
    }
    else{
      oc = new OutG96S();
    }
    
    // loop over all trajectories
    InG96 ic;
    int numFrames = 0;
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);  
      // loop over all frames
      while(!ic.eof()){
	numFrames++;
	ic.select(inc);
	ic >> sys;

	if(writeFrame(numFrames, spec, fnum)){
	  
	  (*pbc.*gathmethod)();
	 
	  if(fit){
	    rf.fit(&sys);
	    PositionUtils::translate(&sys, cog);	  
	  }
	  
	  string file=fileName(numFrames, ext);
	  
	  ofstream os(file.c_str());
	  
	  oc->open (os);       
	  oc->select(inc);
	  oc->writeTitle(file);
	  
	  *oc << sys;
	  os.close();
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

bool writeFrame(int i, std::string const & spec, vector<int> const & fnum)
{
  else if(spec=="SPEC"){
    for(unsigned int j=0; j< fnum.size(); ++j){
      if(fnum[j]==i) return true;
    }
  }
  return false;
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

