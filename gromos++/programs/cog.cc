// Program to write out your trajectories
// in terms of the centre of geometry or
// centre of mass of the solute molecules.
// Solvent molecules are unaffected.
// Cog/Com can either be added,
// or replace the solute coordinates.

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InTopology.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/Box.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace bound;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "traj", "nthframe", "geo_mass", "add_repl"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc      <boundary conditions>\n";
  usage += "\t@traj     <trajectory files>\n";
  usage += "\t@nthframe <write every nth frame> (optional, defaults to 1)\n";
  usage += "\t@cog_com  <calculate centre of geometry (cog) or mass (com); default: cog>\n";
  usage += "\t@add_repl <add (add) the cog/com or replace (repl) the solutes; default: repl>\n";
										
										
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // Boundary and gather method
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // Centre of geometry or mass?
    string inc = "cog";
    if(args.count("cog_com")>0){
      inc = args["cog_com"];
      if(format!="cog" || format!="com")
	throw gromos::Exception("cog","include format "+inc+" unknown.\n");
    }
    
    // Add the cog or com, or replace the solute molecules?
    string inc2 = "repl";
    if(args.count("add_repl")>0){
      inc2 = args["add_repl"];
      if(inc2!="repl" || inc2!="add")
	throw gromos::Exception("cog","include format "+format+" unknown.\n");
    }

    // Construct a system to write out the coordinates
    System osys;
    
    // Create the topology for osys
    // (containing only one-atom molecules of course)
    AtomTopology at;
    at.setName("COG");
    at.setIac(19);
    at.setMass(99);
    at.setCharge(0);
    at.setChargeGroup(1);
    
    MoleculeTopology mt;
    mt.addAtom(at);
    mt.setResNum(1,1);
    mt.setResName(1,"AAP");
    
    osys.addSolvent(sys.sol(0)); // Adds solvent topology info to osysr
    
    // Replace or add?
    if(inc2=="repl"){
      for(int m=0; m<sys.numMolecules(); m++){
	osys.addMolecule(mt);
	osys.mol(m).initPos(); // Creates memory to store coordinates
      }
    }
    else{
      for(int m=0; m<sys.numMolecules(); m++){
	osys.addMolecule(sys.mol(m));
	osys.mol(2*m).initPos();
	osys.addMolecule(mt);
	osys.mol(2*m+1).initPos();
      }
    }
    
    InG96 ic;
    OutG96 oc;
    oc.open(cout);
    
    int nthFrame = 1;
    if(args.count("nthframe")>0)
      nthFrame = atoi(args["nthframe"].c_str());
    
    // loop over all trajectories
    bool isFirstTraj = true;
    int skipFrame = 0;
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      
      ic.open(iter->second);
      
      if (isFirstTraj){
	oc.writeTitle(ic.title());
	isFirstTraj = false;
      }
      
      // loop over all frames
      while(!ic.eof()){
	ic.select("ALL");	
	ic >> sys;
	(*pbc.*gathmethod)();
	if (! skipFrame){
	  // Calculate the centre of geometry or mass of the molecules
	  for(int m=0; m<sys.numMolecules(); m++){
	    Vec cog(0.0,0.0,0.0);
	    double molmass=0;
	    
	    for(int a=0; a<sys.mol(m).numAtoms(); a++){
	      if(inc=="geo"){
		cog+=sys.mol(m).pos(a);
	      }
	      else{
		cog+=sys.mol(m).pos(a) * sys.mol(m).topology().atom(a).mass();
		molmass+=sys.mol(m).topology().atom(a).mass();
	      }
	    }
	    
	    if(inc=="geo"){
	      cog/=sys.mol(m).numAtoms();
	    }
	    else{
	      cog/=molmass;
	    }
	    if(inc2=="repl"){
	      osys.mol(m).pos(0)=cog;
	    }
	    else{
	      for(int a=0; a<sys.mol(m).numAtoms(); a++){
		osys.mol(2*m).pos(a)=sys.mol(m).pos(a);
	      }
	      osys.mol(2*m+1).pos(0)=cog;
	    }
	  }
	  
	  // Add solvent molecules
	  // First, set number of solvent atoms to zero to get rid of solvents
	  // from previous frame
	  osys.sol(0).setNumPos(0);
	  // Add solvent coordinates atom by atom (sol(0) corresponds with solvent of type 1)
	  for(int i=0; i< sys.sol(0).numPos(); ++i){
	    osys.sol(0).addPos(sys.sol(0).pos(i));
	  }
	  
	  // Set box dimensions and write out
	  osys.box()=sys.box();
	  oc.select("ALL");
	  oc << osys;
	}
	
	skipFrame++;
	skipFrame %= nthFrame;
      }    
      
      ic.close();
      
    }
    
    oc.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
