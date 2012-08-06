/**
 * @file cog.cc
 * Calculates centre of geometry or centre of mass for all solute molecules over a 
 * (list of) trajectory file(s)
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cog
 * @section cog Calculates centre of geometry or centre of mass for all solute molecules
 * over a (list of) trajectory file(s)
 * @author @ref dg
 * @date 11-6-07
 *
 * Program cog calculates the centre of geometry (cog) or centre of mass
 * (com) of the solute molecule(s) in the specified frames of the input
 * trajectory file(s), and writes out a single trajectory file in which
 * the position of the cog or com either replaces the atomic coordinates
 * of the solute molecule(s) or are appended directly after the coordinates
 * of the last atom of the solute molecule(s).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary conditions&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;input trajectory file(s)&gt; </td></tr>
 * <tr><td> [\@nthframe</td><td>&lt;write every nth frame (default: 1)&gt;] </td></tr>
 * <tr><td> [\@cog_com</td><td>&lt;calculate centre of geometry (cog) or centre of mass (com) (default: cog)&gt;] </td></tr>
 * <tr><td> [\@add_repl</td><td>&lt;add (add) the cog/com or replace (repl) the solute coordinates with the cog/com (default: repl)&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  cog
    @topo       ex.top
    @pbc        r 
    @traj       ex.tr
    @nthframe   2
    @cog_com    com
    @add_repl   repl
 @endverbatim
 *
 * <hr>
 */


#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/LJException.h"
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

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "traj" << "nthframe" << "cog_com" << "add_repl";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo   <molecular topology file>\n";
  usage += "\t@pbc      <boundary conditions> [<gather method>]\n";
  usage += "\t@traj     <input trajectory files>\n";
  usage += "\t[@nthframe <write every nth frame> (default: 1)]\n";
  usage += "\t[@cog_com  <calculate centre of geometry (cog) or mass (com); default: cog>]\n";
  usage += "\t[@add_repl <add (add) the cog/com or replace (repl) the solutes; default: repl>]\n";
										
										
  try{
    Arguments args(argc, argv, knowns, usage);
    
    
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());
    
    // Boundary and gather method
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // Centre of geometry or mass?
    string inc = "cog";
    if(args.count("cog_com")>0){
      inc = args["cog_com"];
      if(inc!="cog" || inc!="com")
	throw gromos::Exception("cog","include format "+inc+" unknown.\n");
    }
    
    // Add the cog or com, or replace the solute molecules?
    string inc2 = "repl";
    if(args.count("add_repl")>0){
      inc2 = args["add_repl"];
      if(inc2!="repl" || inc2!="add")
	throw gromos::Exception("cog","include format "+inc2+" unknown.\n");
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
    
    int nthFrame = args.getValue<int>("nthframe", false, 1);
    
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
	      if(inc=="cog"){
		cog+=sys.mol(m).pos(a);
	      }
	      else{
		cog+=sys.mol(m).pos(a) * sys.mol(m).topology().atom(a).mass();
		molmass+=sys.mol(m).topology().atom(a).mass();
	      }
	    }
	    
	    if(inc=="cog"){
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
