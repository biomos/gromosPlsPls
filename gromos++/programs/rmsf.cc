/**
 * @file rmsf.cc
 * calculates atom-positional root-mean-square fluctuations
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rmsf
 * @section rmsf atom-positional root-mean-square fluctuations
 * @author @ref mk 
 * @date 26. 7. 2006
 *
 * Program rmsf calculates atom-positional root-mean-square fluctuations 
 * (rmsf) around average positions for selected atoms over a trajectory. A 
 * rotational fit to a reference structure is performed for every structure in
 * the trajectory. Different sets of atoms can be specified for the fitting 
 * procedure and for the calculation of the rmsf.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@atomsrmsf</td><td>&lt;atomspecifier: atoms to consider for rmsf&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;atomspecifier: atoms to consider for fit&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates(if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  rsmf
    @topo       ex.top
    @pbc        r
    @atomsrmsf  1:CA
    @atomsfit   1:CA,C,N
    @ref        exref.coo
    @traj       ex.tr
@endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "atomsfit", "atomsrmsf", "pbc", "ref"};
  int nknowns = 6;

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t@atomsrmsf   <atomspecifier: atoms to consider for rmsf>\n";
  usage += "\t[@atomsfit   <atomspecifier: atoms to consider for fit>]\n";
  usage += "\t@ref         <reference coordinates(if absent, the first frame of @traj is reference)>\n";
  usage += "\t@traj        <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();


    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    // System for calculation
    System sys(refSys);
    //get the atoms for the rmsf
    AtomSpecifier rmsfatoms(sys); 
    {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsf");
      Arguments::const_iterator to = args.upper_bound("atomsrmsf");
      
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	rmsfatoms.addSpecifier(spec);
      }
    }
    if(rmsfatoms.size()==0)
      throw gromos::Exception("rmsf", 
			      "No atoms specified for RMSF calculation");
    rmsfatoms.sort();

    
    Reference ref(&refSys);

    //get the atoms for the fit
    AtomSpecifier fitatoms(refSys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");
      bool a=(iter==to);
      
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	fitatoms.addSpecifier(spec);
      }      
      if(fitatoms.size()==0){
	if(a==false)
	  cout << "# Warning! The specification for atomsfit resulted in an "
	       << "empty set of atoms.\n"
	       << "#          Taking atomsrmsf instead\n";
	
	for(int i=0; i< rmsfatoms.size(); ++i)
	  fitatoms.addAtom(rmsfatoms.mol(i), rmsfatoms.atom(i));
      }
    }
    ref.addAtomSpecifier(fitatoms);
    

    
    
    
    
    
    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);
    RotationalFit rf(&ref);
    
    
    int numFrames = 0;
    
    //vectors to store the positions, average position and eventually the rmsf value
    vector<Vec> pos;
    vector<Vec> apos;
    vector<double> rmsf;
    
    //init apos for averaging
    Vec zero(0.0,0.0,0.0);
    Vec spos(0.0,0.0,0.0);
    for (int i=0; i < rmsfatoms.size(); ++i) {
      apos.push_back(zero);
      rmsf.push_back(0.0);
    }
    
    //loop over trajectory
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
	numFrames++;
	ic >> sys;
	
	(*pbc.*gathmethod)();
	rf.fit(&sys);
	
	
	//push back the postion    
	for (int i=0; i < rmsfatoms.size(); ++i) {
	  spos = *rmsfatoms.coord(i); 
	  pos.push_back(spos);
	  apos[i] += spos;
	}
	
	
      }
      ic.close();
    } //end loop over trajectory
    
    
    //average positions, write them out
    for (int i=0; i < (int) apos.size(); ++i) apos[i] = apos[i]/numFrames;
    
    //calc rmsf
    for (int i=0, j=0; i < (int) pos.size(); ++i) {
      
      Vec diff = pos[i] - apos[j];
      rmsf[j] += diff.abs2();
      
      ++j;
      
      if (j == (int) rmsf.size())  j = 0;
      
    }
    
    //average
    for (int i=0; i < (int) rmsf.size(); ++i)  rmsf[i] = rmsf[i]/numFrames;
    
    
    //spit out results
    cout << "#\n#  at          rmsf name\n";
     
    for (int i=0; i < rmsfatoms.size(); ++i) {
      cout.precision(8);
      cout << setw(5) << i+1
	   << setw(14) << sqrt(rmsf[i])
	   << setw(5) << rmsfatoms.name(i)
	//sys.mol(rmsfatoms.mol(i)).topology().atom(rmsfatoms.atom(i)).name() 
	   << endl;
    }
    
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


