#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Hbond.h"
#include "../src/utils/Dssp.h"


#include <iostream>
#include <iomanip>

/**
 *
 * Dssp within gromos++ defines secondary structure for one single (protein) solute  
 * molecule, according to the DSSP rules defined by W. Kabsch and C. Sander 
 * (Biopolymers, 22, pp2577-2637 (1983)). Within these rules it may occur that one 
 * residue is defined as two different secondary-structure elements. In order to 
 * avoid duplicates in the output, the following (ad hoc) priority rules are applied
 * here: Beta Sheet/Bridge > 5-helix > 4-helix > 3-helix > H-bonded turn > Bend. 
 * As a consequence, there may be, for instance, helices that are shorter than 
 * their minimal length. 
 * 
 */

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;
using utils::AtomSpecifier;

int main(int argc, char **argv){
  
  char *knowns[] = {"topo", "pbc", "protein", "time", "nthframe", "traj"};
  int nknowns = 6;
  
  string usage = argv[0];
  usage += "\n\t@topo     <topology>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t@protein  <atomspecifier for the protein>\n";
  usage += "\t@time     <time and dt>\n";
  usage += "\t@nthframe <write every nth frame> (optional, default is 1)\n";
  usage += "\t@traj     <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    
    int nthFrame = 1;
    try{
      args.check("nthframe", 1);
      nthFrame = atoi(args["nthframe"].c_str());
    }
    catch (const gromos::Exception &e){}
    
    // get the protein atoms
    AtomSpecifier prot(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("protein");
      Arguments::const_iterator to=args.upper_bound("protein");
      if (iter == to)
	throw Arguments::Exception("argument @protein is required");
      
      for(;iter!=to;iter++) {
	prot.addSpecifier(iter->second.c_str());
	cout << iter->second.c_str() << " ";
      }
    }
    cout << "\n# In the output the residues are numbered sequentially"
	 << " from 1 to n and not according to @protein!\n";
    
    Dssp SecStr(sys,args);
    
    //get time
    
    double time=0, dt=1; 
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
	time=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
    }
    
    SecStr.settime(time, dt);
    SecStr.calcnumres(prot);
    SecStr.determineAtoms(prot);
    SecStr.calcHintra_init(prot); 
    
    InG96 ic;
    
    int skipFrame = 0;
    
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter) {
      ic.open((iter->second).c_str());
      while(!ic.eof()) {
	ic >> sys;
	if (! skipFrame) {	
	  SecStr.calcHb_Kabsch_Sander();
	  SecStr.calc_Helices();
	  SecStr.calc_Betas();
	  SecStr.calc_Bends();
	  SecStr.filter_SecStruct();
	  SecStr.writeToFiles(nthFrame);
	  SecStr.keepStatistics();
	}
	skipFrame++;
	skipFrame %= nthFrame;
      }
      ic.close();
    }
    SecStr.writeSummary(cout);
    
  }
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




