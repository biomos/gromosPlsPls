//Radius of Gyration

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using utils::AtomSpecifier;

int main(int argc, char **argv){
  
  char *knowns[] = {"topo", "pbc", "time", "moln", "traj"};
  int nknowns = 5;
  
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@moln <atom specifier for the atoms to consider>\n";
  usage += "\t@traj <trajectory files>\n";
  
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    //   get simulation time
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
    
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    // get molecules
    cout << "# Radius of Gyration for: \n# ";  
    AtomSpecifier moln(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("moln");
      Arguments::const_iterator to=args.upper_bound("moln");
      for(;iter!=to;iter++) {
        moln.addSpecifier(iter->second.c_str());
	cout << iter->second.c_str() << " ";
      }
    }
    cout << endl;  
    
    // define input coordinate
    InG96 ic;

    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();
	
	//calculate cm, rgyr
	double totalMass=0;
	Vec cm;
	cm[0]=cm[1]=cm[2]=0;
	
	for (int i=0;i < moln.size(); i++) {
	  cm += moln.pos(i) * moln.mass(i);
	  totalMass += moln.mass(i);
	}

	cm = (1.0/totalMass)*cm;
	
	double rg=0;   
	for (int i=0;i < moln.size(); i++) {
	  Vec tmp =  (moln.pos(i)-cm);	
	  rg += tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
	}              
	
	rg = sqrt(rg/(moln.size()));
	
	cout << setw(10) << time ;
	cout << setw(10) << rg << "\n";
	time += dt;
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

