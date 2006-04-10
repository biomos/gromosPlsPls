//minimum distance function
//chris

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace std;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "centre", "with", "nsm", "traj"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@centre <atomspecifier>\n";
  usage += "\t@with   <atomspecifier>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

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

  // set centre atoms
  AtomSpecifier centre(sys);
  
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    for(; iter!=to; ++iter){
      centre.addSpecifier(iter->second);
    }
  }
  // set atom to consider
  AtomSpecifier with(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    for(; iter!=to; ++iter){
      with.addSpecifier(iter->second);
    }
  }
  
  // Parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // define input coordinate
  InG96 ic;

  // loop over all trajectories
  int count_frame=0;
  
  // open centre.size() files
  vector<ofstream *> fout(centre.size());
  
  for(int i=0; i<centre.size(); i++){
    ostringstream os;
    if(centre.mol(i)!=-3)
      os << "MIN_" << centre.mol(i)+1 << ":" << centre.atom(i)+1 
	 << ".dat";
    else
      os << "MIN_va_" << i+1 <<".dat";
    
    fout[i] = new ofstream(os.str().c_str());
    
  }
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
   
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      pbc->gather();
      
      // loop over the centre atoms
      int start=0;
      int minat=0;
      
      // now really loop over the centre atoms
      for(int i=start;i<centre.size();i++){
	double min2=sys.box()[0]*sys.box()[0];
	//loop over the atoms to consider
        for(int j=0;j<with.size();j++){
          //calculate distance only if this atom is not the current centre
          if(!(with.mol(j)==centre.mol(i)&&with.atom(j)==centre.atom(i))){
	    Vec tmp=pbc->nearestImage(centre.pos(i),
				      with.pos(j),
				      sys.box());
	    tmp-=centre.pos(i);
	    
            if(tmp.abs2()<min2) {
              min2=tmp.abs2();
              minat=j;
	    }
	    
	  }
	}
	//write out min dist
        (*fout[i]) << time << "\t" << sqrt(min2) << "\t# " << with.toString(minat) << "\t" << with.atom(minat)+1 << endl;
	
      }
      count_frame++;
      time+=dt;
      
    }
    ic.close();
    
  }
  
  //close the files
  for(int i=0;i<centre.size();i++){
    fout[i]->close();
    delete fout[i];
    fout[i] = NULL;
  }
  
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




























