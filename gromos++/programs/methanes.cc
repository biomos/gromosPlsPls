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
#include <set>
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "atoms", "cutoff", "traj"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@atoms  <atoms to consider>\n";
  usage += "\t@cutoff <distance within which atoms are bound>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // read topology
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

  // read in atoms to look at
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms"), 
      to=args.upper_bound("atoms");
    while(iter!=to){
      atoms.addSpecifier(iter->second);
      iter++;
    }
    
  }
  
  // read in cutoff
  double cutoff=atof(args["cutoff"].c_str());
  double cutoff2=cutoff*cutoff;
  
  
  // Parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // define input coordinate
  InG96 ic;

  // loop over all trajectories
  int count_frame=0;
  std::vector<int> counters(atoms.size());
  std::vector<int> counters_sum(atoms.size());
  for(int i=0; i<atoms.size(); i++){
    counters[i]=0;
    counters_sum[i]=0;
  }
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      
      // first we put all the atoms to look at in a set
      set<int> atom_set;

      // set some values
      for(int i=0; i<atoms.size(); i++) {
	atom_set.insert(i);
	counters[i]=0;
      }
      
      while(atom_set.size()!=0){
	set<int>::iterator it1=atom_set.begin();
	int iatom=*it1;
	set<int> cluster;
	cluster.insert(iatom);
	atom_set.erase(iatom);
	
	for(set<int>::iterator it2=atom_set.begin();
	    it2!=atom_set.end() && atom_set.size()!=0; ++it2){
	  int atom2=*it2;
	  
	  for(set<int>::iterator it3=cluster.begin(), to3=cluster.end();
	      it3!=to3; ++it3){
	    int atom3=*it3;
	    Vec distance=sys.mol(atoms.mol(atom3)).pos(atoms.atom(atom3));
	    distance -= pbc->nearestImage(distance,
					  sys.mol(atoms.mol(atom2)).pos(atoms.atom(atom2)),
					  sys.box());
	    if(distance.abs2()<=cutoff2){
	      //	      cout << "found a matching atom " << *it3 << " " << *it2 << endl;
	      
	      cluster.insert(atom2);
	      atom_set.erase(atom2);
	    }
	  }
	}
	/*if(time==27.5){
	  
	  cout << "cluster ";
	  
	  for(set<int>::iterator it3=cluster.begin(), to3=cluster.end();
	      it3!=to3; ++it3)
	    cout << *it3 << " " ;
	  cout << endl;
	}
	*/
        counters[cluster.size()-1]++;
	
      }
      cout << setw(10) << time;
      
      for(int i=0; i<atoms.size(); i++){
	cout << setw(10) << counters[i] << " ";
	counters_sum[i]+=counters[i];
	counters[i]=0;
      }
      cout << endl;
      
    
      
      count_frame++;
      time+=dt;
    }
    ic.close();
    
  }

  cout << endl;
  cout << "average   ";
  for(int i=0; i<atoms.size(); i++)
    cout << setw(10) << double(counters_sum[i])/count_frame;
  cout << endl;
  
 
}
 
 





 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




























