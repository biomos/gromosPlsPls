//diffus calculates diffusion

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;


double calcD(int fr, int nmol, double s, int nd, double t);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "dim", "atom", "ref", "traj"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atom   <atom to follow>\n";
  usage += "\t@ref    <reference frame (r(0))>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // set atom number
  int at=0;
  {
    Arguments::const_iterator iter=args.lower_bound("atom");
    if(iter!=args.upper_bound("atom"))
      at=atoi(iter->second.c_str())-1;
  }
  
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

  // get the relevant dimensions
  int ndim=3;
  int dim[3]={0,1,2};
  {
    Arguments::const_iterator iter=args.lower_bound("dim");
    if(iter!=args.upper_bound("dim")){
      ndim=0;
      
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
      iter++;
    }
    if(iter!=args.upper_bound("dim")){
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
      iter++;
    }
    if(iter!=args.upper_bound("dim")){
      string dum=iter->second.c_str();
      if(dum=="x")      { dim[ndim]=0; ndim++;}
      else if(dum=="y") { dim[ndim]=1; ndim++;}
      else if(dum=="z") { dim[ndim]=2; ndim++;}
    }
  }
  //  cout << ndim << endl;
  //  for(int i=0;i<ndim;i++) cout << dim[i] << endl;
  
  //  read topology
  args.check("topo",1);
  InTopology it(args["topo"]);

  //  read reference coordinates
  System refsys(it.system());
  InG96 ic;

  try{
    ic.open(args["ref"]);
  }
  // if it didn't work, take the first frame of the trajectory
  catch(const Arguments::Exception &){
    Arguments::const_iterator iter=args.lower_bound("traj");
    if(iter!=args.upper_bound("traj"))
      ic.open((iter->second).c_str());
  }
  ic >> refsys;
  ic.close();

  // we always need the old coordinates to take the nearest image
  System oldsys(refsys);
  
  // and the current system
  System sys(refsys);
  
  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // values to store results
  double d,sum=0.0, Dts=0;
  int frames=0;
  
   ofstream ts; ts.open("diffusts.out");      
  // loop over all trajectories
  for(Arguments::const_iterator 
      iter=args.lower_bound("traj"), to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
      
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
   
      // loop over the atoms to consider
      for(int i=0; i<sys.numMolecules(); i++){

        // gather the relevant atoms with respect to the old system,
	sys.mol(i).pos(at)=pbc->nearestImage(oldsys.mol(i).pos(at),
                                             sys.mol(i).pos(at),
                                             sys.box());
        // calculate difference to refsys for the relevant dimensions
        for(int k=0;k<ndim;k++){
	  d=sys.mol(i).pos(at)[dim[k]]-refsys.mol(i).pos(at)[dim[k]];
	  sum+=d*d;
	}

        // copy the current system to oldsys
        oldsys.mol(i).pos(at)=sys.mol(i).pos(at);
      }
      Dts = calcD(frames, sys.numMolecules(), sum, ndim, time);
      ts << setw(5)  << time
         << setw(20) << Dts*(0.01) << endl;

      frames++;
      time+=dt;
    }
    
    ic.close();
  }
  // calculate the diffusion
  // average sum over time and number of molecules;
  double diff = calcD(frames, sys.numMolecules(), sum, ndim, time);
  cout << diff*(0.01) << " cm^2/s\n";  

  ts.close();  
}
 
  
catch (const gromos::Exception &e){
  cerr << e.what() << endl;
  exit(1);
}
return 0;
}

double calcD(int fr, int nmol, double s, int nd, double t) {
 double ave=s/(fr*nmol);
 double diff=ave/(2*nd*t);

 return diff;
}
