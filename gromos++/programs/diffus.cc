//diffus calculates diffusion

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;



double calcD(int nmol, double s, int nd, double t);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "dim", "atoms", "ref", "traj"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@dim    <dimensions to consider>\n";
  usage += "\t@atoms  <atoms to follow>\n";
  usage += "\t@ref    <reference frame (r(0))>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
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
  ic.select("ALL");
  ic >> refsys;
  ic.close();

  // we always need the old coordinates to take the nearest image
  System oldsys(refsys);
  
  // and the current system
  System sys(refsys);
  
  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // set atom number
  AtomSpecifier at(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++)
      at.addSpecifier(iter->second.c_str());
  }
  // for ease of looping, we make three of these atomspecifiers
  // one for each system
  AtomSpecifier ref_at=at;
  ref_at.setSystem(refsys);
  AtomSpecifier old_at=at;
  old_at.setSystem(oldsys);

  // calculate the com of the reference state
  Vec com0(0.0,0.0,0.0);
  for(int i=0; i<ref_at.size(); i++)
    com0+=*ref_at.coord(i);
  com0/=ref_at.size();
  
  // values to store results
  double d,sum=0.0, disp, Dts=0;
  int frames=1;
  Vec comx;
  
   ofstream ts; ts.open("diffusts.out");      
   ofstream dp; dp.open("diffusdp.out");
   ts << "# Time series of the direct diffusion\n";
   dp << "# Time series of the mean square displacement\n";
   vector<double> tdp;
   vector<double> tt;
   
   
  // loop over all trajectories
  for(Arguments::const_iterator 
      iter=args.lower_bound("traj"), to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      comx=Vec(0.0,0.0,0.0);
      sum=0;
      
      // loop over all atoms to gather with respect to their previous position
      for(int i=0; i<at.size(); i++){
	*at.coord(i) =
	  pbc->nearestImage(*old_at.coord(i),
			    *at.coord(i),
			    sys.box());
	comx+=*at.coord(i);
      }
      comx /= at.size();
      
   
      // loop over the atoms to consider
      for(int i=0; i<at.size(); i++){
        // calculate difference to refsys for the relevant dimensions
	// correct for com 
        for(int k=0;k<ndim;k++){
	  //d=sys.mol(m).pos(a)[dim[k]]-comx[dim[k]]
	  //   -refsys.mol(m).pos(a)[dim[k]]+com0[dim[k]];
          d=(*at.coord(i))[dim[k]]-(*ref_at.coord(i))[dim[k]];
	  
	  sum+=d*d;
	}

        // copy the current system to oldsys
        *old_at.coord(i)=*at.coord(i);
      }
      if(time!=0){
	Dts = calcD(at.size(), sum, ndim, time);
	ts << setw(5)  << time
	   << setw(20) << Dts*(0.01) << endl;
      }

      disp=sum/at.size();
      dp << setw(5)  << time
	 << setw(20) << disp << endl;
      tdp.push_back(disp);
      tt.push_back(time);
      
      frames++;
      time+=dt;
    }
    
    ic.close();
  }
  // calculate the diffusion
  // by doing a least square fit to the average displacement
  double sx=0, sy=0, sxx=0, sxy=0;
  int N=tdp.size();
  
  for(int i=0; i<N ; i++){
    sx+=tt[i];
    sy+=tdp[i];
    sxx+=tt[i]*tt[i];
    sxy+=tt[i]*tdp[i];
  }
  double a = (sxy - sx*sy/N)/(sxx-sx*sx/N);
  double b = -(a*sx - sy)/N;
  double diff = calcD(at.size(), sum, ndim, time);

  cout << "Diffusion is calculated from the mean square displacements:\n";
  cout << "  D = <[r0 - r(t)]^2> / (2*ndim*t)  for t -> inf\n";
  cout << endl;
  cout << "Direct application of this relation:\n";
  cout << "  D = " << diff*(0.01) << " cm^2/s\n"; 
  cout << endl;
  cout << "Least square fit of the mean square displacement:" <<endl;
  cout << "  <[r-r(t)]^2> = " << b << " + " << a << " * t " <<endl;
  cout << "  D = " << a/2/ndim*0.01 << " cm^2/s\n";
  cout << endl;

  dp.close();
  ts.close();  
}
 
  
catch (const gromos::Exception &e){
  cerr << e.what() << endl;
  exit(1);
}
return 0;
}

double calcD(int nmol, double s, int nd, double t) {
 double ave=s/nmol;
 double diff=ave/(2*nd*t);

 return diff;
}
