//distributions dist

#include "../src/args/Arguments.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Distribution.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "moln", "atoms", "dist", "traj"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@dist   <lower and upper boundary and number of steps>\n";
  usage += "\t@moln   <molecule numbers of the corresponding atoms [1...X]>\n";
  usage += "\t@atoms  <atom numbers in molecule to calculate quantity from [1..X]>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // set atom number
  int k=0,l=0,m=0,n=0;
  int num=1;
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    if(iter!=args.upper_bound("atoms")){
      k=atoi(iter->second.c_str())-1;
      ++iter; num =1;
    }
    if(iter!=args.upper_bound("atoms")){
      l=atoi(iter->second.c_str())-1;
      ++iter; num=2;
    }
    if(iter!=args.upper_bound("atoms")){
      m=atoi(iter->second.c_str())-1;
      ++iter;num=3;
    }
    if(iter!=args.upper_bound("atoms")) {
      n=atoi(iter->second.c_str())-1; num=4;
    }  
  }
  //   get distribution parameters
  double begin=0, end=0;
  int nsteps=0;
 
  {
    Arguments::const_iterator iter=args.lower_bound("dist");
    if(iter!=args.upper_bound("dist")){
      begin=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      end=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("dist")){
      nsteps=atoi(iter->second.c_str());
    }     
  }

  // set molecule number
  int  molk=0, moll=0, molm=0, moln=0;

  {
    Arguments::const_iterator iter=args.lower_bound("moln");
    if(iter!=args.upper_bound("moln")){
      molk=atoi(iter->second.c_str())-1;
      ++iter;
    }
    if(iter!=args.upper_bound("moln")){
      moll=atoi(iter->second.c_str())-1;
      ++iter;
    }
    if(iter!=args.upper_bound("moln")){
      molm=atoi(iter->second.c_str())-1;
      ++iter;
    }
    if(iter!=args.upper_bound("moln"))
      moln=atoi(iter->second.c_str())-1;
  }

  //  read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
    
  // Parse boundary conditions
  Boundary *pbc;
  try{
    char b=args["pbc"].c_str()[0];
    switch(b){
      case 't':
        pbc=new TruncOct(&sys);
        break;
      case 'v':
        pbc=new Vacuum(&sys);
        break;
      case 'r':
        pbc=new RectBox(&sys);
        break;
      default:
        throw gromos::Exception("Boundary", args["pbc"] + 
				" unknown. Known boundaries are t, r and v");
	
    }
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  // define input coordinate
  InG96 ic;

  // set up distribution arrays
  gmath::Distribution dist(begin, end, nsteps);

  // set pi
  const double pi = 3.1415926535898;
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
       // ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      pbc->gather();
   
      double val = 0;
      if (num ==1) 
        throw gromos::Exception("dist", " at least two atoms are needed.\n");
      else if (num == 2) {  
        Vec tmp = (sys.mol(molk).pos(k)-sys.mol(moll).pos(l)); 
        val = tmp.abs();
      }
      else if (num == 3) {
        Vec tmpA = (sys.mol(molk).pos(k)-sys.mol(moll).pos(l)); 
        Vec tmpB = (sys.mol(molm).pos(m)-sys.mol(moll).pos(l));
        val = acos((tmpA.dot(tmpB))/(tmpA.abs()*tmpB.abs()))*180/pi;
      }          
      else if (num == 4) {
        Vec tmpA = (sys.mol(molk).pos(k)-sys.mol(moll).pos(l));
        Vec tmpB = (sys.mol(moln).pos(n)-sys.mol(molm).pos(m));
        Vec tmpC = (sys.mol(molm).pos(m)-sys.mol(moll).pos(l));
        Vec p1 = tmpA.cross(tmpC);
        Vec p2 = tmpB.cross(tmpC);

        double cosphi = ((p1.dot(p2))/(p1.abs()*p2.abs()));
 
        val = acos(cosphi)*180/pi;
 
        Vec p3 = p1.cross(p2);
        if (p3.dot(tmpC)<0)
          val = 360 - val;
      }
      
      // count this value in the distribution array
      dist.add(val);
    }

  }
  ic.close();
  // print out the distribution, calculate the average and rmsd
  cout << "\nnumber of values calculated: " << dist.nVal() << endl;
  cout << "average value:               "   << dist.ave() << endl;
  cout << "RMSD (from distribution):    "   << dist.rmsd() << endl;
  dist.write(cout);
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
