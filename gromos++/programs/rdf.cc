//radial distribution function
//chris
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/Reference.h"
#include "../src/fit/PositionUtils.h"
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
#include "../src/gmath/Distribution.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace fit;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "centre", "with", "centrecog",
		    "cut", "grid", "traj"};
  int nknowns = 8;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@centre <atomspecifier>\n";
  usage += "\t[@centrecog] take cog for centre atoms\n";
  usage += "\t@with   <atomspecifier> or\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  // set centre atoms
  AtomSpecifier centre(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    for(;iter!=to;iter++)
      centre.addSpecifier(iter->second.c_str());
  }

  bool cog=false;
  if(args.count("centrecog")>=0) cog=true;
  
  // set atom to consider
  AtomSpecifier with(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    for(;iter!=to;iter++)
      with.addSpecifier(iter->second.c_str());
  }

  // read in cut-off distance
  double cut=1.0;
  if(args.count("cut")>0) cut=atof(args["cut"].c_str());
  
  // read in grid number
  int grid=100;
  if(args.count("grid")>0) grid=atoi(args["grid"].c_str());
  
  // parse boundary conditions
  double vol_corr=1;
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
  if(pbc->type()=='t') vol_corr=0.5;

  // define input coordinate
  InG96 ic;

  // set up distribution arrays
  double rdf[grid];
  double correct=4*acos(-1.0)*cut/double(grid);
  double vol,dens, r;
  
  for(int i=0;i<grid;i++) rdf[i]=0;
  
  // loop over all trajectories
  int count_frame=0;
  
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
      
      // just to make absolutely sure that we treat the centre and with
      // always in the same way, sort them
      centre.sort();
      with.sort();

      // calculate the volume
      vol=sys.box()[0]*sys.box()[1]*sys.box()[2]*vol_corr;
      // loop over the centre atoms
      if(!cog){
	for(int i=0; i<centre.size(); i++){
	  gmath::Distribution dist(0,cut,grid);
	  
	  // to know if this atom is also in the with set.
	  int inwith=0;
	  if(with.findAtom(centre.mol(i),centre.atom(i))>-1) inwith=1;
	  
	  // loop over the atoms to consider
	  for(int j=0; j<with.size();j++){
	    if(!(with.mol(j)==centre.mol(i)&&with.atom(j)==centre.atom(i))){
	      Vec tmp;
	      tmp=pbc->nearestImage(*centre.coord(i),
				    *with.coord(j),
				    sys.box());
	      dist.add((tmp-*centre.coord(i)).abs());
	    }
	  }
	  // now calculate the g(r) for this atom
	  dens=(with.size()-inwith)/vol;
	  for(int k=0; k<grid;k++){
	    r=dist.value(k);
	    rdf[k]+=double(dist[k])/(dens*correct*r*r);
	  }
	}
      }
      else{
	// or if we want to do it for the cog
	(*pbc.*gathmethod)();

	Vec vcog(0.0,0.0,0.0);
	for(int i=0; i<centre.size(); i++)
	  vcog+=*centre.coord(i);
	vcog/=centre.size();

	gmath::Distribution dist(0,cut,grid);

	// loop over the atoms to consider
	for(int j=0; j<with.size();j++){
	  Vec tmp;
	  tmp=pbc->nearestImage(vcog,
				*with.coord(j),
				sys.box());
	  dist.add((tmp-vcog).abs());
	}
	
	// now calculate the g(r) for this atom
	dens=with.size()/vol;
	for(int k=0; k<grid;k++){
	  r=dist.value(k);
	  rdf[k]+=double(dist[k])/(dens*correct*r*r);
	}
      }
      count_frame++;
    }
    ic.close();
  }

  //now correct the distribution for the number of frames and the number 
  //of centre atoms
  cout << "# number of frames considered: " << count_frame << endl;
  int divide=count_frame;
  
  if (!cog) divide*=centre.size();
  
  for(int i=0;i<grid;i++){
    double r=(double(i)+0.5)*cut/grid;
    cout << r << "\t" << rdf[i]/double(divide) << endl;
  }
  
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

