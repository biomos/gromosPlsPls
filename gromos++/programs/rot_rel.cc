//Rot_rel calculates the rotational relaxation time

#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "axes", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@axes <two pairs of atoms>\n";
  usage += "\t@traj <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    // get simulation time
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
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // read in the axes to follow
    // In this implementation, we define two axes by giving begin and
    // end point as atom numbers in the molecule. The third axis is defined
    // by the cross-product of these two.
    int ax1[2], ax2[2];
    {
      Arguments::const_iterator iter=args.lower_bound("axes");
      if(iter!=args.upper_bound("axes")){
	ax1[0]=atoi(iter->second.c_str())-1;
	++iter;
      }
      if(iter!=args.upper_bound("axes")){
	ax1[1]=atoi(iter->second.c_str())-1;
	++iter;
      }
      if(iter!=args.upper_bound("axes")){
	ax2[0]=atoi(iter->second.c_str())-1;
	++iter;
      }
      if(iter!=args.upper_bound("axes")){
	ax2[1]=atoi(iter->second.c_str())-1;
      }
    }
    // print out a header
    cout << "# Calculating the autocorrelation function of the legendre"
	 << " polynomials" << endl;
    cout << "# of the dot-product of three molecular axes" << endl;
    cout << "# ax1 points from atom " << ax1[0]+1 << " to atom " << ax1[1]+1
	 << endl;
    cout << "# ax2 points from atom " << ax2[0]+1 << " to atom " << ax2[1]+1
	 << endl;
    cout << "# ax3 is the cross product of these two" << endl;
    cout << "#" << endl;
    
    cout << "#" << setw(9) << "time"
	 << setw(14) << "<p_1(ax1)>"
	 << setw(14) << "<p_2(ax1)>"
	 << setw(14) << "<p_1(ax2)>"
	 << setw(14) << "<p_2(ax2)>"
	 << setw(14) << "<p_1(ax3)>"
	 << setw(14) << "<p_2(ax3)>"
	 << endl;
    
    // prepare vector to store all data.
    vector<vector<Vec> > data(sys.numMolecules()*3);
    Vec v1, v2, v3;
    
    // define input coordinate
    InG96 ic;

    int numFrames=0;
    
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

	//loop over the molecules
	for(int m=0; m<sys.numMolecules(); m++){
	  //calculate the vectors for this molecule
          v1 = sys.mol(m).pos(ax1[1])-sys.mol(m).pos(ax1[0]);
	  v2 = sys.mol(m).pos(ax2[1])-sys.mol(m).pos(ax2[0]);
	  //normalize these vectos
	  v1 = v1.normalize();
	  v2 = v2.normalize();
	  v3 = v1.cross(v2);
	  data[3*m  ].push_back(v1);
	  data[3*m+1].push_back(v2);
	  data[3*m+2].push_back(v3);
	}
	numFrames++;
      }
      ic.close();
    }
    double inp[3], ave1[3], ave2[3];
    
    // now calculate all the autocorrelation functions.
    for(int it=1; it< numFrames; it++){
      double sum1[3]={0.0,0.0,0.0};
      double sum2[3]={0.0,0.0,0.0};
      for(int j=0; j+it < numFrames; j++){
	for(int m=0; m<sys.numMolecules(); m++){
	  for(int k=0; k<3; k++){
	    inp[k]=data[3*m+k][j].dot(data[3*m+k][j+it]);
	    sum1[k]+=inp[k];
	    sum2[k]+=0.5*(3*inp[k]*inp[k]-1);
	  }
	}
      }
      // now print out the information
      cout << setw(10) << time+it*dt;
      for(int k=0; k<3; k++){
        ave1[k]=sum1[k]/numFrames/sys.numMolecules();
        ave2[k]=sum2[k]/numFrames/sys.numMolecules();
	cout << setw(14) << ave1[k]
	     << setw(14) << ave2[k];
      }
      cout << endl;
    }
    
       
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

