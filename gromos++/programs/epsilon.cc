//Dipole - calculates dipole moment and norm with respect to center 

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
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
#include "../src/gmath/physics.h"
#include <vector>
#include <iomanip>
#include <iostream>

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "e_rf", "temp", "traj"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@e_rf <reaction field epsilon>\n";
  usage += "\t@temp <temperature>\n";
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
  
    // read the temperature
    double temp=atof(args["temp"].c_str());
    
    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    //determine net charge
    double ncharge=0;
    double nchargepa=0;
    double totA = 0;
    for (int j=0;j<sys.numMolecules();++j){
      totA += sys.mol(j).topology().numAtoms();
      for (int i = 0; i < sys.mol(j).topology().numAtoms(); i++){
	ncharge+= sys.mol(j).topology().atom(i).charge();
	nchargepa = -(ncharge/totA);
      }
    }
    
    // read e_rf
    double e_rf=1;
    {
      Arguments::const_iterator iter=args.lower_bound("e_rf");
      if(iter!=args.upper_bound("e_rf"))
	e_rf=atof(iter->second.c_str());
    }
    	 
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);

    // prepare the calculation of the average volume
    double vol=0,vsum=0, vave=0, vcorr=1.0;
    Arguments::const_iterator iter=args.lower_bound("pbc");
    if(iter!=args.upper_bound("pbc"))
      if(iter->second[0]=='t') vcorr=0.5;
    
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
    InG96 ic;

    // define a stat class for the size of the dipole
    //    gmath::stat dip_x, dip_y, dip_z;
    // we can of course do it with the statistics class, but in this case
    // there is a) no need for either the average or the error estimate
    // b) no need to store all numbers
    // c) we calculate the msd (and thus loop over the data) every step.
    // so it is probably better to keep the sums ourselves.
 
    int count=0;
    double sum2=0, ave2, fluc;
    Vec sum(0.0,0.0,0.0), ave;
    
    double fac, a, b, eps;
    double f=3.0*EPSILON0*BOLTZ*temp*(2*e_rf + 1.0);

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
        vol=vcorr*sys.box()[0]*sys.box()[1]*sys.box()[2];
	vsum+=vol;
	count++;
        vave=vsum/count;
	Vec dipole(0.0,0.0,0.0);
	for (int j=0; j<sys.numMolecules();++j){
	  for (int i=0; i < sys.mol(j).numAtoms(); i++){
            dipole += (nchargepa+sys.mol(j).topology().atom(i).charge())
	      *sys.mol(j).pos(i);
	  }
	}
        sum2+=dipole.abs2();
	sum+=dipole;
	ave2=sum2/count;
	ave=sum/count;

        fluc=ave2 - ave.abs2();

	fac = vave*f;
        a=2*e_rf*fluc+fac;
	b= -fluc+fac;
	eps = a/b;
	
	cout << setw(10) << time
	  // << setw(14) << dipole.abs()
	  // << setw(14) << ave2
	  // << setw(14) << dipole.abs2()
	  // << setw(14) << ave.abs2()
	  // << setw(14) << fluc
	     << setw(14) << eps
	     << endl;
	
	
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

