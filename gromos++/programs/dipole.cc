//Dipole - calculates dipole moment and norm with respect to center 

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <iostream>

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "moln", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@moln <number of molecules to be considered>\n";
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
  
  // set molecule number
   double  moln=0;
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     moln=atof(iter->second.c_str());
   }
    }
  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());

  //determine net charge
  double ncharge=0;
  double nchargepa=0;
  double totA = 0;
  for (int j=0;j<moln;++j){
     totA += sys.mol(j).numAtoms();
   for (int i = 0; i < sys.mol(j).numAtoms(); i++){
    ncharge+= sys.mol(j).topology().atom(i).charge();
    nchargepa = -(ncharge/totA);
   }
  }
    
    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

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
   
      Vec dipole;
      for (int j=0; j<moln;++j){
       for (int i=0; i < sys.mol(j).numAtoms(); i++){
            dipole += (nchargepa+sys.mol(j).topology().atom(i).charge())*sys.mol(j).pos(i);
            }
      }
      double mag = sqrt(dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2]);

       cout << setw(10) << time ;
       cout << setw(10) << mag << ' ' << dipole[0] << ' ' << dipole[1] << ' ' << dipole[2] <<  endl;
       time += dt;
    }
    ic.close();
    }
    cout << endl;
    cout << "Netcharge of Molecule " << moln << " is " << ncharge << endl;
    cout << "Netcharge per atom of Molecule " << moln << " is " << -nchargepa << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
