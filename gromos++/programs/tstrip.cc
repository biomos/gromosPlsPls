// tstrip.cc

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj"};
  int nknowns = 2;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@traj <trajectory files>\n";



  try{
    Arguments args(argc, argv, nknowns, knowns, usage);


    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

     InG96 ic;
     
     OutCoordinates *oc;

     oc = new OutG96();


    // loop over all trajectories
  int check = 0;
   for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){
    ic.open(iter->second);

if (check == 0){
     oc->open(cout);     
     oc->writeTitle(ic.title());
     oc->close();
          check =1;}

      // loop over all frames
    while(!ic.eof()){
     oc->open(cout);
      
    ic >> sys;

      *oc << sys;
        oc->close();
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
