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

  char *knowns[] = {"topo", "traj", "nthframe"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@nthframe <write every nth frame> (optional, defaults to 1)\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);


    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    InG96 ic;
    OutCoordinates *oc;
    oc = new OutG96();

    int nthFrame = 1;
    try{
      args.check("nthframe", 1);
      nthFrame = atoi(args["nthframe"].c_str());
    }
    catch (const gromos::Exception &e){}


    // loop over all trajectories
    bool isFirstTraj = true;
    int skipFrame = 0;
    oc->open(cout);     
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      if (isFirstTraj){
        oc->writeTitle(ic.title());
        isFirstTraj = false;
      }

      // loop over all frames
      while(!ic.eof()){
        ic >> sys;
        if (! skipFrame){
          *oc << sys;
        }
        skipFrame++;
        skipFrame %= nthFrame;
      }    
    
      ic.close();
    }
    oc->close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
