// rmsdmat.cc

#include "../src/args/Arguments.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/TrajArray.h"

#include <iomanip>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", "mol", "skip", "step"};

  const int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t[@mol <molecules to be considered>] (defaults to all)\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t[@skip <skip this many frames at the beginning>] (default 0)\n";
  usage += "\t[@step <use only every step frame>] (default 1)\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;
    args.check("traj",1);
    ic.open(args["traj"]);
    ic >> refSys;
    ic.close();

    // Adding references
    Reference ref(&refSys);
    ReferenceParser refP(refSys, args, ref);
    refP.add_ref();

    // skip and step numbers
    unsigned int nSkip;
    try{
      args.check("skip", 1);
      nSkip = atoi(args.lower_bound("skip")->second.c_str());
    }
    catch (const gromos::Exception &e){
      nSkip = 0;
    }
    unsigned int nStep;
    try{
      args.check("step", 1);
      nStep = atoi(args.lower_bound("step")->second.c_str());
    }
    catch (const gromos::Exception &e){
      nStep = 1;
    }

    // System for calculation
    System sys(refSys);

    RotationalFit rf(&ref);
    Rmsd rmsd(&ref);

    // Store coordinates of first molecule in trajectories in array
    unsigned int inFrameNum = 0;
    unsigned int storedFrameNum = 0;
    TrajArray ta(sys);

    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()){
	ic >> sys;

        // check if we want this frame
        if (!((inFrameNum - nSkip) % nStep)){
          rf.fit(&sys);

          // store coordinates from sys in ta
          ta.store(sys,storedFrameNum);
          storedFrameNum++;
        }
      inFrameNum++;
      }
      ic.close();
    }

    double r = 0;
    unsigned int stframe, frame;
    // Fit and calculate rmsd for all pairs with saved trajectory data 
    for (stframe = 0; stframe < storedFrameNum - 1; stframe++) {

      // extract stored coordinates and copy back into ref
      ta.extract(ref.sys(), stframe);

      try {
      for (frame=stframe + 1; frame < storedFrameNum; frame++) {
        ta.extract(sys,frame);
  
        rf.fit(&sys);
        try {
          r = rmsd.rmsd(sys);
        } catch (gromos::Exception& e) {
          cerr << e.what() << endl;
          cerr << "Setting rmsd value to 1000000." << endl;
          r = 1000000;  // dirty hack, is there such a thing as MAX_DOUBLE ?
        }

        cout << setw(8) << stframe
             << setw(8) << frame
             << setw(10) << r
             << endl;
      }
      } catch(gromos::Exception& e) {
        cerr << e.what() << endl; 
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
