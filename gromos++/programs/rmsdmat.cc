// rmsdmat.cc

#include "../src/args/Arguments.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/TrajArray.h"

#include <vector>
#include <strstream>
#include <iomanip>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", "mol"};

  const int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@mol <molecules to be considered> (defaults to all)\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@traj <trajectory files>\n";

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
    Arguments::const_iterator iter;
    bool added = false;
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol") == args.upper_bound("mol"))
      for(int i = 0; i < refSys.numMolecules(); ++i)
        mols.push_back(i);
    else{
      for(iter = args.lower_bound("mol"); iter != args.upper_bound("mol");
        ++iter){
        int molNum = atoi(iter->second.c_str());
        if(molNum > refSys.numMolecules()){
          string errmsg = usage;
          errmsg += "\n";
          errmsg += "Supplied molecule index is larger than the ";
          errmsg += "number of molecules in the system.";
          throw Arguments::Exception(errmsg);
        }
        mols.push_back(molNum - 1);
      }
    }
    // add classes
    for(iter = args.lower_bound("class"); iter != args.upper_bound("class");
      ++iter){
      vector<int>::const_iterator mol;
      for(mol = mols.begin(); mol != mols.end(); ++mol)
        ref.addClass(*mol,iter->second);
      added = true;
    }
    // add single atoms
    for(iter = args.lower_bound("atoms"); iter != args.upper_bound("atoms");
      ++iter){
      int atom = atoi(iter->second.c_str())-1;
      int mol = 0;
      while(atom >= refSys.mol(mol).numAtoms()){
        atom -= refSys.mol(mol).numAtoms();
        ++mol;
        if(mol == refSys.numMolecules()){
          string errmsg = usage;
          errmsg += "\n";
          errmsg += "Supplied atom index is larger than the ";
          errmsg += "number of atoms in the system.";
          throw Arguments::Exception(errmsg);
        }
      }
      ref.addAtom(mol,atom);
      added = true;
    }
    // did we add anything at all?
    if(!added){
      string errmsg = usage;
      errmsg += "\n";
      errmsg += "Either \"class\" or \"atom\" must be non-empty.\n";
      throw Arguments::Exception(errmsg);
    }

    // System for calculation
    System sys(refSys);

    RotationalFit rf(&ref);
    Rmsd rmsd(&ref);

    // Store coordinates of first molecule in trajectories in array
    unsigned int frame_number=0;
    TrajArray ta(sys);

    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);

      // loop over all frames
      while(!ic.eof()){
	ic >> sys;

	rf.fit(&sys);

        // store coordinates from sys in ta
        ta.store(sys,frame_number);
        frame_number++;
      }
      ic.close();
    }

    double r = 0;
    unsigned int stframe, frame;
    // Fit and calculate rmsd for all pairs with saved trajectory data 
    for (stframe = 0; stframe < frame_number - 1; stframe++) {

      // extract stored coordinates and copy back into ref
      ta.extract(ref.sys(), stframe);

      for (frame=stframe + 1; frame < frame_number; frame++) {
        ta.extract(sys,frame);
  
        rf.fit(&sys);

        r = rmsd.rmsd(sys);

        cout << setw(5) << stframe
             << setw(5) << frame
             << setw(11) << r
             << endl;
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
