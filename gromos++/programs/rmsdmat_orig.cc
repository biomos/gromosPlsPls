/**
 * @file rmsdmat_orig.cc
 * create a rmsd matrix over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rmsdmat_orig
 * @section rmsdmat_orig create an rmsd matrix over a trajectory
 * @author @ref vk
 * @date 22. 11. 2004
 *
 * This program has been renamed to @ref rmsdmat and will no longer be 
 * maintained in this form
 *
 * <hr>
 */


// rmsdmat.cc

#include <cassert>
#include <sstream>

#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/TrajArray.h"
#include "../src/bound/Boundary.h"

#include <iomanip>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;
using namespace bound;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "pbc", "type", "prop", "atomsrmsdpos", "atomsfit", "skip", "step"};

  const int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@atomsrmsdpos <atomspecifier: atoms to consider for positional rmsd and/or fit>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t[@atomsfit  <atomspecifier: atoms to consider for fit>]\n";
  usage += "\t[@type <either: posrmsd or property>]\n";
  usage += "\t[@prop <property to calculate RMSD from>]\n";
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

    int type = 0;
    string t;
    try{
      args.check("type", 1); 
      t = args.lower_bound("type")->second.c_str(); 
      if (t == "property") type = 1;
      else if (t == "posrmsd") type = 0;
      else if (t != "property" && t != "posrmsd") type = 2;
    }
    catch (const gromos::Exception &e){
      type = 0;
    }

    if (type == 2) throw gromos::Exception("rmsdmat", "Option '"+t+"' for type not known! Abort!\n");
    
        
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // gather reference system
    (*pbc.*gathmethod)();

    delete pbc;

    Reference reffit(&refSys);
    // System for calculation
    System sys(refSys);
    Reference refrmsd(&refSys);
    AtomSpecifier fitatoms(refSys);
    AtomSpecifier rmsdatoms(sys);
    

    //get rmsd atoms
    {
       Arguments::const_iterator iter = args.lower_bound("atomsrmsdpos");
       Arguments::const_iterator to = args.upper_bound("atomsrmsdpos");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        rmsdatoms.addSpecifier(spec);
       }
    }  

    refrmsd.addAtomSpecifier(rmsdatoms);
    
    //try for fit atoms
    try{
      args.check("fitatoms",1);
     
      {
       Arguments::const_iterator iter = args.lower_bound("atomsfit");
       Arguments::const_iterator to = args.upper_bound("atomsfit");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        fitatoms.addSpecifier(spec);
       }
      }

    }
    catch(const Arguments::Exception &){
      fitatoms = rmsdatoms;
    }
    
    reffit.addAtomSpecifier(fitatoms);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);


    // System for calculation
    PropertyContainer props_sys(sys, pbc);
    PropertyContainer props_ref(refSys, pbc);

    if (type != 0) {
     // what properties. 
    
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
        to=args.upper_bound("prop");
      if(iter==to)
        throw gromos::Exception("rmsdmat", 
                                "no property specified");
      for(; iter!=to; ++iter) {
        props_sys.addSpecifier(iter->second.c_str());
        props_ref.addSpecifier(iter->second.c_str());
      }
    }

    if (props_sys.size() == 0) throw gromos::Exception("rmsdmat", "No property specified! Abort!\n");

    }

    RotationalFit rf(&reffit);
 
    Rmsd rmsd(&refrmsd);
    
    Rmsd::MemPtr rmsdmethod;

    if (type == 0) rmsdmethod = &Rmsd::rmsd;
    else { 
     rmsdmethod = &Rmsd::rmsdproperty;
     rmsd.addproperty(&props_ref, &props_sys);
    }

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

	//pbc call
        (*pbc.*gathmethod)();
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
      ta.extract(refrmsd.sys(), stframe);

      try {
      for (frame=stframe + 1; frame < storedFrameNum; frame++) {
        ta.extract(sys,frame);
  
        rf.fit(&sys);
        try {
          r = (rmsd.*rmsdmethod)(sys);
        } catch (gromos::Exception& e) {
          cerr << e.what() << endl;
          cerr << "Setting rmsd value to 1000000." << endl;
          r = 1000000;  // dirty hack, is there such a thing as MAX_DOUBLE ?
        }

        cout << setw(14) << stframe << " "
             << setw(14) << frame << " "
             << setw(18) << r
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
