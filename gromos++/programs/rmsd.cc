/**
 * @file rmsd.cc
 * calculates atom-positional root-mean-square deviations
 */
/**
 * @page programs Program Documentation
 *
 * @anchor rmsd
 * @section rmsd calculates atom-positional root-mean-square deviations
 * @author @ref rb
 * @date 26.7.06
 *
 * The structural deformation of a molecule with respect to a reference
 * structure can be expressed in terms of a root-mean-square deviation (rmsd)
 * of the position of selected atoms. Program rmsd calculates the rmsd over a 
 * molecular trajectory. If requested a least-square rotational fit is performed before
 * to the rmsd calculation. The fit
 * can be performed using a different set of atoms than the calculation of the 
 * rmsd. If no fit is required "no" should be given.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gathermethod&gt;] </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomsrmsd</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsd&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref Atomspecifier "atoms" to consider for fit&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates (if absent, the first frame of \@traj is reference)&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  ditrans
    @topo       ex.top
    @pbc        r
    @time       0 0.1
    @atomsrmsd  1:CA
    @atomsfit   1:CA,C,N
    @ref        exref.coo
    @traj       ex.tr

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <functional>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/args/ReferenceParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/OutG96S.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/groTime.h"
#include "../src/gmath/Vec.h"
#//include "FROGGY/include/gromos++/Vec.h"
///#include "FROGGY/include/gromos++/Vec.h"


using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;
using namespace gmath;

int main(int argc, char **argv){
  Argument_List knowns; 
  knowns << "topo" << "traj" << "atomsfit" << "atomsrmsd" << "pbc" << "ref" 
         << "time" << "list" << "debug" << "fit";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc        <boundary type> [<gathermethod>]\n";
  usage += "\t[@list      <atom_list for gathering>]\n";
  usage += "\t@time       <time and dt>\n";
  usage += "\t@atomsrmsd  <atoms to consider for rmsd>\n";
  usage += "\t[@atomsfit  <atoms to consider for fit>]\n";
  usage += "\t[@ref        <reference coordinates (if absent, the first frame of @traj is reference)>]\n";
  usage += "\t@traj       <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    //System sys(refSys);
    System sys(it.system());
 
    // Parse boundary conditions
    //Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    //Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
   
    // read reference coordinates...
    InG96 ic;
    if(args.count("ref")>0)
      ic.open(args["ref"]);
    else
      if(args.count("traj")>0)
	ic.open(args.lower_bound("traj")->second);

    ic.select("ALL");
    ic >> refSys;
    ic.close();
    //cout << " ref mol " << refSys.numMolecules() << " solv " << refSys.sol(0).numPos() << endl;

    int debug=0;
    if(args.count("debug")>0)
        debug=1;
    int fit=0;
    if(args.count("fit")>0)
        //fit=1;
        fit=atoi(args.lower_bound("fit")->second.c_str());

    // this always goes wrong. check that there is a box block in the refSys
    //if(refSys.hasBox == false && pbc->type()!='v')
    //  throw gromos::Exception("rmsd",
//			      "If pbc != v you have to give a box block "
//			      "in the reference system as well.");
    // and that we managed to read coordinates from it
    if(!refSys.hasPos)
      throw gromos::Exception("rmsd",
                              "Unable to read POSITION(RED) block from "
			      "reference positions file.");

    // DW : read in the atom list for gathering if requested
    Arguments::const_iterator pbciter = args.lower_bound("pbc");
    ++pbciter;

    string gath = pbciter->second;
    cout << "# gather option : " << gath << endl;

    //if(pbciter->second == "1" || pbciter->second == "4"){
    if(gath=="1" || gath == "4"){
        if(args.count("list") <= 0){
            /*throw gromos::Exception("gathering",
                              "request for gathering based on an atom list: "
			      "give the atom list.");
            */
            cout << " ###############  WARNING  ############### "
                    << "# Gathering : You have requested to gather the system based on " << endl
                    << "# an atom list, while you didn't define such a list, therefore "<< endl
                    << "# the gathering will be done according to the 1st atom of the previous molecule" << endl
                    << "# BUT BE AWARE that ++++++ this could be UN-REASONABLE ++++++" << endl;

        } else {
            AtomSpecifier gathlist(sys);
//            AtomSpecifier refelist(refSys);

            if(args.count("list") > 0){
                Arguments::const_iterator iter = args.lower_bound("list");
                Arguments::const_iterator to = args.upper_bound("list");

                for(;iter!=to;iter++){
                    string spec=iter->second.c_str();
                    gathlist.addSpecifierStrict(spec);
                }
                for(int j=0;j<gathlist.size()/2;++j){
                    int i=2*j;
                    sys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                    sys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                    sys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);

                    refSys.primlist[gathlist.mol(i)][0]=gathlist.atom(i);
                    refSys.primlist[gathlist.mol(i)][1]=gathlist.mol(i+1);
                    refSys.primlist[gathlist.mol(i)][2]=gathlist.atom(i+1);
                    //if(debug)
                        //cout << "# updated prim : mol " << gathlist.mol(i) << " atom " << gathlist.atom(i)
                        // << "# refe : mol " << sys.primlist[gathlist.mol(i)][1] << " atom " << sys.primlist[gathlist.mol(i)][2] << endl;
                }
            }
        }
    }
    // end here

    if(gath=="2" || gath=="4"){
        ifstream refframe("REFERENCE.g96");
        if(!refframe){
              gio::OutCoordinates *oref;
              oref = new gio::OutG96S();
              string reffile="REFERENCE.g96";
              ofstream ofile;
              ofile.open(reffile.c_str());
              oref->open(ofile);
              oref->select("ALL");
              oref->writeTitle(reffile);
              *oref << refSys;
              ofile.close();
        }
    }
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    //Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    //if(gath!="2"){
    // gather reference system
    if(gath=="3" || gath=="5"){
        cout << "# reference won't be gathered." << endl;
    }else
    (*pbc.*gathmethod)();
    //}

    delete pbc;

    //cout << "# gathering method : " << gathmethod << endl;
    
    // System for calculation
    Reference refrmsd(&refSys);
    AtomSpecifier fitatoms(refSys);
    AtomSpecifier rmsdatoms(sys);

    //get rmsd atoms
    {
       Arguments::const_iterator iter = args.lower_bound("atomsrmsd");
       Arguments::const_iterator to = args.upper_bound("atomsrmsd");

       for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        rmsdatoms.addSpecifier(spec);
       }
    }  
    if(rmsdatoms.size()==0)
      throw gromos::Exception("rmsd", "No rmsd-atoms specified!");
    
    refrmsd.addAtomSpecifier(rmsdatoms);
    
    //try for fit atoms
    if(args.count("atomsfit") > 0){
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");

      for(;iter!=to;iter++){
        string spec=iter->second.c_str();
        fitatoms.addSpecifier(spec);
      }
    } else {
      cout << "# @atomsrmsd atoms are taken for fit." << endl;
      const vector<string> & spec = rmsdatoms.toString();

      for(vector<string>::const_iterator it = spec.begin(), to = spec.end();
              it != to; ++it) {
        fitatoms.addSpecifier(*it);
      }
    }

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    //Vec cog=PositionUtils::cog(refSys, reffit);

    RotationalFit * rf = NULL;
    if (fitatoms.size()) {
      Reference * ref = new Reference(&refSys);
      ref->addAtomSpecifier(fitatoms);
      rf = new RotationalFit(ref);
    }
    
    Rmsd rmsd(&refrmsd);
    
    int numFrames = 0;

    //cout << rmsd.rmsd(refSys) << endl;
    
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
          
          numFrames++;
          ic.select("ALL");
	ic >> sys >> time;
        if(!sys.hasPos)
          throw gromos::Exception("rmsd",
                             "Unable to read POSITION(RED) block from "
                              "trajectory file.");
	
	(*pbc.*gathmethod)();

        if (fitatoms.size())
          rf->fit(&sys);
                
        double r = rmsd.rmsd(sys);

        // create the reference frame
        /*if((gath=="2" || gath=="4") && numFrames==1){
              cout << "# this frame defined as reference for next frame if any "<< endl;

              gio::OutCoordinates *oref;
              oref = new gio::OutG96S();
              string reffile="REFERENCE.g96";
              ofstream ofile;
              ofile.open(reffile.c_str());
              oref->open(ofile);
              oref->select("ALL");
              oref->writeTitle(reffile);
              *oref << sys;
              ofile.close();

              // the assignment below only for checking whether
              // turning on gathering refering to previous frame
              sys.primlist[0][0] = 31415926;
          }*/
        if((gath=="4") && numFrames==1){
            sys.primlist[0][0] = 31415926;
        }


	cout.precision(2);
	cout << time;
	cout.precision(5);
	cout << setw(10) << r << endl;
         
      }
      ic.close();
    }

    if (rf != NULL) {
      delete rf->getReference();
      delete rf;
    }
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

