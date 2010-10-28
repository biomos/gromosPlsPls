/**
 * @file check_box.cc
 * Check box dimensions over a trajectory
 */

/**
 * @page programs Program Documentation
 *
 * @anchor check_box
 * @section check_box Check box dimensions over a trajectory
 * @author @ref th
 * @date 7-6-07
 *
 * To check for the distances between atoms and periodic copies of the other
 * atoms in the system, program check_box can be used. Check_box calculates and
 * writes out the minimum distance between any atom in the central box of
 * the system and any atom in the periodic copies (rectangular box and
 * truncated octahedron are supported).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms" to include in calculation (default: all solute)&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;input coordinate (trajectory) files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  check_box
    @topo  ex.top
    @pbc   r
    @atoms 1:1-30
    @traj  ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/TrajArray.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace fit;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "traj" << "atoms" << "pbc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pbc   <periodic boundary conditions>\n";
  usage += "\t@atoms <atoms to include in calculation (default: all solute)>\n";
  usage += "\t@traj  <input coordinate (trajectory) files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    
    // Systems for calculation
    System sys(it.system());

    System refSys(it.system());

    // expand the coordinates already for tis system
    for(int m=0; m<sys.numMolecules(); m++){
      sys.mol(m).initPos();
    }
    //get the atoms to be included
    utils::AtomSpecifier atoms(sys);
    if(args.count("atoms")>0)
      for(Arguments::const_iterator iter=args.lower_bound("atoms"), 
	    to=args.upper_bound("atoms"); iter!=to; ++iter){
	atoms.addSpecifier(iter->second);
      }
    else
      for(int m=0; m<sys.numMolecules(); m++)
	for(int a=0; a<sys.mol(0).numAtoms(); a++)
	  atoms.addAtom(m,a);
 
    // Periodic images of system
    const int max_num_of_images=14;
    System *compare_periodic[max_num_of_images];
    utils::AtomSpecifier compare_atoms[max_num_of_images];
    
    for (int ii=0; ii<max_num_of_images; ii++) {
      compare_periodic[ii]= new System(sys);
      compare_atoms[ii].setSystem(*compare_periodic[ii]);
      for(int i=0; i<atoms.size(); i++)
	compare_atoms[ii].addAtom(atoms.mol(i), atoms.atom(i));
      
    }
    gmath::Vec displace[max_num_of_images];

    // Coordinates input file
    InG96 ic;

    // Coordinates temporary storage
    TrajArray ta(sys);

  
    // Parse boundary conditions for sys
    int num_of_images=0;
    {
      Arguments::const_iterator iter=args.lower_bound("pbc");
      char b=iter->second.c_str()[0];
      switch(b){
	case 't': 
          num_of_images = 14;
          break;
        case 'r':
          num_of_images = 6;
        break;
        default:
          stringstream msg;
          msg << "Periodic boundary of type " << iter->second << " not supported.";
          throw gromos::Exception("check_box", msg.str());
          break;
      }
    }
    
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);
    
    
    // Distance calculation variables
    gmath::Vec distvec;
    double overall_min_dist2=9.9e12;
    double mindist2=9.9e12;
    
    int frame=0;
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
	ic >> sys; frame++;
	
	(*pbc.*gathmethod)();
        ta.store(sys,0); 
        for (int ii=0; ii<num_of_images; ii++) {
          ta.extract(*compare_periodic[ii],0);
          displace[ii][0]=sys.box().K().abs();
          displace[ii][1]=sys.box().L().abs();
          displace[ii][2]=sys.box().M().abs();
        }
        if (num_of_images >= 6) {
          displace[ 0][0] *= 1.0; displace[ 0][1] *= 0.0; displace[ 0][2] *= 0.0;
          displace[ 1][0] *=-1.0; displace[ 1][1] *= 0.0; displace[ 1][2] *= 0.0;
          displace[ 2][0] *= 0.0; displace[ 2][1] *= 1.0; displace[ 2][2] *= 0.0;
          displace[ 3][0] *= 0.0; displace[ 3][1] *=-1.0; displace[ 3][2] *= 0.0;
          displace[ 4][0] *= 0.0; displace[ 4][1] *= 0.0; displace[ 4][2] *= 1.0;
          displace[ 5][0] *= 0.0; displace[ 5][1] *= 0.0; displace[ 5][2] *=-1.0;
        }
        if (num_of_images == 14) {
          displace[ 6][0] *= 0.5; displace[ 6][1] *= 0.5; displace[ 6][2] *= 0.5;
          displace[ 7][0] *= 0.5; displace[ 7][1] *= 0.5; displace[ 7][2] *=-0.5;
          displace[ 8][0] *= 0.5; displace[ 8][1] *=-0.5; displace[ 8][2] *= 0.5;
          displace[ 9][0] *= 0.5; displace[ 9][1] *=-0.5; displace[ 9][2] *=-0.5;
          displace[10][0] *=-0.5; displace[10][1] *= 0.5; displace[10][2] *= 0.5;
          displace[11][0] *=-0.5; displace[11][1] *= 0.5; displace[11][2] *=-0.5;
          displace[12][0] *=-0.5; displace[12][1] *=-0.5; displace[12][2] *= 0.5;
          displace[13][0] *=-0.5; displace[13][1] *=-0.5; displace[13][2] *=-0.5;
        }
        for (int ii=0; ii<num_of_images; ii++) {
          PositionUtils::translate(compare_periodic[ii],displace[ii]);
        }
	
        mindist2=9.9e12;
	int minatom1=-1, minatom2=-1;
	

	for(int i1=0; i1<atoms.size(); i1++){
	  for(int i2=0; i2<atoms.size(); i2++){
	    for(int image=0; image<num_of_images; image++){
	      distvec = atoms.pos(i1) - compare_atoms[image].pos(i2);
	      if(distvec.abs2() < mindist2){
		mindist2=distvec.abs2();
		minatom1=i1;
		minatom2=i2;
	      }
	    }
	  }
	}

        cout << setw(11) << frame;
        cout.precision(7);
        cout << setw(15) << sqrt(mindist2) 
	     << " # " << atoms.toString(minatom1) << " - " 
	     << atoms.toString(minatom2) << endl;
	
        if (mindist2 < overall_min_dist2) {
          overall_min_dist2=mindist2;
        }
	
      }
      ic.close();
    }
    cout << "OVERALL MIN";
    cout.precision(7);
    cout << setw(15) << sqrt(overall_min_dist2) << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

