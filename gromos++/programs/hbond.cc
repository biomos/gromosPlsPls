/**
 * @file hbond.cc
 * Monitors the occurrence of hydrogen bonds
 */

/**
 * @page programs Program Documentation
 *
 * @anchor hbond
 * @section hbond monitors the occurrence of hydrogen bonds
 * @author @ref MK
 * @date 9-8-2006
 *
 * Program hbond monitors the occurrence of hydrogen bonds over a molecular 
 * trajectory file. It can monitor conventional hydrogen bonds, as well as 
 * three-centered hydrogen bonds through geometric criteria.
 *
 * A hydrogen bond is considered to be present if the distance between a 
 * hydrogen atom, H, connected to a donor atom D, is within a user specified 
 * distance (typically 0.25 nm) from an acceptor atom A and the D-H-A angle is 
 * larger than another user specified value (typically 135 degree). Occurrences
 * of three centered hydrogen bonds are defined for a donor atom D, hydrogen
 * atom H and two acceptor atoms A1 and A2 if (i) the distances H-A1 and H-A2 
 * are within a user specified value (typically 0.27 nm); (ii) the angles 
 * D-H-A1 and D-H-A2 are larger than a second user specified value (typically 
 * 90 degree); (iii) the sum of the angles D-H-A1, D-H-A2 and A1-H-A2 is larger
 * than a third user specified value (typically 340 degree); and (iv) the 
 * dihedral angle defined by the planes through the atoms D-A1-A2 and H-A1-A2 
 * is smaller than a fourth user specified value (typically 15 degree).
 *
 * The user can specify two groups of atoms (A and B) between which the 
 * hydrogen bonds are to be monitored. If hydrogen bond donors and acceptors 
 * are not explicitly specified, these can be filtered based on their masses, 
 * as can be specified in a so-called massfile. If a reference structure is 
 * given, only hydrogen bonds that are observed in the reference structure will
 * be monitored. 
 * 
 * The program calculates average angles, distances and occurrences for all 
 * observed hydrogen bonds over the trajectories and prints out a time series 
 * of the observed hydrogen bonds.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> [\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;] </td></tr>
 * <tr><td> \@DonorAtomsA</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@AcceptorAtomsA</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@DonorAtomsB</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@AcceptorAtomsB</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@Hbparas</td><td>&lt;distance [nm] and angle; default: 0.25, 135&gt; </td></tr>
 * <tr><td> [\@threecenter</td><td>&lt;distances [nm]&gt; &lt;angles&gt; &lt;sum&gt; &lt;dihedral&gt]; </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates for native H-bonds&gt;] </td></tr>
 * <tr><td> [\@massfile</td><td>&lt;massfile&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  hbond
    @topo             ex.top
    @pbc              r
    @time             0 1
    @DonorAtomsA      1:a
    @AcceptorAtomsA   1:a
    @DonorAtomsB      s:a
    @AcceptorAtomsB   s:a
    @Hbparas          0.25 135
    @threecenter      0.27 90 340 15
    @massfile         ../data/hbond.massfile
    @ref              exref.coo
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iostream>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Hbond.h"
#include "../src/utils/Hbond3c.h"
#include "../src/utils/Hbondcalc.h"
#include "../src/utils/Time.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "ref" << "DonorAtomsA" << "AcceptorAtomsA" 
         << "DonorAtomsB" << "AcceptorAtomsB" << "Hbparas" << "threecenter"
         << "time" << "massfile" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t[@time           <time and dt>]\n";
  usage += "\t@DonorAtomsA    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsA <atomspecifier>\n";
  usage += "\t@DonorAtomsB    <atomspecifier>\n";
  usage += "\t@AcceptorAtomsB <atomspecifier>\n";
  usage += "\t@Hbparas        <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t[@threecenter   <distances [nm]> <angles> <sum> <dihedral>]\n";
  usage += "\t[@massfile      <massfile>]\n";
  usage += "\t[@ref           <reference coordinates for native H-bonds>]\n";
  usage += "\t@traj           <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);
    
    InTopology it(args["topo"]);
    System sys(it.system());
    
    Hbondcalc HB(sys,args);
    
    //get time
    Time time(args);
    
    // get the paras
    double maxdist=0.25, minangle = 135;
    {
      Arguments::const_iterator iter=args.lower_bound("Hbparas");
      if(iter!=args.upper_bound("Hbparas")){
	maxdist=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("Hbparas"))
        minangle=atof(iter->second.c_str());
    }
    // and for three-center hydrogen bonds
    double maxdist3c=0.27, minangle3c = 90;
    double minanglesum3c = 340, maxdihedral3c = 15;
    bool hbond3c=false;
    if(args.count("threecenter")>=0){
      hbond3c=true;
      Arguments::const_iterator iter=args.lower_bound("threecenter");
      if(args.count("threecenter")>0){
	maxdist3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>1){
	minangle3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>2){
	minanglesum3c = atof(iter->second.c_str());
	++iter;
      }
      if(args.count("threecenter")>3){
	maxdihedral3c = atof(iter->second.c_str());
      }
    }
    
    HB.setmaxdist(maxdist);
    HB.setminangle(minangle);
    HB.setmaxdist3c(maxdist3c);
    HB.setminangle3c(minangle3c);
    HB.setminanglesum3c(minanglesum3c);
    HB.setmaxdihedral3c(maxdihedral3c);
    HB.determineAtoms();
    
    // print out the information
    if(hbond3c){
      std::cout << "#\n"
		<< "# 2-Centered hydrogen bond D-H..A counted if:\n"
		<< "#     Distance H..A is at most " << maxdist << "\n"
		<< "#     Angle    D-H..A is at least " << minangle << "\n";
      std::cout << "#\n"
		<< "# 3-Centered hydrogen bond D-H..A1\n"
		<< "#                              \\A2 counted if:\n"
		<< "#     Distance H..A1 is at most " << maxdist3c << "\n"
		<< "#     Distance H..A2 is at most " << maxdist3c << "\n"
		<< "#     Angle    D-H..A1 is at least " << minangle3c << "\n"
		<< "#     Angle    D-H..A2 is at least " << minangle3c << "\n"
		<< "#     Sum of angles D-H..A1, D-H..A2, A1..H..A2 is at least "
		<< minanglesum3c << "\n"
		<< "#     Dihedral angle D..A1..A2..H is at most " << maxdihedral3c << "\n";
      std::cout << "#\n#\n";
    }
      
    //check for massfile
    string mfile;
    if(args.count("massfile")>0){
      Arguments::const_iterator iter=args.lower_bound("massfile");
      if(iter!=args.upper_bound("massfile")){
	mfile=(iter->second.c_str());
      }
      HB.readinmasses(mfile);
      HB.determineAtomsbymass();
    }
    

    // initialize the calculation
    HB.init();
    
    // do native?
    bool do_native=false;
    if(args.count("ref")>0){
      InG96 ic(args["ref"]);
      ic.select("ALL");
      ic >> sys;
      ic.close();
      
      // calculate the hb.
      HB.calc();
      // and clear the statistics and reset the time
      HB.clear();
      do_native=true;
    }
    
    InG96 ic;
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      // loop over single trajectory
      
      while(!ic.eof()){
	
	ic >> sys >> time;
        HB.settime(time.time());
        
	if(do_native)
	  HB.calc_native();
	else if(hbond3c)
	  HB.calc3c();
	else
	  HB.calc();
	
	HB.writets();
      }
      ic.close();
      
    }
    
    std::cout << "# Statistics of the run:" << endl;
    if(hbond3c)
      std::cout << endl << "# Two-centered hydrogen bonds:\n";
    std::cout << "#" 
	      << setw(8) << "HB"  
	      << setw(18)<< "Donor" 
	      << setw(15)<< "Acceptor" 
	      << setw(14)<<   "D -"    
	      << setw(14)<< "H ..."
	      << setw(10) << "A" 
	      << setw(8) << "DIST"    
	      << setw(8) << "ANGLE"   
	      << setw(8) << "OCCUR"
	      << setw(8) << "%"       << endl;
    
    
    HB.printstatistics();

    if(hbond3c){
      std::cout << endl << endl
		<< "# Three-centered hydrogen bonds:" << endl
		<< "#" 
		<< setw(82) << " "
		<< setw(12) << "2-CENTER"
		<< setw(30) << "3-CENTER" << endl
		<< "#"
		<< setw(8) << "HB"  
		<< setw(18)<< "Donor" 
		<< setw(15)<< "Acceptor" 
		<< setw(14)<<   "D -"    
		<< setw(14)<< "H ..."
		<< setw(10) << "A" 
		<< setw(8) << "DIST"    
		<< setw(8) << "ANGLE"   
		<< setw(8) << "OCCUR"
		<< setw(8) << "%"       
		<< setw(10) << "SUM"
		<< setw(8) << "DIHED."
		<< setw(8) << "OCCUR"
		<< setw(8) << "%" << endl;
      
      HB.printstatistics3c();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
