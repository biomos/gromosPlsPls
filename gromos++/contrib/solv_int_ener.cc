/**
 * @file solv_int_ener.cc
 * Recalculates interaction energies of specified atoms with the solvent
 */

/**
 * @page programs Program Documentation
 *
 * @anchor solv_int_ener
 * @section solv_int_ener Recalculates interaction energies of specified atoms with the solvent
 * @author @ref ja
 * @date 22-11-2004
 *
 * This is a copy, paste and hack of ener with a few modifications. It can 
 * recalculate interaction energies of the selected atoms with the solvent
 * over molecular trajectory files using the interaction parameters specified
 * in the molecular topology file.
 *
 * The programme prints out the total bonded and nonbonded energies separately,
 * as well as the overall total energy. It is easily modified to print out more
 * detailed energy contributions as well.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; [&lt;gather method&gt;] </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier : atoms for nonbonded interaction&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;time&gt; &lt;dt&gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance&gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field contribution&gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa for reaction field contribution&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  solv_int_ener
    @topo    ex.top
    @pbc     r
    @atoms   1:3-13
    @time    0 0.2
    @cut     1.4
    @eps     61
    @kap     0.0
    @traj    ex.tr
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pbc" << "atoms"  << "time" << "cut"
         << "eps" << "kap" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type> [<gather method>]\n";
  usage += "\t@atoms   <atomspecifier>\n";
  usage += "\t@time    <time> <dt>\n";
  usage += "\t@cut     <cut-off distance>\n";
  usage += "\t@eps     <epsilon for reaction field correction>\n";
  usage += "\t@kap     <kappa for reaction field correction>\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

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

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // declare the energy class
  Energy en(sys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  en.setAtoms(atoms);
  
  // set non-bonded parameters
  //   get cut-off distance
  {
    Arguments::const_iterator iter=args.lower_bound("cut");
    if(iter!=args.upper_bound("cut"))
      en.setCutOff(atof(iter->second.c_str()));
  }
  //  get epsilon and kappa
  {
    double eps=0.0, kap=0.0;
    Arguments::const_iterator iter=args.lower_bound("eps");
    if(iter!=args.upper_bound("eps"))
      eps=atof(iter->second.c_str());
    iter=args.lower_bound("kap");
    if(iter!=args.upper_bound("kap"))
      kap=atof(iter->second.c_str());
    en.setRF(eps, kap);
  }
 
  // define input coordinate
  InG96 ic;
  
  
  // print titles
  cout << "# Time"
       << "                   vdw"
       << "                  elec"
       << "                 Total"
       << endl;

  // declare some variables for averaging
  int num_frames=0;
  double ave_vdw=0.0;
  double ave_elec=0.0;
  double ave_tot=0.0;
  
  
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
      // declare variables for this frame
      double vdw=0.0;
      double elec=0.0;
      
      ic >> sys;
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      pbc->gathergr();

      // calculate the energies
      en.calc();
      
      // loop over selected atoms
      for ( int i=0; i < atoms.size(); ++i) {
        vdw+=en.vdw_s(i);
        elec+=en.el_s(i);
        
      }
      
      // total interaction of selected atoms with solvent
      double tot = vdw + elec;

      // print any ouput you like
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << time
	   << setw(22) << vdw  
           << setw(22) << elec
           << setw(22) << tot
	   << endl;

      //store some averages
      ave_vdw+=vdw;
      ave_elec+=elec;
      ave_tot+=tot;
      
      time+=dt;
      num_frames++;
    }
  }
  // print out averages
  if(num_frames>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << endl << "# ave."
         << setw(22) << ave_vdw/num_frames 
         << setw(22) << ave_elec/num_frames
         << setw(22) << ave_tot/num_frames
         << endl;
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

