/**
 * @file build_box.cc
 * Generate a condensed phase system on a grid
 */

/**
 * @page programs Program Documentation
 *
 * @anchor build_box
 * @section build_box Generate a condensed phase system on a grid
 * @author @ref co
 * @date 7-6-07
 *
 * When simulating a molecular liquid, a starting configuration for the solvent
 * molecules has to be generated. Program build box generates a cubic box
 * filled with identical solvent molecules which are put on an evenly spaced
 * grid such that the density of the box matches the specified value. Note that
 * to generate a starting configuration for the simulation of a binary mixture,
 * the program @ref bin_box can be used (see section V-2.11). Alternatively, 
 * program @ref ran_box (see section V-2.10) generates a starting configuration for the
 * simulation of mixtures consisting of an unlimited number of components, in
 * which the molecules are oriented randomly.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file for a single molecule&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file for a single molecule&gt; </td></tr>
 * <tr><td> \@nsm</td><td>&lt;number of molecules per dimension&gt; </td></tr>
 * <tr><td> \@dens</td><td>&lt;density of liquid (kg/m^3)&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  build_box
    @topo  ex.top
    @pos   exref.coo
    @nsm   5
    @dens  1000
 @endverbatim
 *
 * <hr>
 */

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "pos" << "nsm" << "dens";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file for a single molecule>\n";
  usage += "\t@pos <input coordinate file for a single molecule>\n";
  usage += "\t@nsm <number of molecules per dimension>\n";
  usage += "\t@dens <density of liquid (kg/m^3)>\n";

  try{
    Arguments args(argc, argv, knowns, usage);
    // set some values
    args.check("nsm",1);
    Arguments::const_iterator iter=args.lower_bound("nsm");
    int nsm3=atoi(iter->second.c_str());
    int nsm=nsm3*nsm3*nsm3;
    args.check("dens",1);
    iter=args.lower_bound("dens");
    double densit=atof(iter->second.c_str());

    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System smol(it.system());

    // and calculate some more values
    double weight=0;
    for(int i=0; i<smol.numMolecules();i++)
      for(int j=0; j< smol.mol(i).numAtoms();j++)
        weight+=smol.mol(i).topology().atom(j).mass();
    
    double vtot=nsm*(weight*1.66056)/densit;
    double box=pow(vtot,1.0/3.0);
    double box3=box/nsm3;
    Vec box32(box3/2.0, box3/2.0, box3/2.0);
     
    // read singe atom coordinates...
    InG96 ic;
    args.check("pos",1);
    ic.open(args["pos"]);
    ic >> smol;
    ic.close();
    
    Vec rc=PositionUtils::com(smol)-box32;
    
    PositionUtils::translate(&smol, -rc);

     // new system, cubic box only!
    System sys;
    sys.box().K()[0] = box;
    sys.box().L()[1] = box;
    sys.box().M()[2] = box;

    // box is by default vacuum
    sys.box().setNtb(gcore::Box::rectangular);
        
    for(int i=0;i<nsm3;i++){
      for(int j=0;j<nsm3;j++){
	for(int k=0;k<nsm3;k++){
	  Vec shift(i*box3,j*box3, k*box3);
	  PositionUtils::translate(&smol,shift);
          for(int q=0;q<smol.numMolecules();q++)
	    sys.addMolecule(smol.mol(q));
	  PositionUtils::translate(&smol, -shift);
	}
      }
    }
    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    os << "build_box: " << nsm << " copies of "<<args["pos"]<<endl;
    os << "Density : " << densit << " kg/m^3\t";
    os << "Molecular weight : " << weight << " u";
    
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << sys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




