/**
 * @file gridbox.cc
 * @page programs Program Documentation
 *
 * @anchor gridbox
 * @section gridbox generating box
 * @author @ref co
 * @date 30. 11. 2004
 *
 * puts nsm^3 molecules on a grid in a cubic box at the specified density in kg m^-3
 * 
 * arguments:
 * - topo topologies
 * - pos coordinate files
 * - nsm number of molecules (per dimension)
 * - dens system density in kg m^-3 
 * 
 *
 * Example:
 * @verbatim
 gridbox
 @topo ex.top
 @pos exref.coo
 @nsm 8
 @dens 1000
 
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

  char *knowns[] = {"topo", "pos", "nsm", "dens"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pos <coordinates for a single molecule>\n";
  usage += "\t@nsm <number of molecules per dimension>\n";
  usage += "\t@dens <density of liquid (kg/m^3)>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
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

     // new system
    System sys;
    for(int i=0;i<3;i++){
      double *tmp = (double *) &sys.box()[i];
      *tmp = box;
    }
        
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
    os << "Buildbox: " << nsm << " copies of "<<args["insx"]<<endl;
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




