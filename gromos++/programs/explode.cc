// buildbox.cc
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/bound/RectBox.h"
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
using namespace bound;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "insx", "nsm", "dist"};  

  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@insx <coordinates for the molecules>\n";
  usage += "\t@nsm <number of molecules per dimension>\n";
  usage += "\t@dist <distance to put between molecules>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    // set some values
    args.check("nsm",1);
    Arguments::const_iterator iter=args.lower_bound("nsm");
    int nsm3=atoi(iter->second.c_str());
    args.check("dist",1);
    iter=args.lower_bound("dist");
    double box3=atof(iter->second.c_str());
    
    int nsm=nsm3*nsm3*nsm3;
    double box=nsm3*box3;
    Vec box32(box3/2.0, box3/2.0, box3/2.0);

    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System smol(it.system());
    for(int i=1;i<nsm;i++) smol.addMolecule(it.system().mol(0));
    
    
    // read singe atom coordinates...
    InG96 ic;
    args.check("insx",1);
    ic.open(args["insx"]);
    ic >> smol;
    ic.close();
    Molecule onemol=smol.mol(1);
    System sonemol;
    sonemol.addMolecule(onemol);
    
    Vec rc=PositionUtils::com(sonemol)-box32;
    
    PositionUtils::translate(&sonemol, -rc);

     // new system
    System sys;
    Boundary *pbc;
    pbc=new RectBox(&smol);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    //gather
    (*pbc.*gathmethod)();
    
    for(int i=0;i<3;i++){
      double *tmp = (double *) &sys.box()[i];
      *tmp = box;
    }
    int count=0;
    
    for(int i=0;i<nsm3;i++){
      for(int j=0;j<nsm3;j++){
	for(int k=0;k<nsm3;k++){
	  Vec shift(i*box3,j*box3, k*box3);    
          Molecule onemol=smol.mol(count++);
          System sonemol;
          sonemol.addMolecule(onemol);
          Vec rc=PositionUtils::com(sonemol)-box32;
          PositionUtils::translate(&sonemol, -rc+shift);
	  sys.addMolecule(sonemol.mol(0));
	}
      }
    }
    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    os << "Explode : " << nsm << " molecules put at intermolecular distance ";
    os << box3 << " nm\n";
    os << "Taken from: " <<args["insx"];
    
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




