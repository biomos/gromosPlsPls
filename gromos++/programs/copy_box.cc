// buildbox.cc

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
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

  char *knowns[] = {"topo", "insx", "nsm", "dir"};  

  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@insx <coordinates for the molecules>\n";
  usage += "\t@nsm <number of molecules in coordinates>\n";
  usage += "\t@dir <coordinate to duplicate: x/y/z>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    ostringstream title;
    
    // set some values
    args.check("nsm",1);
    Arguments::const_iterator iter=args.lower_bound("nsm");
    int nsm=atoi(iter->second.c_str());
    args.check("dir",1);
    iter=args.lower_bound("dir");
    char dir=iter->second.c_str()[0];
    
    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    int numMol=it.system().numMolecules();
    
    for(int i=1;i<nsm;i++) 
      for(int j=0;j<numMol;j++)
        sys.addMolecule(it.system().mol(j));
    
    
    // read singe atom coordinates...
    InG96 ic;
    args.check("insx",1);
    ic.open(args["insx"]);
    ic.select("ALL");
    ic >> sys;
    title << ic.title();    
    ic.close();
  
    // The following lines gave me a core dump. Probably due to an old
    // boundary parser
    //Boundary *pbc;
    //pbc=new RectBox(&sys);
    //pbc->gather();
    
    
    //calculate shifting vector
    Vec shift;
    if(dir=='x') shift = Vec(sys.box()[0],0,0);
    else if(dir=='y') shift = Vec(0,sys.box()[1],0);
    else if(dir=='z') shift = Vec(0,0,sys.box()[2]);
    else throw gromos::Exception("copy_box", 
         "invalid direction specified, select x,y or z");
    
    //copy and move all molecules
    System sy2(sys);
    PositionUtils::translate(&sy2, shift);
    for(int i=0; i< sy2.numMolecules();i++){
      sys.addMolecule(sy2.mol(i));
    }
    for(int i=0; i< sy2.sol(0).numCoords();i++){
      sys.sol(0).addCoord(sy2.sol(0).pos(i)+shift);
    }
    
    for(int i=0;i<3;i++){
      double *tmp = (double *) &sys.box()[i];
      *tmp = sys.box()[i]+shift[i];
    }

    // Print the new set to cout
    OutG96S oc;
    title << "\nCopy_box: " << args["insx"] << " duplicated in " 
          << dir << "-direction";
    
    oc.open(cout);
    oc.select("ALL");
    
    oc.writeTitle(title.str());
    oc << sys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




