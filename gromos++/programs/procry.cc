// buildbox.cc

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Matrix.h"
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

void read_spec(std::string name, 
	       vector<Matrix> & rotation, 
	       vector<Vec> &translation, 
	       double factor);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "coord", "spec", "factor"};  

  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@coord  <coordinates for the molecules>\n";
  usage += "\t@spec   <specification file for the translation matrices>\n";
  usage += "\t@factor <conversion factor for distances>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    cout << "read the topology" << endl;
    
    // read the conversion factor
    double factor=1;
    if(args.count("factor")>0) factor=atof(args["factor"].c_str());
 
    cout << "read the factor" << endl;
    
    // read the specification file
    vector<Matrix> rotation;
    vector<Vec> translation;
    read_spec(args["spec"], rotation, translation, factor);

    cout << "read the file" << endl;
    
    // create a final system to work on here.
    // tell it already about the solvent that we have
    System finalsys;
    finalsys.addSolvent(sys.sol(0));
    
    // read single coordinates
    InG96 ic;
    args.check("coord",1);
    ic.open(args["coord"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();    

    // create one more system to keep the coordinates
    System refsys(sys);
    
    int num=rotation.size();
    for(int i=0; i< num; ++i){
      //solute
      for(int m=0; m<sys.numMolecules(); ++m){
	for(int a=0; a<sys.mol(m).numAtoms(); ++a){
	  sys.mol(m).pos(a) = rotation[i]*refsys.mol(m).pos(a) 
	    + translation[i];
	}
	finalsys.addMolecule(sys.mol(m));
      }
      //solvent
      for(int s=0; s<sys.sol(0).numPos(); ++s){
	finalsys.sol(0).addPos(rotation[i]*refsys.sol(0).pos(s)
			       + translation[i]);
      }
    }

    // Print the new set to cout
    OutG96S oc;
    ostringstream title;
    title << "PROCRY did something for you!";

    oc.open(cout);
    oc.select("ALL");
    
    oc.writeTitle(title.str());
    oc << finalsys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




void read_spec(std::string name, 
	       vector<Matrix> & rotation, 
	       vector<Vec> &translation, 
	       double factor)
{
  Ginstream file(name);
  vector<string> buffer;
  file.getblock(buffer);
  file.close();
  
  if(buffer[0]!="TRANSFORM")
    throw gromos::Exception("procry","Could not read TRANSFORM block in specification file");
  if(buffer[buffer.size()-1].find("END")!=0)
    throw gromos::Exception("procry", "Specification file " + file.name() +
			    " is corrupted. No END in "+buffer[0]+
			    " block. Got\n"
			    + buffer[buffer.size()-1]);
  int num=0;
  vector<string>::iterator iter=buffer.begin()+1;

  istringstream is(*iter);
  ++iter;

  if(!(is >> num) || num <=0)
    throw gromos::Exception("procry", "Need some transformations");
  if(buffer.size() - 3 != unsigned (num * 3))
    throw gromos::Exception("procry", "Line count wrong in " +file.name());
  
  Matrix rot(3,3);
  Vec v;

  for(int i=0; i< num; i++){
    for(int j=0; j < 3; j++, ++iter){
      is.clear();
      is.str(*iter);
      for(int k=0; k< 3; ++k){
	if(!(is >> rot(j,k)))
	   throw gromos::Exception("procry", "error reading file");
      }
      if(!(is >> v[j]))
	throw gromos::Exception("procry", "error reading file");
    }
    rotation.push_back(rot);
    translation.push_back(factor*v);
  }
  
}

  
    
  
