// comtop.cc

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include <fstream>
#include <strstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "param", "solv"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo <topologies>\n";
  usage += "\t@param <number of topology to take parameters from>\n";
  usage += "\t@solv <number of topology to take solvent from>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    // set some values
    args.check("param",1);
    
    Arguments::const_iterator iter=args.lower_bound("param");
    int parnum=atoi(iter->second.c_str());
    iter=args.lower_bound("solv");
    int solnum=atoi(iter->second.c_str());
    int topnum=0;
   
    System sys;

    //ugly solution to the not yet implemented '=' for force fields
    string paramname;
    ostrstream title;
    title << "COMTOP: Combined topology using:\n";
    
       
    for(Arguments::const_iterator iter=args.lower_bound("topo"),
         to=args.upper_bound("topo"); iter!=to; ++iter){
      topnum++;
      
      // read topology
      InTopology it(iter->second.c_str());
    
      for(int j=0;j<it.system().numMolecules();j++)
        sys.addMolecule(it.system().mol(j));
      if(topnum==solnum)
        sys.addSolvent(it.system().sol(0));
      if(topnum==parnum)
        paramname=iter->second.c_str();
      title << topnum << ". " << iter->second.c_str() << endl;
    }
    InTopology it(paramname);
    title << "Parameters from " << parnum 
          << ", solvent from " << solnum;
    title << '\0';
    
    
    OutTopology ot(cout);
    
    ot.setTitle(title.str());
    ot.write(sys,it.forceField());
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




