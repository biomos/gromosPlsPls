// comtop.cc

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutTopology.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "param", "solv"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo  <topologies>\n";
  usage += "\t@param <number of topology to take parameters from>\n";
  usage += "\t@solv  <number of topology to take solvent from>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    // set some values
    int parnum=1;
    if(args.count("param")>0) parnum=atoi(args["param"].c_str());
    int solnum=1;
    if(args.count("solv")>0)  solnum=atoi(args["solv"].c_str());
    if(args.count("topo")<=0)
      throw gromos::Exception("COMTOP", "needs at least one topology\n"+usage);
    
    System sys;

    //ugly solution to the not yet implemented '=' for force fields
    string paramname, toponame, s;
    std::string::size_type s_it;
    int repeat = 1, topnum=0, oldtopnum=-1;
    
    ostringstream title;
    title << "COMTOP: Combined topology using:\n";
    
       
    for(Arguments::const_iterator iter=args.lower_bound("topo"),
         to=args.upper_bound("topo"); iter!=to; ++iter){
      oldtopnum=topnum+1;
      
      s=iter->second;
      s_it=s.find(':');
      
      if(s_it == string::npos){
	toponame=s;
	repeat=1;
      }
      else{
	toponame=s.substr(s_it+1,s.size());
	repeat=atoi(s.substr(0,s_it).c_str());
      }
      topnum+=repeat;
      
      // read topology
      InTopology it(toponame);
      for(int i=0; i<repeat; i++)
	for(int j=0;j<it.system().numMolecules();j++)
	  sys.addMolecule(it.system().mol(j));

      if(solnum <= topnum && solnum >= oldtopnum)
	sys.addSolvent(it.system().sol(0));
      if(parnum <= topnum && parnum >= oldtopnum)
        paramname=toponame;
      if(topnum!=oldtopnum)
	title << setw(4) << oldtopnum << " .. " << setw(4) << topnum;
      else
	title << setw(12) << topnum;
      
      title << " : " << toponame << endl;
    }
    InTopology it(paramname);
    title << "Parameters from " << parnum 
          << ", solvent from " << solnum;
    
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




