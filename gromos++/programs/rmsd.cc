// rmsd.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"

#include <vector>
#include <iostream>
#include <iomanip>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", "pbc", "ref", "mol", "time"};
  int nknowns = 8;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@mol <molecules to be considered>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // get simulation time
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

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();


    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // gather reference system
    (*pbc.*gathmethod)();
 
    delete pbc;
    
    Reference ref(&refSys);

    // Adding references
    int added=0;
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<refSys.numMolecules();++i)
	mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
	  it!=args.upper_bound("mol");++it){
	if(atoi(it->second.c_str())>refSys.numMolecules())
	  throw Arguments::Exception(usage);
	mols.push_back(atoi(it->second.c_str())-1);
      }
    // add classes
    for(Arguments::const_iterator it=args.lower_bound("class");
	it != args.upper_bound("class"); ++it){
      for(vector<int>::const_iterator mol=mols.begin();
	  mol!=mols.end();++mol)
	ref.addClass(*mol,it->second);
      added=1;
    }
    // add single atoms
    for(Arguments::const_iterator it=args.lower_bound("atoms");
	it != args.upper_bound("atoms"); ++it){
      int atom=atoi(it->second.c_str())-1, mol=0;
      while(atom >= refSys.mol(mol).numAtoms()){
	atom-=refSys.mol(mol).numAtoms();
	++mol;
	if(mol==refSys.numMolecules())
	  throw Arguments::Exception(usage);
      }
      ref.addAtom(mol,atom);
      added=1;
    }
    // did we add anything at all?
    if(!added)
      throw Arguments::Exception(usage);

    // System for calculation
    System sys(refSys);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&ref);
    Rmsd rmsd(&ref);

    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
	ic >> sys;
	
	(*pbc.*gathmethod)();

	rf.fit(&sys);

	double r = rmsd.rmsd(sys);
	cout.precision(2);
	cout << setw(10) << time;
	cout.precision(5);
	cout << setw(10) << r << endl;
	time+=dt;
      }
      ic.close();
    }

    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

