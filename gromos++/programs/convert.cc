// convert.cc

#include "../src/args/Arguments.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Vacuum.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;
using namespace std;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "atoms", 
		    "pbc", "ref", "mol", "out", "outformat"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to be considered for fit>\n";
  usage += "\t@class <classes of atoms to consider for fit>\n";
  usage += "\t@atoms <atoms to consider to consider for fit>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@outformat <output format>\n";
  usage += "\t@out <output file names>\n";


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

    Reference ref(&refSys);
    try{
      // Adding references
      args.check("ref",1);
      ic.open(args["ref"]);
      ic >> refSys;
      ic.close();
      
      ofstream os("ref.pdb");
      OutPdb oc(os);
      oc << refSys;

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
	for(vector<int>::const_iterator mol=mols.begin();
	    mol!=mols.end();++mol)
	  ref.addClass(*mol,"ALL");

    }
    catch(Arguments::Exception &){}

    // System for calculation
    System sys(refSys);

    // Parse boundary conditions
    Boundary *pbc;
    try{
      char b=args["pbc"].c_str()[0];
      switch(b){
      case 'r':
        pbc=new RectBox(&sys);
        break;
      case 't':
	pbc=new TruncOct(&sys);
	break;
      case 'v':
	pbc=new Vacuum(&sys);
      default:
	throw gromos::Exception("Boundary", args["pbc"] + 
				" unknown. Known boundaries are r, t and v");
      }
    }
    catch(Arguments::Exception &e){
      pbc = new Vacuum(&sys);
    }

    OutCoordinates *oc;
    
    try{
      string format = args["outformat"];
      if(format == "pdb")
	oc = new OutPdb();
      else if(format == "g96")
	oc = new OutG96S();
      else 
	throw gromos::Exception("Convert","output format "+format+" unknown.\n");
    }
    catch(Arguments::Exception &){
      oc = new OutG96S();
    }

    RotationalFit rf(&ref);
    
    
    // loop over all trajectories
    Arguments::const_iterator oit=args.lower_bound("out");

    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
	if(oit==args.upper_bound("out"))
	  throw gromos::Exception("convert","not all frames considered;\ntoo few output file names.");
	cerr << "writing to " << oit->second.c_str() << endl;
	
	ofstream os(oit->second.c_str());
	oc->open(os);
	oc->writeTitle("Converted from " + iter->second);

	ic >> sys;
	pbc->gather();
	
	try{
	  args.check("ref");
	  //	  tf.fitToCog(&sys);
       	  rf.fit(&sys);
	}
	catch (Arguments::Exception &){}
	
	*oc << sys;
	oc->close();
	os.close();
	++oit;
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

