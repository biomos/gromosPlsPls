#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/physics.h"
#include <fstream>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <map>
#include <cmath>

using namespace args;
using namespace gio;
using namespace gcore;

#include "ene_ana.h"


int main(int argc, char **argv){

  char *knowns[] = {"topo", "time", "en_files", "fr_files", "prop", "library"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@time    <t and dt>\n";
  usage += "\t@en_files   <energy files>\n";
  usage += "\t@fr_files  <free energy files>\n";
  usage += "\t@prop    <properties to monitor>\n";
  usage += "\t@library <library for property names> [print]\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // calculate total mass of the system
    double mass=0;
    for(int m=0; m<sys.numMolecules(); m++){
      for(int a=0; a< sys.mol(m).topology().numAtoms(); a++){
	mass += sys.mol(m).topology().atom(a).mass();
      }
    }

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

    // read a library file?
    string library;
    int do_library=0;
    int print_library=0;
    
    {
      Arguments::const_iterator iter=args.lower_bound("library"), 
	to=args.upper_bound("library");
      if(iter!=to){
	do_library=1;
	library=iter->second;
	++iter;
      }
      if(iter!=to && iter->second == "print") print_library=1;
    }
    
    // which properties do we follow
    vector<string> prop;
    int num_prop=0;
    
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
	to=args.upper_bound("prop");
      while(iter!=to){
	prop.push_back(iter->second);
	++iter;
      }
      num_prop=prop.size();
    }

    // define an energy trajectory
    energy_trajectory etrj;

    // learn about the variable names how they map to the elements
    set_standards(etrj, mass);
    if(do_library) read_library(library, etrj);
    if(print_library) etrj.write_map();

    // prepare for the statistical information
    gmath::stat s[num_prop];

    // loop over the files
    // as we can either specify an energy or a free energy file or both 
    // (and then they need not be the same length (well?)) this is a bit
    // different

    // define two input streams
    Ginstream gin_en;
    Ginstream gin_fr;
    bool do_energy_files     =(args.count("en_files")>0);
    bool do_free_energy_files=(args.count("fr_files")>0);
    
    Arguments::const_iterator it_en=args.lower_bound("en_files"),
      to_en=args.upper_bound("en_files"),
      it_fr=args.lower_bound("fr_files"),
      to_fr=args.upper_bound("fr_files");
    int cont=0, en_cont=0, fr_cont=0;
    if(do_energy_files) {
      gin_en.open(it_en->second.c_str()); 
      ++it_en; 
      en_cont=1;
    }
    
    if(do_free_energy_files) {
      gin_fr.open(it_fr->second.c_str());
      ++it_fr;
      fr_cont=1;
    }
    
    cont=en_cont+fr_cont;
    
    while(cont){
      
      // read the numbers into the energy trajectory
      if(do_energy_files) etrj.read_energy_frame(gin_en);
      if(do_free_energy_files) etrj.read_free_energy_frame(gin_fr);
      
      // calculate and store the necessary number in the stat-classes
      for(int i=0; i<num_prop; i++)
	s[i].addval(etrj[prop[i]]);
      
      // check if we continue or possibly have to open a new files
      if(do_energy_files && gin_en.eof()){
	if(it_en!=to_en){
	  gin_en.close();
	  gin_en.clear();
	  gin_en.open(it_en->second.c_str());
	  ++it_en;
	}
	else en_cont =0;
      }
      if(do_free_energy_files && gin_fr.eof()){
	if(it_fr!=to_fr){
	  gin_fr.close();
	  gin_fr.clear();
	  gin_fr.open(it_fr->second.c_str());
	  ++it_fr;
	}
	else fr_cont =0;
      }
      cont=en_cont + fr_cont;
      
      if(do_energy_files && do_free_energy_files)
	if(en_cont != fr_cont)
	  cerr << "Energy files and free energy files do not contain the "
	       << "same number of frames!" << endl;
    }
    //print out the statistical information
    cout << setw(10) << "property"
	 << setw(14) << "average"
	 << setw(14) << "rmsd"
	 << setw(14) << "error est."
	 << endl;
    for(int i=0; i<num_prop; i++)
      print(s[i], prop[i], time, dt);
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
