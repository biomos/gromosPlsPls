//ener calculates (non-bonded) interaction energies for specific atoms

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Energy.h"
#include "../src/gmath/Vec.h"
#include <iomanip>
#include <iostream>
#include <strstream>
#include <fstream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "atoms", "props", "time", "cut", 
                    "eps", "kap", "soft", "softpar", "traj"};
  int nknowns = 11;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <atomspecifier>\n";
  usage += "\t@props <propertyspecifier>\n";
  usage += "\t@time <time> <dt>\n";
  usage += "\t@cut <cut-off distance>\n";
  usage += "\t@eps <epsilon for reaction field correction>\n";
  usage += "\t@kap <kappa for reaction field correction>\n";
  usage += "\t@soft <atom specifier for soft atoms>\n";
  usage += "\t@softpar <lam> <a_lj> <nkt>\n";
  usage += "\t@traj  <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //   get simulation time
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

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // declare the energy class
  Energy en(sys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  en.setAtoms(atoms);
  
  // set properties
  PropertyContainer props(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("props");
    Arguments::const_iterator to=args.upper_bound("props");
    for(;iter!=to;iter++){
      string p=iter->second.c_str();
      props.addSpecifier(p);
    }
  }
  en.setProperties(props);

  // set non-bonded parameters
  //   get cut-off distance
  {
    Arguments::const_iterator iter=args.lower_bound("cut");
    if(iter!=args.upper_bound("cut"))
      en.setCutOff(atof(iter->second.c_str()));
  }
  //  get epsilon and kappa
  {
    double eps=0.0, kap=0.0;
    Arguments::const_iterator iter=args.lower_bound("eps");
    if(iter!=args.upper_bound("eps"))
      eps=atof(iter->second.c_str());
    iter=args.lower_bound("kap");
    if(iter!=args.upper_bound("kap"))
      kap=atof(iter->second.c_str());
    en.setRF(eps, kap);
  }
  // get soft atom list
  AtomSpecifier soft(sys);
  {
    int lsoft=0;
    Arguments::const_iterator iter=args.lower_bound("soft");
    Arguments::const_iterator to=args.upper_bound("soft");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      soft.addSpecifier(spec);
      lsoft=1;
    }
    //  get al2
    double lam=0, alj=0, nkt=0;
    iter=args.lower_bound("softpar");
    if(iter!=args.upper_bound("softpar")){
      lam=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("softpar")){
      alj=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("softpar"))
      nkt=atof(iter->second.c_str());
    else if(lsoft)
      throw gromos::Exception("Ener", 
	 "soft atoms indicated, but not all parameters defined.\n");
    
    en.setSoft(soft, lam, alj, nkt);
  }
 
  // define input coordinate
  InG96 ic;
  
  
  // print titles
  cout << "# Time"
       << "              covalent"
       << "            non-bonded"
       << "                 Total"
       << endl;

  // declare some variables for averaging
  int num_frames=0;
  double cov=0.0;
  double nb=0.0;
  double tot=0.0;
  
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      pbc->gathergr();

      // calculate the energies
      en.calc();

      // print any ouput you like
      cout.precision(10);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << time
	   << setw(22) << en.cov()
           << setw(22) << en.nb()
           << setw(22) << en.tot()
	   << endl;

      //store some averages
      cov+=en.cov();
      nb+=en.nb();
      tot+=en.tot();
      
      time+=dt;
      num_frames++;
    }
  }
  // print out averages
  if(num_frames>1){
    cout.precision(10);
    cout.setf(ios::right, ios::adjustfield);
    cout << endl << "# ave."
         << setw(22) << cov/num_frames 
         << setw(22) << nb/num_frames
         << setw(22) << tot/num_frames
         << endl;
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







