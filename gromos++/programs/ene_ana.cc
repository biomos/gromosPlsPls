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
#include <cmath>

using namespace args;
using namespace gio;
using namespace gcore;


void print(gmath::stat &p, string s, double time, double dt);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "time", "files", "gasener", "temp"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@time    <t and dt>\n";
  usage += "\t@files   <energy files>\n";
  usage += "\t@gasener <average gasphase energy>\n";
  usage += "\t@temp    <temperature>\n";
  

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

    // read in the gas phase energy
    double gasphase=0;
    {
      Arguments::const_iterator iter=args.lower_bound("gasener");
      if(iter!=args.upper_bound("gasener"))
	gasphase=atof(iter->second.c_str());
    }
    
    // and the temperature
    double temp=atof(args["temp"].c_str());
    
    // setup the energy arrays
    double ener[22], eneres[6], volprt[20];
    double engrp[4];

    // and the statistical information
    gmath::stat totene;
    gmath::stat totkin;
    gmath::stat totpot;
    gmath::stat pressu;
    gmath::stat densit;
    gmath::stat boxvol;
    gmath::stat hvap;
    
    // loop over the files
    Ginstream gin;
 
    Arguments::const_iterator iter=args.lower_bound("files"),
      to=args.upper_bound("files");
    for(;iter!=to; ++iter){
      gin.open((iter->second).c_str());
      string dum;
      int num;
      while(dum!="ENERGY") gin >> dum;
     
      while(!gin.eof()){
	//while(dum!="ENERGY") gin >> dum;
	for(int i=0; i<22; i++) gin >> ener[i];
	for(int i=0; i<6 ; i++) gin >> eneres[i];
	gin >> num;
	for(int i=0; i<num; i++)
	  for(int j=0; j<4; j++)
	    gin >> engrp[j];
	while(dum!="VOLUMEPRESSURE") gin >> dum;
	for(int i=0; i<20; i++) gin >> volprt[i];
        
	// fill the stat-blocks
	// !!!! NOTE gromos++ starts counting at 0 !!!!
        totene.addval(ener[0]);
	totkin.addval(ener[1]);
	totpot.addval(ener[8]);
	pressu.addval(volprt[11]*16.388453);
	densit.addval(mass*1.66056/volprt[7]);
	boxvol.addval(volprt[7]);
	hvap.addval(gasphase-ener[8]/sys.numMolecules()+BOLTZ*temp);
	
	
	while(dum!="ENERGY"&&!gin.eof()) gin >> dum;

      }
      gin.close();
      
    }
    //print out the statistical information
    cout << setw(10) << "property"
	 << setw(14) << "average"
	 << setw(14) << "rmsd"
	 << setw(14) << "error est."
	 << endl;
    print(totene, "totene", time, dt);
    print(totkin, "totkin", time, dt);
    print(totpot, "totpot", time, dt);
    print(pressu, "pressu", time, dt);
    print(densit, "densit", time, dt);
    print(boxvol, "boxvol", time, dt);
    print(hvap  , "hvap"  , time, dt);
      
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void print(gmath::stat &p, string s, double time, double dt)
{
  ostrstream os;
  os << s << ".out" << ends;
  ofstream fout(os.str());
  fout << "#"
       << setw(9) << "time"
       << setw(14) << s
       << endl;
  for(int i=0; i< p.n(); i++){
    fout << setw(10) << time+i*dt
	 << setw(14) << p.val(i)
	 << endl;
  }
  fout.close();
  // and print the averages etc to cout
  cout << setw(10) << s
       << setw(14) << p.ave()
       << setw(14) << p.rmsd()
       << setw(14) << p.ee()
       << endl;
}




