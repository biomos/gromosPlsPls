//Radius of Gyration

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include <vector>
#include <iomanip>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "mol", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@mol <number of molecules to consider>\n";
  usage += "\t@traj <trajectory files>\n";
  
 
    try{
    Arguments args(argc, argv, nknowns, knowns, usage);


  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());

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
  
  // set molecule number
  vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i< sys.numMolecules();++i)
        mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
          it!=args.upper_bound("mol");++it){
        if(atoi(it->second.c_str())> sys.numMolecules())
          throw Arguments::Exception(usage);
        mols.push_back(atoi(it->second.c_str())-1);
      }
   
      
    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
  InG96 ic;

    double totalMass=0;
    double totA=0;
    for (int i=0; i< (int) mols.size();++i){
     totA += sys.mol(mols[i]).numAtoms();
     for (int j=0; j < sys.mol(mols[i]).numAtoms(); j++) {      
      totalMass += sys.mol(mols[i]).topology().atom(j).mass();
    }
  }

    
    // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
      
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      (*pbc.*gathmethod)();
   
      //calculate cm, rgyr
  
  Vec cm;cm[0]=cm[1]=cm[2]=0;
  for (int i=0; i< (int) mols.size();++i){
    for (int j=0; j < sys.mol(mols[i]).numAtoms(); j++) {
      cm += sys.mol(mols[i]).pos(j) * sys.mol(mols[i]).topology().atom(j).mass();
    }
  }
      cm = (1.0/totalMass)*cm;


       double rg=0;   
       
      for (int i=0; i< (int) mols.size();++i){       
       for (int j=0; j<sys.mol(mols[i]).numAtoms();++j){ 
        Vec tmp =  (sys.mol(mols[i]).pos(j)-cm);	
	rg += tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];
       }              
      }
       rg = sqrt(rg/(totA));
    
    cout << setw(10) << time ;
    cout << setw(10) << rg << "\n";
       time += dt;
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
