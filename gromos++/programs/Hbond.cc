//Hydrogen bond analysis Hbond
//dedicated to wilfred: "you know, mika: proahb IS my favorite program"
//--mika

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/Neighbours.h"
#include "../src/utils/Hintra.h"
#include "../src/utils/Hinter.h"
#include "../src/utils/HintraNative.h"
#include "../src/utils/HinterNative.h"
#include "../src/utils/Hsolusolv.h"
#include "../src/utils/Hsolvsolv.h"


using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "type",  "moln", "ref", "Hbparas", "time", "traj"};
  int nknowns = 8;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@type   <1 2 3 4 5 6>\n"
    "\t\tChoices are:\n"
    "\t\t1: Intramolecular\n"
    "\t\t2: Intermolecular\n"
    "\t\t3: native Intramolecular\n"
    "\t\t4: native Intermolecular\n"
    "\t\t5: Solute <-> Solvent\n"
    "\t\t6: Solvent <-> Solvent (SLOW...)\n";
  usage += "\t@moln   <molecule numbers in topology 1..X>\n";
  usage += "\t@ref <reference coordinates for native H-bonds>\n";
  usage += "\t@Hbparas <distance [nm] and angle; default: 0.25, 135>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
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
  



// set molecule numbers
   int  m1=0, m2=1;
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     m1=atoi(iter->second.c_str())-1;
     ++iter;
   }
   if(iter!=args.upper_bound("moln"))
        m2=atoi(iter->second.c_str())-1;
   }

   //get type
   int type=0;
     {
   Arguments::const_iterator iter=args.lower_bound("type");
   if(iter!=args.upper_bound("type")){
     type=atoi(iter->second.c_str());
   }
    }
    
  // get the paras
  double maxdist=0.25, minangle = 135;
  {
    Arguments::const_iterator iter=args.lower_bound("Hbparas");
    if(iter!=args.upper_bound("Hbparas")){
      maxdist=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("Hbparas"))
        minangle=atof(iter->second.c_str());
  }


  //switch statement

  switch(type){
  case 1: Intramolecular(args, m1, time, dt, maxdist, minangle);
  break;
  case 2: Intermolecular(args, m1, m2, time, dt, maxdist, minangle);
  break;
  case 3: natIntramolecular(args, m1, time, dt, maxdist, minangle);
  break;
  case 4: natIntermolecular(args, m1, m2, time, dt, maxdist, minangle);
  break;
  case 5: SoluSolv(args, m1, time, dt, maxdist, minangle);
  break;
  case 6: SolvSolv(args, time, dt, maxdist, minangle);
  break;
  default:
  throw gromos::Exception("Hond", args["type"] + 
                                " unknown. Known types are 1, 2, 3, 4, 5 and 6");
  }


    
    }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
