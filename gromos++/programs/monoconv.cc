// monoconv

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Vec.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "beta", "boxformat"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@beta <monoclinic angle>\n";
  usage += "\t@boxformat <format of the box: box96, triclinicbox, genbox>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);


    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    InG96 ic;
    OutG96 oc(cout);
    ostringstream os;
    os << "Monoconv converted box block in\n";
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){
      os << "  " << iter->second << endl;
    }
    oc.writeTitle(os.str());
    oc.select("ALL");
    
    bool monoclinic=false;
    double beta=M_PI/2;
    if(args.count("beta")>0){
      monoclinic=true;
      beta=atof(args["beta"].c_str())/180.0*M_PI;
    }
    
    const double cosb=cos(beta);
    const double sinb=sin(beta);

    gcore::Box::boxformat_enum bf;
    if(args.count("boxformat")>0){
      if(args["boxformat"] == "box96") bf=gcore::Box::box96;
      else if(args["boxformat"] == "triclinicbox") bf=gcore::Box::triclinicbox;
      else if(args["boxformat"] == "genbox") bf=gcore::Box::genbox;
      else
	throw gromos::Exception("monoconv", "Boxformat not known");
    }
    
    // loop over all trajectories
    for(Arguments::const_iterator iter=args.lower_bound("traj");
      iter!=args.upper_bound("traj"); ++iter){

      ic.open(iter->second);
      ic.select("ALL");
      
      // loop over all frames
      while(!ic.eof()){
        ic >> sys;

	if(monoclinic){
	  sys.box().K()=gmath::Vec(sys.box()[0],0.0,0.0);
	  sys.box().L()=gmath::Vec(0.0,sys.box()[1],0.0);
	  sys.box().M()=gmath::Vec(sys.box()[2]*cosb, 0.0, sys.box()[2]*sinb);
	  sys.box().setNtb(gcore::Box::triclinic);
	}
	sys.box().boxformat()=bf;
	oc << sys;
	
      }    
    
      ic.close();
    }
    oc.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
