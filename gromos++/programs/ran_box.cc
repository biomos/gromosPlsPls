// ran_box.cc -- should generate a random box -- d 

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>

using namespace std;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;


const double pi = acos(-1.0);


int main(int argc, char **argv){
  
  char *knowns[] = {"topo", "pbc", "insx", "nsm", "densit", "thresh"};
  int nknowns = 6;
  
  string usage = argv[0];
  usage += "\n\t@topo   <topologies of single molecule for each molecule type: topo1 topo2 ...>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@insx   <coordinates of single molecule for each molecule type: insx1 insx2 ...>\n";
  usage += "\t@nsm    <number of molecules for each molecule type: nsm1 nsm2 ...>\n";
  usage += "\t@densit <density of liquid (kg/m^3)>\n";
  usage += "\t@thresh <threshold distance in overlap check; default: 0.22 nm>";
  
  srand(time(NULL));
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    // reading input and setting some values
    if ( args.count("topo") != args.count("insx") ||  args.count("topo") != args.count("nsm") ) {
      throw gromos::Exception("ran_box", "Check the number of arguments for @topo, @insx and @nsm");
    }
    
    args.check("nsm",1);
    vector<int> nsm;
    Arguments::const_iterator iter=args.lower_bound("nsm");
    while(iter!=args.upper_bound("nsm")){
      nsm.push_back(atoi(iter->second.c_str()));
      ++iter;
    }

    args.check("topo",1);
    vector<string> tops;
    iter=args.lower_bound("topo");
    while(iter!=args.upper_bound("topo")){
      tops.push_back(iter->second.c_str());
      ++iter;
    }
    
    args.check("insx",1);
    vector<string> insxs;
    iter=args.lower_bound("insx");
    while(iter!=args.upper_bound("insxs")){
      insxs.push_back(iter->second.c_str());
      ++iter;
    }    
    
    args.check("densit",1);
    iter=args.lower_bound("densit");
    double densit=atof(iter->second.c_str());
    
    iter=args.lower_bound("thresh");
    double thresh = (iter!=args.upper_bound("thresh")) ? thresh=atof(iter->second.c_str()) : 0.22; 
    thresh *=thresh;
    
    // read all topologies only to get the box length (via the mass)
    // daan has suggested this to make the program even unglier than it is!!!!
    double weight=0; 
    for(unsigned int tcnt=0; tcnt<tops.size(); tcnt++) {
      InTopology it(tops[tcnt]);
      System smol(it.system());
      for(int i=0; i<smol.numMolecules();i++)
	for(int j=0; j< smol.mol(i).numAtoms();j++)
	  weight+=nsm[tcnt]*smol.mol(i).topology().atom(j).mass(); 
    }
    double vtot=(weight*1.66056)/densit;
    // we need the volume, correct for truncated octahedron!!
    if(args["pbc"] == "t") vtot*=2;
    double box=pow(vtot,1.0/3.0);
    
    // getting a reference vector in het midder van de box (pbc correction)
    Vec box_mid(box/2.0, box/2.0, box/2.0);
    
    // printing the box size
    cerr << "Cell volume: " << vtot << endl;
    cerr << "Length: " << box << endl;
    // sleep(2);

    // now we do the whole thing
    // new system and parse the box sizes
    System sys;
    for(int i=0;i<3;i++){
      sys.box()[i] = box;
    }

    // parse boundary conditions
    Boundary *pbc; 
    if(args["pbc"] == "t")
      pbc = new TruncOct(&sys);
    else 
      pbc = new RectBox(&sys);
    
    // loop over the number of topologies. i has to be created here. inthebox counts
    // the number of molecules of tcnt topology already in the box
    for(unsigned int tcnt=0, i=0, inthebox=0; tcnt<tops.size(); tcnt++) {
      
      //read topologies again
      InTopology it(tops[tcnt]);
      System smol(it.system());
      
      // read single molecule coordinates...
      InG96 ic;
      ic.open(insxs[tcnt]);
      ic >> smol;
      ic.close();
     
      // single molecule coordinates are translated to reference frame 
      //  with origin at cog
      fit::PositionUtils::shiftToCog(&smol);

      //loop over the number of desired mol    
      for(i=0; i<unsigned(nsm[tcnt]); i++){
	sys.addMolecule(smol.mol(0));
	
	// get three random numbers between the box dimensions
	// get also the rotation axes and angles
	Vec rpos;
      UGLY_GOTO:
	for(int d=0; d<3; d++){
	  int r=rand();
	  rpos[d]=double(box*r)/double(RAND_MAX);
	}
	// correcting rpos for pbc
        Vec rpos2=pbc->nearestImage(box_mid, rpos, sys.box());
	if(! (rpos2 == rpos)){ goto UGLY_GOTO;  }

	int r=rand();
	int r_axis = int(double(3.0*r)/(double(RAND_MAX)+1))+1;
	double phi = 2.0*pi*double(r)/(double(RAND_MAX)+1);
	double cosp = cos(phi), sinp = sin(phi);
	
	// rotate and put the molecule in the box      
	if (r_axis==1)
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    sys.mol(i+inthebox).pos(j)[0] = smol.mol(0).pos(j)[0] + rpos[0];
	    sys.mol(i+inthebox).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp -
					     smol.mol(0).pos(j)[2]*sinp) + rpos[1];
	    sys.mol(i+inthebox).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp +
					     smol.mol(0).pos(j)[1]*sinp) + rpos[2];
	  } 
	else if (r_axis==2) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    sys.mol(i+inthebox).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp +
					     smol.mol(0).pos(j)[2]*sinp) + rpos[0];
	    sys.mol(i+inthebox).pos(j)[1] = smol.mol(0).pos(j)[1] + rpos[1];
	    sys.mol(i+inthebox).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp -
					     smol.mol(0).pos(j)[0]*sinp) + rpos[2];
	  }
	else if (r_axis==3) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    sys.mol(i+inthebox).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp -
					     smol.mol(0).pos(j)[1]*sinp) + rpos[0];
	    sys.mol(i+inthebox).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp +
					     smol.mol(0).pos(j)[0]*sinp) + rpos[1]; 
	    sys.mol(i+inthebox).pos(j)[2] =  smol.mol(0).pos(j)[2] + rpos[2];
	  }
	
	//checking overlap
	if ((i+inthebox)!=0) {
	  for(int k=0; k<sys.numMolecules()-1; k++) {
	    for(int l=0; l<sys.mol(k).numAtoms(); l++) {
	      for(int m=0; m<smol.mol(0).numAtoms(); m++) {
		// next line is for no pbc check!
		// if ((sys.mol(k).pos(l)-sys.mol(i+inthebox).pos(m)).abs2()<thresh)
		if ( (sys.mol(k).pos(l) -
		      (pbc->nearestImage(sys.mol(k).pos(l), sys.mol(i+inthebox).pos(m), sys.box()))).abs2()
		     < thresh)
		  goto UGLY_GOTO;
	      }
	    }
	  }
	}
	cerr << (i+1) << " of " << nsm[tcnt] << " copies of molecule " << tcnt+1 
	     << " already in the box. (Total number of molecules = " << i+1+inthebox << ")." << endl;
      }
      inthebox+=i;
      cerr << "Box now with: " << inthebox << " molecules" << endl;  
    }
    
    // Print the new set to cout
    OutG96S oc;
    ostringstream os;
    oc.open(cout);
    oc.writeTitle(string(os.str()));
    oc << sys;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


