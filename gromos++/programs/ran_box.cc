/**
 * @file ranbox.cc
 * @page programs Program Documentation
 *
 * @anchor ranbox
 * @section ranbox generating box
 * @author @ref dt
 * @date 29. 11. 2004
 *
 * creates a box by inserting molecules randomly.
 * 
 * arguments:
 * - topo topologies
 * - pbc [v,r,t,c] [gathermethod]
 * - pos coordinate files
 * - nsm number of molecules
 * - dens system density in kg m^-3 
 * - thresh minimum atom-atom distance
 * 
 *
 * Example:
 * @verbatim
 ranbox
 @topo ic4.top urea.top h2o.top
 @pos  ic4.gsf urea.gsf h2o.gsf
 @nsm  1       153      847
 @dens 1000
 #@thresh 0.19
 
 @endverbatim
 *
 * <hr>
 */

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

#include <unistd.h>

using namespace std;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;


const double pi = acos(-1.0);
const double fac_amu2kg = 1.66056;

int main(int argc, char **argv){
  
  char *knowns[] = {"topo", "pbc", "pos", "nsm", "dens", "thresh", "layer", "boxsize"};
  int nknowns = 8;
  
  string usage = argv[0];
  usage += "\n\t@topo     <topologies of single molecule for each molecule type: topo1 topo2 ...>\n";
  usage += "\t@pbc      <boundary type>\n";
  usage += "\t@pos      <coordinates of single molecule for each molecule type: pos1 pos2 ...>\n";
  usage += "\t@nsm      <number of molecules for each molecule type: nsm1 nsm2 ...>\n";
  usage += "\t@dens     <density of liquid (kg/m^3)>\n";
  usage += "\t@thresh   <threshold distance in overlap check; default: 0.20 nm>\n";
  usage += "\t@layer    <create molecules in layers (along z axis)>\n";
  usage += "\t@boxsize  <boxsize>\n";
  
  srand(time(NULL));
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    // reading input and setting some values
    if ( args.count("topo") != args.count("pos") ||  args.count("topo") != args.count("nsm") ) {
      throw gromos::Exception("ran_box", "Check the number of arguments for @topo, @pos and @nsm");
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
    
    args.check("pos",1);
    vector<string> insxs;
    iter=args.lower_bound("pos");
    while(iter!=args.upper_bound("pos")){
      insxs.push_back(iter->second.c_str());
      ++iter;
    }    
    
    double box = 0.0;
    double vtot = 0.0;
    double densit = 0.0;

    // read all topologies only to get the box length (via the mass)
    double weight=0; 
    for(unsigned int tcnt=0; tcnt<tops.size(); tcnt++) {
      InTopology it(tops[tcnt]);
      System smol(it.system());
      for(int i=0; i<smol.numMolecules();i++)
	for(int j=0; j< smol.mol(i).numAtoms();j++)
	  weight+=nsm[tcnt]*smol.mol(i).topology().atom(j).mass(); 
    }
    
    if (args.count("boxsize") > 0){
      if (args.count("dens") >=0)
	throw Arguments::Exception("don't specify boxsize and density!");
      
      std::istringstream is(args["boxsize"]);
      if (!(is >> box))
	throw Arguments::Exception("could not read boxsize");
      
      vtot = pow(box, 3);
      if (args["pbc"] == "t") vtot /= 2;
      densit = weight * 1.66056 / vtot;
    }
    else{
      args.check("dens",1);
      iter=args.lower_bound("dens");
      densit=atof(iter->second.c_str());

      vtot=(weight*1.66056)/densit;
      // we need the volume, correct for truncated octahedron!!
      if(args["pbc"] == "t") vtot*=2;
      box=pow(vtot,1.0/3.0);
    }
    
    iter=args.lower_bound("thresh");
    double thresh = (iter!=args.upper_bound("thresh")) ? thresh=atof(iter->second.c_str()) : 0.20; 
    thresh *=thresh;
    
    bool layer = false;
    if (args.count("layer") >= 0) layer = true;
    
    std::cerr << "creating molecules in layers" << std::endl;
    
    // getting a reference vector in het midder van de box (pbc correction)
    Vec box_mid(box/2.0, box/2.0, box/2.0);
    
    // printing the box size
    cerr << "Cell volume: " << vtot << endl
         << "Total mass: " << weight * fac_amu2kg << endl
	 << "Input density: " << densit << endl
         << "Cell length: " << box << endl 
	 << "PBC: " << args["pbc"] << endl;
    
    // cerr << "Now, sleeping for 5 seconds... " 
    // << "there is still time for a ctrl-C!" << endl;
    // sleep(5);  

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
    
    // loop over the number of topologies.
    for(unsigned int tcnt=0; tcnt<tops.size(); tcnt++) {
      
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
      for(unsigned int i=0; i<unsigned(nsm[tcnt]); i++){
	sys.addMolecule(smol.mol(0));
	
	// get three random numbers between the box dimensions
	// get also the rotation axes and angles
	Vec rpos;
      UGLY_GOTO:
	for(int d=0; d<3; d++){
	  int r=rand();

	  if (d == 2 && layer)
	    rpos[d] = tcnt * (box/tops.size()) + double(box/tops.size() * r)/double(RAND_MAX);
	  else
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
	    sys.mol(sys.numMolecules()-1).pos(j)[0] = smol.mol(0).pos(j)[0] + rpos[0];
	    sys.mol(sys.numMolecules()-1).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp -
						       smol.mol(0).pos(j)[2]*sinp) + rpos[1];
	    sys.mol(sys.numMolecules()-1).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp +
						       smol.mol(0).pos(j)[1]*sinp) + rpos[2];
	  } 
	else if (r_axis==2) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    sys.mol(sys.numMolecules()-1).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp +
						       smol.mol(0).pos(j)[2]*sinp) + rpos[0];
	    sys.mol(sys.numMolecules()-1).pos(j)[1] = smol.mol(0).pos(j)[1] + rpos[1];
	    sys.mol(sys.numMolecules()-1).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp -
						       smol.mol(0).pos(j)[0]*sinp) + rpos[2];
	  }
	else if (r_axis==3) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    sys.mol(sys.numMolecules()-1).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp -
						       smol.mol(0).pos(j)[1]*sinp) + rpos[0];
	    sys.mol(sys.numMolecules()-1).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp +
						       smol.mol(0).pos(j)[0]*sinp) + rpos[1]; 
	    sys.mol(sys.numMolecules()-1).pos(j)[2] =  smol.mol(0).pos(j)[2] + rpos[2];
	  }
	
	//checking overlap
	if (sys.numMolecules() != 1) {
	  for(int k=0; k<sys.numMolecules()-1; k++) {
	    for(int l=0; l<sys.mol(k).numAtoms(); l++) {
	      for(int m=0; m<smol.mol(0).numAtoms(); m++) {
		if ( (sys.mol(k).pos(l) -
		      (pbc->nearestImage(sys.mol(k).pos(l), sys.mol(sys.numMolecules()-1).pos(m), sys.box()))).abs2()
		     < thresh)
		  goto UGLY_GOTO;
	      }
	    }
	  }
	}
	cerr << (i+1) << " of " << nsm[tcnt] << " copies of molecule " << tcnt+1 
	     << " already in the box. (Total number of molecules = " << sys.numMolecules() << ")." << endl;
      }
      cerr << "Box now with: " << sys.numMolecules() << " molecules" << endl;  
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


