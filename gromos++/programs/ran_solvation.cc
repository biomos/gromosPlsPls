// written by gee, chris and d
// a terrible mixture of probox and ran_box.
// "cleaned" by d, whatever that means...
//
//
// - protein placed in centre of box
// - solvent molecules placed with COG at random
//   positions within box at given density
// - threshold checks for overlap
//

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InTopology.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gmath/physics.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include <vector>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;
using namespace gmath;
using namespace bound;
using namespace utils;

using namespace std;
using namespace bound;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace args;

// setting some constants
const double pi = acos(-1.0);
const double fac_amu2kg = 1.66056;

// defining some functions
void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as);
double calc_max_size(System &sys, AtomSpecifier &as, int dim);
double calc_mol_radius(System &sys);
double calc_mass_solv(const vector<string> &tops);


int main(int argc, char **argv){
  
  char *knowns[] = { "topo_pro", "insx_pro", "sev", "topo_solv", "insx_solv", "densit",
		     "molf_solv", "pbc", "min_wall", "boxsize", "thresh_pro", "thresh_solv"};
  
  int nknowns = 12;
  
  string usage = argv[0];
  usage += "\n";
  usage += "\t@topo_pro    <topology of main solute (protein)>\n";
  usage += "\t@insx_pro    <co-ordinates of main solute (protein)>\n";
  usage += "\t@sev         <solvent-excluded volume of main solute (protein)>\n";
  usage += "\n";
  usage += "\t@topo_solv   <list of (single molecule) topologies of solvents      : topo_solv1 topo_solv2 ...>\n";
  usage += "\t@insx_solv   <list of (single molecule) coordinate files of solvents: insx_solv1 insx_solv2 ...>\n";
  usage += "\t@molf_solv   <mole fraction of each solvent                         : molf_solv1 molf_solv2 ...>\n";
  usage += "\t@densit      <mass density of solvent mixture (kg/m^3)>\n";
  usage += "\n";
  usage += "\t@pbc         <boundary type>\n";
  usage += "\t@min_wall    <minimum distance from protein to box-face>\n";
  usage += "\t@boxsize     <length of box-edge>\n";
  usage += "\n";
  usage += "\t@thresh_pro  <threshold distance in overlap check (protein - solvent) ; default: 0.4 nm>\n";
  usage += "\t@thresh_solv <threshold distance in overlap check (solvent - solvent) ; default: 0.23 nm>\n";

  srand(time(NULL));  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);    
    Arguments::const_iterator iter ;
    
    // topo_pro:
    args.check("topo_pro",1);
    string topo_pro;
    iter=args.lower_bound("topo_pro");
    topo_pro = iter->second.c_str();
    // insx_pro:
    args.check("insx_pro",1);
    string insx_pro;
    iter=args.lower_bound("insx_pro");
    insx_pro = iter->second.c_str();
    // sev:
    args.check("sev",1);
    double sev;
    iter=args.lower_bound("sev");
    sev = atof(iter->second.c_str());

    // consistency check of solvent specification:
    if ( args.count("topo_solv") != args.count("insx_solv") ) {
      throw gromos::Exception("ran_solvation", "Check the number of arguments for @topo_solv, @insx_solv");
    }

    // topo_solv:
    args.check("topo_solv",1);
    vector<string> tops;
    iter=args.lower_bound("topo_solv");
    while(iter!=args.upper_bound("topo_solv")){
      tops.push_back(iter->second.c_str());
      ++iter;
    }
    // insx_solv:
    args.check("insx_solv",1);
    vector<string> insxs;
    iter=args.lower_bound("insx_solv");
    while(iter!=args.upper_bound("insx_solv")){
      insxs.push_back(iter->second.c_str());
      ++iter;
    }
    // mol_frac:
    args.check("molf_solv",1);
    vector<double> molfsv;
    iter=args.lower_bound("molf_solv");
    while(iter!=args.upper_bound("molf_solv")){
      molfsv.push_back(atof(iter->second.c_str()));
      ++iter;
    }
    // densit:    
    args.check("densit",1);
    iter=args.lower_bound("densit");
    double densit=atof(iter->second.c_str());
    //
    //
    // min_wall:    
    // read the minimum solute to wall distances.
    // three possibilities:
    // 1. nothing specified: the user will specify the box size
    //    using the @boxsize flag
    // 2. one number specified: the box will be cubic (pbc: r) or a
    //    trunc. oct. (pbc: t) with the size defined by maximum solute
    //    atom-atom distance plus twice this number
    // 3. three numbers specified: the box will be rectangular. The
    //    solute is rotated and the specified numbers are added to the
    //    maximum distances in the resulting dimensions
    vector<double> min_wall;
    int calculate_dimensions=0;
    {
      Arguments::const_iterator iter=args.lower_bound("min_wall"), 
	to=args.upper_bound("min_wall");
      while(iter!=to){
	min_wall.push_back(atof(iter->second.c_str()));
	++iter;
      }
      if(min_wall.size()) calculate_dimensions=1;
    }    
    // read the boxsize
    // three posibilities:
    // 1. nothing specified: the box size is determined from solute 
    //    geometry and min_wall specifications
    // 2. one number specified: the box will be cubic (pbc: r) or a trunc.
    //    oct. (pbc:t)
    // 3. three numbers specified: the box will be rectangular
    vector<double> boxsize;
    {
      Arguments::const_iterator iter=args.lower_bound("boxsize"),
	to=args.upper_bound("boxsize");
      while(iter!=to){
	boxsize.push_back(atof(iter->second.c_str()));
	++iter;
      }
      if(boxsize.size() == 1 || boxsize.size() == 2) 
	boxsize.resize(3, boxsize[0]);
    }
    if(boxsize.size() && min_wall.size())
      throw(gromos::Exception("ran_solvation", 
			      "Cannot specify both boxsize and min_wall."));
    // pbc:
    InTopology it(args["topo_pro"]);
    System solu(it.system());
    Boundary *pbc = BoundaryParser::boundary(solu, args);
    int truncoct=0, rectbox=0, cubic=0;
    double size_corr=1.0;
    double fac_vol=1.0;
    switch(pbc->type()){
    case('t'):
	truncoct=1;
	rectbox=0;
	size_corr = 2.0*sqrt(3.0)/3.0;
        fac_vol=0.5;
	if(min_wall.size()>1)
	  throw(gromos::Exception("ran_solvation", 
				  "For truncated octahedral boxes you can only specify one number for @min_wall"));
	if(min_wall.size()==0)
	  if(boxsize[0]!=boxsize[1] || boxsize[0]!=boxsize[2])
	    throw(gromos::Exception("ran_solvation", 
				    "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));
	break;
    case('r'):
      truncoct=0;
	rectbox=1;
        fac_vol=1.0;
	if(min_wall.size()==1) cubic=1;
	else if (min_wall.size()==0 &&  
		 (boxsize[0]==boxsize[1] && boxsize[0]==boxsize[2])) 
	  cubic=1;
	break;
    case('v'):
      throw(gromos::Exception("ran_solvation", 
			      "Why are you running this program if @pbc is vacuum?"));
      break;
    }

    // reading the threshold values
    //
    // thresh_pro    
    iter=args.lower_bound("thresh_pro");
    double thresh_pro = (iter!=args.upper_bound("thresh_pro")) ? thresh_pro=atof(iter->second.c_str()) 
      : 0.4; 
    thresh_pro *=thresh_pro;
    // thresh_solv
    iter=args.lower_bound("thresh_solv");
    double thresh_solv = (iter!=args.upper_bound("thresh_solv")) ? thresh_solv=atof(iter->second.c_str()) 
      : 0.23; 
    thresh_solv *=thresh_solv;
    
    // Compute vol_box:
    //  
    InG96 ic;
    ic.open(args["insx_pro"]);
    ic >> solu;
    ic.close();

    // shift molecule to cog
    double radius_pro = calc_mol_radius(solu);    
    Vec shiftcog=fit::PositionUtils::shiftToCog(&solu);   

    // there is only need to rotate if minwall has three elements
    // store the maximum dimensions in max_dim
    vector<double> max_dim;
    AtomSpecifier as(solu);
    if(min_wall.size()==3){
      rotate_solute(solu, max_dim, as);
      // calculate the box dimensions
      for(int i=0; i<3; i++)
	boxsize.push_back(max_dim[i]+2*min_wall[i]);
    }
    else if(boxsize.size()==0){
      max_dim.push_back(calc_max_size(solu, as, 3));
      // calculate box dimensions
      for(int i=0; i<3; i++)
	boxsize.push_back(size_corr*(max_dim[0] + 2*min_wall[0]));
    }
    
    for(int i=0;i<3;i++){ solu.box()[i] = boxsize[i]; }
    Vec box_mid(solu.box()[0]/2.0, solu.box()[1]/2.0, solu.box()[2]/2.0);
    PositionUtils::translate(&solu, box_mid);
    double vol_cell = fac_vol * boxsize[0] * boxsize[1] * boxsize[2] ;
    
    // Compute vol_solv:
    double vol_solv = vol_cell - sev;    

    // computing the masses of the solvents and
    vector<double> mass_solv;
    for(unsigned int i=0; i<tops.size(); i++)
      mass_solv.push_back(calc_mass_solv(tops));
    
    // Computing the number of molecules of each solvent required to match densit:
    double tmp=0;
    for(unsigned int i=0; i<tops.size(); i++)
      tmp+=molfsv[i]*mass_solv[i];
    vector<double> nr_solv; // number of molecules for each type of slv
    vector<double> tot_mass_solv;
    for(unsigned int i=0; i<tops.size(); i++) {
      nr_solv.push_back(int(densit*vol_solv*molfsv[i]/tmp));
      tot_mass_solv.push_back(mass_solv[i]* nr_solv[i]) ;
    }
    
    //
    // report computational specs:
    // 
    cerr << "Cell volume     :" << vol_cell  << endl;
    cerr << "SEV of protein  :" << sev       << endl;
    cerr << "Solvent volume  :" << vol_solv  << endl;
    cerr << "Solvent density :" << densit    << endl;
    for(unsigned int i=0; i<tops.size(); i++) {
      cerr << "# solvent " << i+1 << ": "  << nr_solv[i] << endl;
      cerr << "Total mass of solvent " << i+1 << ": " << tot_mass_solv[i] << endl;
    }
    
    //
    // MAIN PROCESS:
    //
    // loop over the number of solvent topologies. i has to be created here. inthebox counts
    // the number of molecules of tcnt topology already in the box
    for(unsigned int tcnt=0, i=0, inthebox=1; tcnt<tops.size(); tcnt++) {
      
      //read topologies again
      InTopology it(tops[tcnt]);
      System smol(it.system());

      // read single molecule coordinates...
      InG96 ic;
      ic.open(insxs[tcnt]);
      ic >> smol;
      ic.close();
      
      // single molecule coordinates are translated to reference frame 
      //  with origin at centre of mass
      fit::PositionUtils::shiftToCog(&smol);   
      
      //loop over the number of desired mol
      double nr_mol=nr_solv[tcnt];      
      for(i=0; i<nr_mol; i++){
	solu.addMolecule(smol.mol(0));
      UGLY_GOTO:
	
	// trial position of solvent molecule
	//
	Vec rpos;                       
	for(int d=0; d<3; d++){
	  int r=rand();
	  rpos[d]=double(solu.box()[0]*r)/double(RAND_MAX);
	}
	// if outside box then discard
        Vec rpos2=pbc->nearestImage(box_mid, rpos, solu.box());
        if(! (rpos2 == rpos)){ goto UGLY_GOTO;	}
	
	// trial orientation of solvent molecule
	//
	int r=rand();
	int r_axis = int(double(3.0*r)/(double(RAND_MAX)+1))+1;  // i.e. [1,2,3]
	double phi = 2.0*pi*double(r)/(double(RAND_MAX)+1);
	double cosp = cos(phi), sinp = sin(phi);
	
	// rotate and put the molecule in the box      
	if (r_axis==1)
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    solu.mol(i+inthebox).pos(j)[0] = smol.mol(0).pos(j)[0] + rpos[0];
	    solu.mol(i+inthebox).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp -
					      smol.mol(0).pos(j)[2]*sinp) + rpos[1];
	    solu.mol(i+inthebox).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp +
					      smol.mol(0).pos(j)[1]*sinp) + rpos[2];
	  } 
	else if (r_axis==2) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    solu.mol(i+inthebox).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp +
					     smol.mol(0).pos(j)[2]*sinp) + rpos[0];
	    solu.mol(i+inthebox).pos(j)[1] = smol.mol(0).pos(j)[1] + rpos[1];
	    solu.mol(i+inthebox).pos(j)[2] = (smol.mol(0).pos(j)[2]*cosp -
					     smol.mol(0).pos(j)[0]*sinp) + rpos[2];
	  }
	else if (r_axis==3) 
	  for(int j=0;j<smol.mol(0).numAtoms();j++) {
	    solu.mol(i+inthebox).pos(j)[0] = (smol.mol(0).pos(j)[0]*cosp -
					      smol.mol(0).pos(j)[1]*sinp) + rpos[0];
	    solu.mol(i+inthebox).pos(j)[1] = (smol.mol(0).pos(j)[1]*cosp +
					      smol.mol(0).pos(j)[0]*sinp) + rpos[1]; 
	    solu.mol(i+inthebox).pos(j)[2] =  smol.mol(0).pos(j)[2] + rpos[2];
	  }
	
	//checking overlap:
        // first with main solute (protein)
	radius_pro += thresh_pro;
        radius_pro *= radius_pro;      
	if( (box_mid - rpos).abs2() <= radius_pro){
	  for(int l=0; l<solu.mol(0).numAtoms(); l++) {
	    for(int m=0; m<smol.mol(0).numAtoms(); m++) {
	      // no pbc check because protein certainly within box
	      if ((solu.mol(0).pos(l)-solu.mol(i+inthebox).pos(m)).abs2()<thresh_pro){		
		goto UGLY_GOTO;
	      }
	      
	    }
	  }
	}//if rpos in sphere around protein
	
	// then with all other solvent molecules
	for(int k=1; k<solu.numMolecules()-1; k++) {
	  for(int l=0; l<solu.mol(k).numAtoms(); l++) {
	    for(int m=0; m<smol.mol(0).numAtoms(); m++) {
	      if ( (solu.mol(k).pos(l) -
		    (pbc->nearestImage(solu.mol(k).pos(l), solu.mol(i+inthebox).pos(m), solu.box()))).abs2()
		   < thresh_solv){
		goto UGLY_GOTO;
	      }
	    }
	  }
	}
	
	cerr << (i+1) << " of " << nr_mol << " copies of solvent molecule of type " << tcnt+1 
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
    oc << solu;
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


// function bodies

void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as)
{
  max.resize(3);
  
  // this will be a two stage function.
  // First calculate the maximum distance between any two solute atoms
  max[2] =calc_max_size(sys, as, 3);
  
  // rotate the solute such that the atoms in as are along the z-axis
  Vec v=sys.mol(as.mol(0)).pos(as.atom(0)) 
    - sys.mol(as.mol(1)).pos(as.atom(1));
  double r = v.abs();
  double r_yz = sqrt(v[1]*v[1] + v[2]*v[2]);  
  
  // the rototation matrix is the product of two rotations
  // 1. around the x-axis by theta
  //    with sin(theta) = v[2]/r_yz; cos(theta) = align[2]/r_yz
  // 2. around the y-axis by phi
  // with sin(phi) = align[0]/r; cos(phi) = r_yz / r
  Matrix rot1(Vec( r_yz / r         ,  0         , v[0]/r ),
	      Vec(-v[0]*v[1]/r/r_yz ,  v[2]/r_yz , v[1]/r ),
	      Vec(-v[0]*v[2]/r/r_yz , -v[1]/r_yz , v[2]/r ));
  
  for(int m=0; m<sys.numMolecules(); m++)
    for(int a=0; a<sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot1*sys.mol(m).pos(a);
  // and take along any solvent
  for(int i=0; i<sys.sol(0).numPos(); i++){
    sys.sol(0).pos(i) = rot1*sys.sol(0).pos(i);
  }
  
  // calculate the maximum distance in the x-y-plane
  max[1] =calc_max_size(sys, as, 2);
  
  // rotate the solute around the z-axis, such that the atoms in as are
  // along the y-axis, this is done by a rotation around psi with
  // sin(psi) = x/r_xy; cos(psi) = y/r_xy;
  v = sys.mol(as.mol(2)).pos(as.atom(2)) - sys.mol(as.mol(3)).pos(as.atom(3));
  double r_xy= sqrt(v[0]*v[0] + v[1] * v[1]);
  
  Matrix rot2(Vec( +v[1]/r_xy ,  v[0]/r_xy , 0),
	      Vec( -v[0]/r_xy ,  v[1]/r_xy , 0),
	      Vec( 0         ,  0         , 1));
  
  for(int m=0; m<sys.numMolecules(); m++)
    for(int a=0; a<sys.mol(m).numAtoms(); a++)
      sys.mol(m).pos(a) = rot2*sys.mol(m).pos(a);
  
  // and take along any solvent
  for(int i=0; i<sys.sol(0).numPos(); i++){
    sys.sol(0).pos(i) = rot2*sys.sol(0).pos(i);
  }
  
  // finally we calculate the maximum distance in the x-direction
  max[0] =calc_max_size(sys, as, 1);
}


double calc_max_size(System &sys, AtomSpecifier &as, int dim)
{
  // calculate the longest distance between solute atoms considering
  // the first dim dimensions.
  double max2=0.0;
  double d2=0;
  int max_m1=0, max_m2=0, max_a1=0, max_a2=0;
  
  for(int m1=0; m1 < sys.numMolecules(); m1++){
    for(int a1=0; a1 < sys.mol(m1).numAtoms(); a1++){
      Vec current=sys.mol(m1).pos(a1);
      for(int m2=m1; m2 < sys.numMolecules(); m2++){
	int start=0;
	if(m1==m2) start=a1;
	for(int a2=start; a2 < sys.mol(m2).numAtoms(); a2++){
	  d2=0.0;
	  for(int i=0; i<dim; i++){
	    d2 += (current[i] - sys.mol(m2).pos(a2)[i]) *
	      (current[i] - sys.mol(m2).pos(a2)[i]);
	  }
	  if(d2>max2){
	    max2=d2;
	    max_m1=m1;
	    max_m2=m2;
	    max_a1=a1;
	    max_a2=a2;
	  }
	}
      }
    }
  }
  as.addAtomStrict(max_m1, max_a1);
  as.addAtomStrict(max_m2, max_a2);
  
  return sqrt(max2);
}


double calc_mol_radius(System &sys)
{
  // calculate the radius of the solute molecule
  double d2=0;
  
  Vec cog(0.0,0.0,0.0);
  int counter=0;  
  for(int m2=0; m2 < sys.numMolecules(); m2++){

    for(int a2=0; a2 < sys.mol(m2).numAtoms(); a2++){
      cog+=sys.mol(m2).pos(a2); 
    }
    counter += sys.mol(m2).numAtoms() ;    
  }
  cog/=counter;
  for(int m2=0; m2 < sys.numMolecules(); m2++){
    
    for(int a2=0; a2 < sys.mol(m2).numAtoms(); a2++){
      double t=(cog-sys.mol(m2).pos(a2)).abs2();
      
      if(d2< t) d2=t;
    }
  }
  return sqrt(d2);
}


double calc_mass_solv(const vector<string> &tops)
{
  // calculates the total mass of type i solvent
  InTopology itA(tops[0]);
  System solv_A(itA.system());

  double mass_A=0;
  if(solv_A.numMolecules() != 1){     
    throw(gromos::Exception("ran_solvation", 
			    "Topology of first solvent has more than one molecule in SOLUTEATOMS block."));
  }
  else{
    for(int j=0; j< solv_A.mol(0).numAtoms();j++){
      mass_A+=solv_A.mol(0).topology().atom(j).mass(); 
    }
    mass_A*=fac_amu2kg;
  } 
  return mass_A;
}