// probox: put a solute in a box of solvent

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


void rotate_solute(System &sys, vector<double> &max, AtomSpecifier &as);
double calc_max_size(System &sys, AtomSpecifier &as, int dim);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "solute", "solvent", "minwall", "minsol", 
		    "boxsize"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo    <topology>\n";
  usage += "\t@pbc     <periodic boundary conditions>\n";
  usage += "\t@solute  <solute coordinates>\n";
  usage += "\t@solvent <solvent coordinates>\n";
  usage += "\t@minwall <minimum solute to wall distance>\n";
  usage += "\t@minsol  <minimum solvent-solute distance (default 0.23 nm)>\n";
  usage += "\t@boxsize <specified boxsize>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology make 2 systems 
    InTopology it(args["topo"]);
    System solu(it.system());
    System solv;
    solv.addSolvent(solu.sol(0));
    
    // set the definition of a hydrogen, for the heavy atom criterion
    for(int m=0; m<solu.numMolecules(); m++){
      solu.mol(m).topology().setHmass(1.008);
    }
    solu.sol(0).topology().setHmass(1.008);
    
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
    vector<double> minwall;
    int calculate_dimensions=0;
    {
      Arguments::const_iterator iter=args.lower_bound("minwall"), 
	to=args.upper_bound("minwall");
      while(iter!=to){
	minwall.push_back(atof(iter->second.c_str()));
	++iter;
      }
      if(minwall.size()) calculate_dimensions=1;
    }
    
    // read the boxsize
    // three posibilities:
    // 1. nothing specified: the box size is determined from solute 
    //    geometry and minwall specifications
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
    if(boxsize.size() && minwall.size())
      throw(gromos::Exception("probox", 
			      "Cannot specify both boxsize and minwall."));

    // read the minimum solvent-solute distance
    double minsol=0.23;
    if(args.count("minsol")>0) minsol=atof(args["minsol"].c_str());
    double minsol2 = minsol * minsol;
    
    // read coordinates into the systems
    InG96 ic;
    ic.open(args["solute"]);
    // we also read in any solvent that is already in the file
    ic.select("ALL");
    ic >> solu;
    ic.close();
   
    ic.open(args["solvent"]);
    ic.select("SOLVENT");
    ic >> solv;
    ic.close();
    
    int num_solv_atoms_per_box = solv.sol(0).numCoords();
    
    // determine box shape, calculate relevant things and check the 
    // input for consistency
    Boundary *pbc = BoundaryParser::boundary(solu, args);
    int truncoct=0, rectbox=0, cubic=0;
    double size_corr=1.0;
    switch(pbc->type()){
      case('t'):
	truncoct=1;
	rectbox=0;
	size_corr = 2.0*sqrt(3.0)/3.0;
	
	if(minwall.size()>1)
	  throw(gromos::Exception("probox", 
	     "For truncated octahedral boxes you can only specify one number for @minwall"));
	if(minwall.size()==0)
	  if(boxsize[0]!=boxsize[1] || boxsize[0]!=boxsize[2])
	    throw(gromos::Exception("probox", 
               "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));
	break;
      case('r'):
	truncoct=0;
	rectbox=1;
	if(minwall.size()==1) cubic=1;
	else if (minwall.size()==0 &&  
		 (boxsize[0]==boxsize[1] && boxsize[0]==boxsize[2])) 
	  cubic=1;
	break;
      case('v'):
	throw(gromos::Exception("probox", 
           "Why are you running this program if @pbc is vacuum?"));
	break;
    }

    // get the number of original solvent atoms    
    int numSolventAtoms=solu.sol(0).numCoords();
      
    // move the solute to the centre of geometry
    Vec shiftcog=fit::PositionUtils::shiftToCog(&solu);
    // shift the solvent atoms in there as well
    for(int i=0; i< numSolventAtoms; i++){
      solu.sol(0).pos(i)+= shiftcog;
    }
    
    // calculate the solvent cog
    Vec solv_cog(0.0,0.0,0.0);
    for(int i=0; i<solv.sol(0).numCoords(); i++)
      solv_cog+=solv.sol(0).pos(i);
    solv_cog/=solv.sol(0).numCoords();
    // move to solvent cog
   for(int i=0; i<solv.sol(0).numCoords(); i++)
     solv.sol(0).pos(i) -= solv_cog;
   
    // there is only need to rotate if minwall has three elements
    // store the maximum dimensions in max_dim
    vector<double> max_dim;
    AtomSpecifier as(solu);
    if(minwall.size()==3){
      rotate_solute(solu, max_dim, as);
      // calculate the box dimensions
      for(int i=0; i<3; i++)
	boxsize.push_back(max_dim[i]+2*minwall[i]);
    }
    else if(boxsize.size()==0){
      max_dim.push_back(calc_max_size(solu, as, 3));
      // calculate box dimensions
      for(int i=0; i<3; i++)
	boxsize.push_back(size_corr*(max_dim[0] + 2*minwall[0]));
    }

    // set the solute box size, calculate how many solvent boxes we should
    // use; calculate where to move the initial box
    vector<int> needed_boxes;
    Vec move_solvent;
    for(int i=0; i<3; i++){
      solu.box()[i]=boxsize[i];
      needed_boxes.push_back(int(solu.box()[i]/solv.box()[i])+1);
      move_solvent[i]= -0.5*(needed_boxes[i]-1)*solv.box()[i];
    }

    // move the initial box so that after multiplication the complete box
    // is centered around the origin.
    for(int i=0; i< num_solv_atoms_per_box; i++)
      solv.sol(0).pos(i)+=move_solvent;
    
    // do the multiplications
    for(int ix=0; ix<needed_boxes[0]; ix++){
      
      for(int iy=0; iy<needed_boxes[1]; iy++){
	
	for(int iz=0; iz<needed_boxes[2]; iz++){
	  if(ix!=0 || iy!=0 || iz!=0){
	    Vec shift(ix*solv.box()[0], iy*solv.box()[1], iz*solv.box()[2]);
	    for(int atom = 0; atom < num_solv_atoms_per_box; atom++){
	      solv.sol(0).addCoord(solv.sol(0).pos(atom)+shift);
	    }
	  }
	}
      }
    }
    int num_atoms_per_solvent=solv.sol(0).topology().numAtoms();
    int num_solvent_molecules=solv.sol(0).numCoords() / num_atoms_per_solvent;
    

    // now we have to keep only those waters that are inside the box and 
    // far enough away from the solute
    // we look at the centre of geometry of the solvent molecule
    Vec o(0.0,0.0,0.0);
    double min_init = solu.box()[0] * solu.box()[0]
      + solu.box()[1] * solu.box()[1]
      + solu.box()[2] * solu.box()[2];
    
    for(int i=0; i< num_solvent_molecules; i++){
      // calculate the centre of geometry of this solvent
      Vec sol_i(0.0,0.0,0.0);
      for(int j=0; j< num_atoms_per_solvent; j++)
	sol_i+=solv.sol(0).pos(num_atoms_per_solvent*i+j);
      sol_i/= num_atoms_per_solvent;
      
      // are we inside the box
      Vec check=pbc->nearestImage(o, sol_i, solu.box());
      if(check[0]==sol_i[0] && 
	 check[1]==sol_i[1] && 
	 check[2]==sol_i[2]){
	// yes we are in the box
	// calculate the closest distance to any solute
	double min2 = min_init;
	for(int m=0; m< solu.numMolecules(); m++){
	  for(int a=0; a<solu.mol(m).numAtoms(); a++){
	    if(!solu.mol(m).topology().atom(a).isH() && 
	       (check - solu.mol(m).pos(a)).abs2() < min2)
	      min2 = (check - solu.mol(m).pos(a)).abs2();
	  }
	}
	// or original solvent
	for(int j=0; j<numSolventAtoms; j++){
	  if(!solu.sol(0).topology().atom(j%num_atoms_per_solvent).isH() &&
	     (check - solu.sol(0).pos(j)).abs() < min2)
	    min2 = (check - solu.sol(0).pos(j)).abs2();
	}
	
	if(min2>minsol2){
	  // yes! we keep this solvent 
	  for(int k=0; k< num_atoms_per_solvent; k++)
	    solu.sol(0).addCoord(solv.sol(0).pos(num_atoms_per_solvent*i+k));
	}
      }
    }

    ostringstream title;
    title << "Solvating " << args["solute"] << " in " << args["solvent"] 
	  << endl;
    title << "Box dimensions (";
    if(truncoct) title << "truncated octahedron";
    else if(cubic) title << "cubic";
    else title << "rectangular";
    title << ") were ";
    if(calculate_dimensions){
      title << "calculated from maximum" << endl;
      title << "solute atom-atom distance";
      if(max_dim.size()==1)
	title << " (not rotated):" << endl;
      else
	title << "s (after rotation):" << endl;
      for(unsigned int i=0; i<max_dim.size(); i++){
	int index = 2*(max_dim.size() - i - 1);
	title << "\t" << max_dim[i] << " between atoms "  
	      << as.mol(index)+1 << ":" << as.atom(index)+1 << " and " 
	      << as.mol(index+1)+1 << ":" << as.atom(index+1)+1 << endl;
      }
    }
    else 
      title <<"specified by user" << endl;
    if(numSolventAtoms){
      title << "System contained " << numSolventAtoms/num_atoms_per_solvent
	    << " solvent molecules" << endl;
    }
    title << "Added " << (solu.sol(0).numCoords()-numSolventAtoms)
      /num_atoms_per_solvent
	  << " solvent molecules";
    
    OutG96S oc(cout);      
    oc.select("ALL");
    oc.writeTitle(title.str());
    oc << solu;
    oc.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
    
  }
  return 0;
}

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
  for(int i=0; i<sys.sol(0).numCoords(); i++){
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
  for(int i=0; i<sys.sol(0).numCoords(); i++){
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
  as.addAtom(max_m1, max_a1);
  as.addAtom(max_m2, max_a2);
  
  return sqrt(max2);
}

