/**
 * @file probox.cc
 * @page programs Program Documentation
 *
 * @anchor probox
 * @section probox generating box
 * @author @ref co
 * @date 30. 11. 2004
 *
 * This program has been renamed to @ref sim_box and will no longer be 
 * maintained under this name
 *
 * <hr>
 */

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
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
		    "boxsize", "gather", "rotate"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo    <topology>\n";
  usage += "\t@pbc     <periodic boundary conditions>\n";
  usage += "\t@solute  <solute coordinates>\n";
  usage += "\t@solvent <solvent coordinates>\n";
  usage += "\t[@minwall <minimum solute to wall distance>]\n";
  usage += "\t[@minsol  <minimum solvent-solute distance (default 0.23 nm)>]\n";
  usage += "\t[@boxsize <use specified boxsize (in solute coordinates)>]\n";
  usage += "\t[@gather  <gather solute>]\n";
  usage += "\t[@rotate  <rotate solute: biggest axis along z, second along y>]\n";


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
    // 3. three numbers specified: the box will be rectangular.
    //    the specified numbers are added to the
    //    maximum distances in x, y and z (after possibly rotating the solute)

    vector<double> minwall;
    {
      Arguments::const_iterator iter=args.lower_bound("minwall"), 
	to=args.upper_bound("minwall");
      while(iter!=to){
	minwall.push_back(atof(iter->second.c_str()));
	++iter;
      }
    }

    // read the minimum solvent-solute distance
    double minsol=0.23;
    if(args.count("minsol")>0) minsol=atof(args["minsol"].c_str());
    double minsol2 = minsol * minsol;

    // check for the boxsize flag
    // if it is given, the box from the solute coordinates is used
    bool boxsize = false;
    if (args.count("boxsize") != -1)
      boxsize = true;
    
    // read coordinates into the systems
    InG96 ic;
    ic.open(args["solute"]);
    // we also read in any solvent that is already in the file
    ic.select("ALL");
    ic >> solu;
    ic.close();

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(solu, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    if (args.count("gather") != -1){
      // should not be used if boxsize is specified
      // ('cause you probably changed the boxsize by hand)
      if (boxsize)
	throw gromos::Exception("probox",
				"don't gather if you specified a boxsize in the solute coordinates\n"
				"they should be gathered before!");
      if (!solu.hasBox)
	throw gromos::Exception("probox",
				"solute gathering requested, but solute does not contain a box block!");
      (*pbc.*gathmethod)();
    }

    vector<double> max_dist;
    AtomSpecifier as(solu);

    if (args.count("rotate") != -1){
      
      rotate_solute(solu, max_dist, as);
      
    }
    else{
      if (minwall.size() == 3)
	throw gromos::Exception("probox",
				"solute should be rotated in order to align largest extension with z axis");
      
      max_dist.push_back(calc_max_size(solu, as, 3));
    }

    // determine box shape, calculate relevant things and check the 
    // input for consistency
    enum boundary_enum { vacuum, rectangular, cubic, truncoct, triclinic } boundary = vacuum;

    double size_corr=1.0;

    switch(pbc->type()){
      case('t'):
	boundary = truncoct;
	size_corr = 2.0 * sqrt(3.0)/3.0;
	
	if(minwall.size()>1)
	  throw(gromos::Exception("probox", 
	     "For truncated octahedral boxes you can only specify one number for @minwall"));

	if(minwall.size()==0)
	  if(solu.box()[0] != solu.box()[1] ||
	     solu.box()[0] != solu.box()[2])

	    throw(gromos::Exception("probox", 
               "For truncated octahedral boxes, the specified boxsize should be the same in all dimensions"));

	break;

      case('r'):
	boundary = rectangular;
	
	if(minwall.size()==1) boundary = cubic;
	else if (minwall.size()==0 &&  
		 (solu.box()[0] == solu.box()[1] &&
		  solu.box()[0] == solu.box()[2]))
	  boundary = cubic;

	break;

      case('c'):
	if(!boxsize)
	  throw gromos::Exception("probox",
				  "boxsize has to be specified for triclinic");
	break;
	
      case('v'):
	throw(gromos::Exception("probox", 
           "Why are you running this program if @pbc is vacuum?"));
	break;
    }
    
    // size of the solvent box
    vector<double> solvent_box(3, 0.0);
    
    if (boxsize){
      if (minwall.size())
	throw(gromos::Exception("probox", 
				"Cannot specify both boxsize and minwall."));

      if(!solu.hasBox)
	throw gromos::Exception("probox",
				"If you specify boxsize, the box dimensions should be in "
				"the BOX block of the solute");

      if(pbc->type() == 't'){
	for(int i=0; i<3; ++i){
	  solvent_box[i] = solu.box()[0];
	}
      }
      else if (pbc->type() == 'r'){
	for(int i=0; i<3; ++i){
	  solvent_box[i] = solu.box()[i];
	}
      }
      else if(pbc->type() == 'c'){
	// calculate the dimension of a large enough box to encompass the triclinic box
	// first construct the four diagonals
	vector<Vec> diag(4);
	diag[0] = 0.5*(solu.box().K() + solu.box().L() + solu.box().M());
	diag[1] = 0.5*(solu.box().K() + solu.box().L() - solu.box().M());
	diag[2] = 0.5*(solu.box().K() - solu.box().L() + solu.box().M());
	diag[3] = 0.5*(solu.box().K() - solu.box().L() - solu.box().M());
	
	// find the maximum and minimum x,y,z, store these in boxsize
	for(int j=0; j<3; ++j){
	  double maxdim=diag[0][j], mindim=diag[0][j];
	  for(int i=1; i<4; ++i){
	    if(diag[i][j] > maxdim) maxdim=diag[i][j];
	    if(diag[i][j] < mindim) mindim=diag[i][j];
	    if(-diag[i][j] > maxdim) maxdim= -diag[i][j];
	    if(-diag[i][j] < mindim) mindim= -diag[i][j];
	  }
	  solvent_box[j] = maxdim - mindim;
	}
      } // triclinic
      else{
	throw gromos::Exception("probox",
				"unknown boundary condition");
      }

    } // boxsize
    else{
      if (minwall.size() == 0)
	throw gromos::Exception("probox",
				"either use a specified boxsize "
				"or give a minimum distance from solute to the walls");
      
      if (boundary == truncoct){
	for(int i=0; i<3; i++){
	  solu.box()[i] = size_corr*(max_dist[0] + 2 * minwall[0]);
	  solvent_box[i] = solu.box()[i];
	}
      }
      else{
	if (minwall.size() == 1){
	  for(int i=0; i<3; i++){
	    solu.box()[i] = max_dist[0]+2*minwall[0];
	    solvent_box[i] = solu.box()[i];
	  }
	}
	else{
	  // the solute has been rotated, 3 max_dist are known, 3 minwall distances are given
	  for(int i=0; i<3; i++){
	    solu.box()[i] = max_dist[2-i]+2*minwall[i];
	    solvent_box[i] = solu.box()[i];
	  }
	}
      }
    }
    
    // read in the solvent coordinates. 
    // to make absolutely sure that there is a box block, check this    
    ic.open(args["solvent"]);
    ic.select("SOLVENT");
    ic >> solv;
    ic.close();

    if(!solv.hasBox)
      throw gromos::Exception("probox", 
			      "Could not read BOX block from solvent "
			      "coordinates");
    
    int num_solv_atoms_per_box = solv.sol(0).numPos();
	 
    // get the number of original solvent atoms    
    int numSolventAtoms=solu.sol(0).numPos();
      
    // move the solute to the centre of geometry
    Vec shiftcog=fit::PositionUtils::shiftToCog(&solu);
    // shift the solvent atoms in there as well
    // Is done in shiftcog nowadays!!
    //for(int i=0; i< numSolventAtoms; i++){
    //  solu.sol(0).pos(i)+= shiftcog;
    //}
    
    // calculate the solvent cog
    Vec solv_cog(0.0,0.0,0.0);
    for(int i=0; i<solv.sol(0).numPos(); i++)
      solv_cog+=solv.sol(0).pos(i);
    solv_cog/=solv.sol(0).numPos();
    // move to solvent cog
    for(int i=0; i<solv.sol(0).numPos(); i++)
      solv.sol(0).pos(i) -= solv_cog;

    // set the solute box size, calculate how many solvent boxes we should
    // use; calculate where to move the initial box
    vector<int> needed_boxes;
    Vec move_solvent;
    for(int i=0; i<3; i++){
      needed_boxes.push_back(int(solvent_box[i] / solv.box()[i]) + 1);
      move_solvent[i]= -0.5 * (needed_boxes[i] - 1) * solv.box()[i];
    }

    // move the initial box so that after multiplication the complete box
    // is centered around the origin.
    for(int i=0; i< num_solv_atoms_per_box; i++)
      solv.sol(0).pos(i) += move_solvent;
    
    // do the multiplications
    for(int ix=0; ix < needed_boxes[0]; ix++){
      
      for(int iy=0; iy < needed_boxes[1]; iy++){
	
	for(int iz=0; iz < needed_boxes[2]; iz++){

	  if(ix != 0 || iy != 0 || iz != 0){
	    Vec shift(ix * solv.box()[0], iy * solv.box()[1], iz * solv.box()[2]);

	    for(int atom = 0; atom < num_solv_atoms_per_box; atom++){
	      solv.sol(0).addPos(solv.sol(0).pos(atom) + shift);
	    }
	  }
	}
      }
    }

    int num_atoms_per_solvent=solv.sol(0).topology().numAtoms();
    int num_solvent_molecules=solv.sol(0).numPos() / num_atoms_per_solvent;

    // now we have to keep only those waters that are inside the box and 
    // far enough away from the solute
    // we look at the centre of geometry of the solvent molecule

    Vec o(0.0,0.0,0.0);
    double min_init =
      solvent_box[0] * solvent_box[0] +
      solvent_box[1] * solvent_box[1] +
      solvent_box[2] * solvent_box[2];
    
    for(int i=0; i< num_solvent_molecules; i++){

      // calculate the centre of geometry of this solvent
      Vec sol_i(0.0,0.0,0.0);
      for(int j=0; j< num_atoms_per_solvent; j++)
	sol_i += solv.sol(0).pos(num_atoms_per_solvent * i + j);
      sol_i /= num_atoms_per_solvent;
      
      // are we inside the box
      Vec check = pbc->nearestImage(o, sol_i, solu.box());

      if(check[0]==sol_i[0] && 
	 check[1]==sol_i[1] && 
	 check[2]==sol_i[2]){
	// yes we are in the box
	// calculate the closest distance to any solute
	double min2 = min_init;
	for(int m=0; m < solu.numMolecules(); m++){
	  for(int a=0; a < solu.mol(m).numAtoms(); a++){
	    if(! solu.mol(m).topology().atom(a).isH() && 
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
	    solu.sol(0).addPos(solv.sol(0).pos(num_atoms_per_solvent * i + k));
	}
      }
    }

    ostringstream title;
    title << "Solvating " << args["solute"];
    if (args.count("gather") != -1)
      title << " (gathered) ";
    
    title << " in " << args["solvent"] 
	  << endl;

    title << "Box dimensions (";
    if (boundary == truncoct) title << "truncated octahedron";
    else if (boundary == cubic) title << "cubic";
    else if (boundary == triclinic) title << "triclinic";
    else if (boundary == rectangular) title << "rectangular";
    else throw gromos::Exception("probox", "wrong boundary!!!");
    
    title << ") were ";
    if(!boxsize){
      title << "calculated from maximum" << endl;
      title << "solute atom-atom distance";
      if(args.count("rotate") == -1)
	title << " (not rotated):" << endl;
      else{
	if (boundary == rectangular)
	  title << "s (x, y, z) (after rotation):" << endl;
	else
	  title << " (after rotation):" << endl;
      }
      
      for(unsigned int i=0; i<max_dist.size(); ++i){
	int index = 2*(max_dist.size() - i - 1);
	title << "\t" << max_dist[i] << " between atoms "  
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
    title << "Added " << (solu.sol(0).numPos()-numSolventAtoms)
      / num_atoms_per_solvent
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
  max[0] =calc_max_size(sys, as, 3);
  
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
  max[2] =calc_max_size(sys, as, 1);
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

