#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/bound/TruncOct.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Neighbours.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <set>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p);
int in_property(Property &p, int i);
int in_atoms(AtomSpecifier &as, int i);
void check_existing_property(System &sys, Property &p);
int findBond(System &sys, utils::Property &pp);
int findAngle(System &sys, utils::Property &pp);
int findDihedral(System &sys, utils::Property &pp);
void move_atoms(System &sys, AtomSpecifier &as, Vec v);
void rotate_atoms(System &sys, AtomSpecifier &as, gmath::Matrix rot, Vec v);
gmath::Matrix rotate_around(Vec v, double d);



int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "prop", "coord"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@prop <properties to change>\n";
  usage += "\t@coord <coordinate file>\n";
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    //  read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // parse boundary conditions
    Boundary *pbc= BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // input 
    InG96 ic(args["coord"]);

    // read in the coordinates and gather
    ic >> sys;
    (*pbc.*gathmethod)();
    ic.close();
    
    // output
    OutG96S oc(cout);
    
    //ic.close();
    // prepare the title
    ostringstream stitle;
    stitle << "progca has modified coordinates in " <<args["coord"]
	   << " such that";
    
    // what properties. 
    PropertyContainer props(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("prop"), 
	to=args.upper_bound("prop");
      if(iter==to)
	throw gromos::Exception("progca", 
				"no property specified");
      for(; iter!=to; ++iter)
	props.addSpecifier(iter->second.c_str());
    }

    // loop over the properties
    for(unsigned int i=0; i< props.size(); i++){
      
      // check whether this is actually a property we know from the topology
      // For the program, this is not necessary, but if you define arbitrary
      // properties, the decision on which atoms will move also gets 
      // arbitrary
      check_existing_property(sys, *props[i]);
      
      // store the target value
      double target=props[i]->getZValue();

      // calculate what the current value is      
      double value=props[i]->calc();
      
      // write this to the title
      stitle << endl << props[i]->toTitle()  << "\tgoes from " << value 
	     << "\tto " << target;
      // determine which atoms to move
      AtomSpecifier as(sys);
      atoms_to_change(sys, as, *props[i]);
      int m=props[i]->mols()[0];
      
      // now handle the three different cases
      switch(props[i]->toTitle()[0]){
	case 'd':{
	  Vec v = ( sys.mol(m).pos(props[i]->atoms()[1])
		    - sys.mol(m).pos(props[i]->atoms()[0])).normalize();
	  v *= (target-value);
	  move_atoms(sys, as, v);
	  break;
	}
	case 'a':{
	  Vec v1 = (sys.mol(m).pos(props[i]->atoms()[0])
		    - sys.mol(m).pos(props[i]->atoms()[1]));
	  Vec v2 = (sys.mol(m).pos(props[i]->atoms()[2])
		    - sys.mol(m).pos(props[i]->atoms()[1]));
	  Vec v3 = v1.cross(v2);
	  gmath::Matrix rot=rotate_around(v3, target-value);
	  rotate_atoms(sys, as, rot, sys.mol(m).pos(props[i]->atoms()[1]));
	  break;
	}
	case 't':{
	  Vec v = (sys.mol(m).pos(props[i]->atoms()[2])
		   - sys.mol(m).pos(props[i]->atoms()[1]));
	  gmath::Matrix rot=rotate_around(v, target-value);
	  rotate_atoms(sys, as, rot, sys.mol(m).pos(props[i]->atoms()[2]));
	  break;
	}
      }
    }
    oc.writeTitle(stitle.str());

    oc << sys;
    oc.close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

int in_property(Property &p, int i)
{
  // checks whether the atom i, is part of the property definition
  for(unsigned int j=0; j< p.atoms().size(); j++){
    if(i==p.atoms()[j]) return 1;
  }
  return 0;
}

int in_atoms(AtomSpecifier &as, int i)
{
  // checks whether atom i appears in the atom specifier
  // As our properties are always real bonds, angles or dihedrals, we do
  // not care about the molecule number here.
  // The AtomSpecifier should maybe get a function like this on its own.
  for(int j=0; j< as.size(); j++){
    if(i==as.atom(j)) return 1;
  }
  return 0;
}

void atoms_to_change(System &sys, AtomSpecifier &as, Property &p)
{
  // Determines which atoms of the whole molecule are influenced by changing
  // the value of the property

  // there are always two ends to a property. Regardless of the order of the
  // atoms we change the last end.
  int end2=p.atoms()[p.atoms().size()-1];

  // all atoms bounded to this atom but not part of the property 
  // will be modified
  set<int> a, new_a;
  a.insert(end2);
  // for bonds and angles we start with the last atom
  // for torsions we start with the third atom
  if(p.toTitle()[0]=='t') a.insert(p.atoms()[2]);
  else a.insert(end2);
  
  int m=p.mols()[0];
  
  while(a.begin()!=a.end()){
    set<int>::const_iterator b=a.begin(), e=a.end();
    new_a.clear();
    for(; b!=e; ++b){
      as.addAtom(m, *b);
      Neighbours n(sys, m, *b);
      for(unsigned int i=0; i<n.size(); i++)
	if(!in_property(p, n[i]) && !in_atoms(as, n[i]))
	   new_a.insert(n[i]);
    }
    a.clear();
    a=new_a;
  }
}

void move_atoms(System &sys, AtomSpecifier &as, Vec v)
{
  // Move the atoms in the atom-specifier by the vector v
  int m, a;
  
  for(int i=0; i<as.size(); i++){
    m=as.mol(i);
    a=as.atom(i);
    sys.mol(m).pos(a) += v;
  }
}
void rotate_atoms(System &sys, AtomSpecifier &as, gmath::Matrix rot, Vec v)
{
  // Rotate the atoms in the atom-specifyer according to the Matrix rot.
  // before rotation the atoms are moved so that v is the origin, after
  // rotation the atoms are moved back by v again.
  int m, a;
  Vec t;
  for(int i=0; i<as.size(); i++){
    m = as.mol(i);
    a = as.atom(i);
    t=sys.mol(m).pos(a) - v;
    t = rot*t;
    sys.mol(m).pos(a) = t + v;
  }
}
gmath::Matrix rotate_around(Vec v, double d)
{
  // Determine the rotation matrix that represents a rotation around the v
  // by a value of d (degrees)
  d*=M_PI/180.0;
  
  // first we create the matrix that brings v along z.
  // this is the product of two matrices
  // 1. around the x-axis by theta
  //    with sin(theta) = v[2]/r_yz; cos(theta) = align[2]/r_yz
  // 2. around the y-axis by phi
  // with sin(phi) = align[0]/r; cos(phi) = r_yz / r

  double r=v.abs();
  double r_yz = sqrt(v[1]*v[1] + v[2]*v[2]);
  gmath::Matrix rot1(Vec( r_yz / r         ,  0         , v[0]/r ),
		     Vec(-v[0]*v[1]/r/r_yz ,  v[2]/r_yz , v[1]/r ),
		     Vec(-v[0]*v[2]/r/r_yz , -v[1]/r_yz , v[2]/r ));

  // we also need its transpose to rotate back
  gmath::Matrix rot2(3,3);
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      rot2(i,j)=rot1(j,i);
  
  // and a Matrix that rotates around z by the angle d
  double cosd=cos(d);
  double sind=sin(d);
  gmath::Matrix m3(Vec(cosd, sind, 0.0),
		   Vec(-sind, cosd, 0.0),
		   Vec(0.0, 0.0, 1.0));
  
  // The matrix that we want is now rot2 * m3 * rot1
  return rot2*m3*rot1;
}

void check_existing_property(System &sys, Property &p)
{
  // Checks whether the property p is defined in the topology.
  char type=p.toTitle()[0];
  switch(type){
    case 'd':
      findBond(sys, p);
      break;
    case 'a':
      findAngle(sys, p);
      break;
    case 't':
      findDihedral(sys, p);
      break;
      
    default:
      throw(gromos::Exception("progca", 
			      "Unknown property type"+p.toTitle()));
  }
}

int findBond(System &sys, utils::Property &pp){
  // searches for the bond, defined by the property in the topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b,f=0;
  if(pp.mols()[0]==pp.mols()[1])
      m=pp.mols()[0];
  else
    throw gromos::Exception("progca",
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[0]<pp.atoms()[1]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
  }
  else{
    a=pp.atoms()[1];
    b=pp.atoms()[0];
  }
  BondIterator bi(sys.mol(m).topology());
  while(bi && f==0)
    if(bi()[0]==a && bi()[1]==b) return 1;
    else ++bi;
  throw gromos::Exception("progca", 
			  " Bond not found in topology: "+pp.toTitle());
}

 int findAngle(System &sys, utils::Property &pp){
  // searches for the angle, defined by the property in the topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b,c,f=0;
  if(pp.mols()[0]==pp.mols()[1]&&pp.mols()[0]==pp.mols()[2])
      m=pp.mols()[0];
  else
    throw gromos::Exception("progca",
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[0]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
  }
  else{
    a=pp.atoms()[2];
    b=pp.atoms()[1];
    c=pp.atoms()[0];
  }
  AngleIterator ai(sys.mol(m).topology());
  while(ai&&f==0)
      if(ai()[0]==a && ai()[1]==b && ai()[2]==c) return 1;
    else ++ai;
  throw gromos::Exception("progca", 
        " Angle not found in topology: "+pp.toTitle());
}
    
int findDihedral(System &sys, utils::Property &pp){
  // searches for the (improper) dihedral, defined by the property in the 
  // topology
  // returns 1 if found. Otherwise throws an exception.
  // Copied from the energy class, probably nicer to throw the exception 
  // from check_existing_property or from main
  int m, a, b, c, d, f=0;
  if(pp.mols()[0]==pp.mols()[1] &&
     pp.mols()[0]==pp.mols()[2] &&
     pp.mols()[0]==pp.mols()[3])
      m=pp.mols()[0];
  else
    throw gromos::Exception("progca", 
       " Covalent interactions are always within one molecule: "+pp.toTitle());
  if(pp.atoms()[1]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
    d=pp.atoms()[3];
  }
  else{
    a=pp.atoms()[3];
    b=pp.atoms()[2];
    c=pp.atoms()[1];
    d=pp.atoms()[0];
  }
  DihedralIterator di(sys.mol(m).topology());
  while(di&&f==0)
    if(di()[0]==a && di()[1]==b && di()[2]==c && di()[3]==d) return 1;
    else ++di;
  //Maybe we have an improper
  ImproperIterator ii(sys.mol(m).topology());
  while(ii&&f==0) 
    if(ii()[0]==a && ii()[1]==b && ii()[2]==c && ii()[3]==d) return 1;
    else ++ii;
  throw gromos::Exception("progca", 
			  " (improper) Dihedral not found in topology: "+pp.toTitle());
}

