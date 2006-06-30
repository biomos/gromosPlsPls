// progch generate coordinates for explicit hydrogens

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/Neighbours.h"
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
using namespace utils;
using namespace gmath;
using namespace bound;



int generate_coordinate(System *sys, GromosForceField *gff, int m, int a,
			 vector<int> h, vector<int> nh, int geom, double eps);

double find_bond(System *sys, GromosForceField *gff, int m, Bond b, double guess);

double find_angle(System *sys, GromosForceField *gff, int m, Angle a, double guess);
int find_dihedral(System *sys, int m, int i, int j, vector<int> h);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "coord", "tol", "pbc"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@coord <coordinates>\n";
  usage += "\t@tol   <tolerance (default 0.1 %)>\n";
  usage += "\t@pbc   <boundary conditions>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology make system and force field
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());
    
    // read in coordinates
    InG96 ic(args["coord"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();
    
    // read in the accuracy
    double eps=0.001;
    if(args.count("tol") > 0)
      eps=atof(args["tol"].c_str())/100.0;

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // gather the system!
    (*pbc.*gathmethod)();

    // a bit of an ugly hack for solvent molecules. The whole program has
    // been written for solutes, but sometimes we would have crystallographic
    // waters for which we want to do the same thing. 
    if(sys.sol(0).numPos()){
      MoleculeTopology mt;
      for(int a=0; a< sys.sol(0).topology().numAtoms(); a++){
	mt.addAtom(sys.sol(0).topology().atom(a));
	mt.setResNum(a,0);
      }
      mt.setResName(0,"SOLV");
      
      ConstraintIterator ci(sys.sol(0).topology());
      for(;ci; ++ci){
	gff.addBondType(BondType(1.0, ci().dist()));
	Bond b(ci()[0], ci()[1]);
	b.setType(gff.numBondTypes()-1);
	mt.addBond(b);
      }
      
      // add every solvent molecule as a solute
      int numSolvent=sys.sol(0).numPos()/sys.sol(0).topology().numAtoms();
      for(int i=0; i<numSolvent; i++){
	Molecule m(mt);
	m.initPos();
	for(int a=0; a< mt.numAtoms(); a++){
	  m.pos(a)=sys.sol(0).pos(i*sys.sol(0).topology().numAtoms() + a);
	}
	sys.addMolecule(m);
      }
      // and remove the original solvent
      sys.sol(0).setNumPos(0);
    }
    


    // initialize two counters
    int replaced=0, kept=0;
    
    // loop over all atoms
    for(int m=0; m< sys.numMolecules(); m++){

      // flag the atoms with mass 1.008 as hydrogens
      sys.mol(m).topology().setHmass(1.008);
      for(int a=0; a< sys.mol(m).numAtoms(); a++){

	// find the neighbours of this atom
	utils::Neighbours n(sys, m, a);

	// divide into hydrogens and non-hydrogens
	vector<int> h;
	vector<int> nh;

	for(unsigned int i=0; i< n.size(); i++){
	  if(sys.mol(m).topology().atom(n[i]).isH()) {
	    h.push_back(n[i]);
	  }
	  else
	    nh.push_back(n[i]);
	}
	
	// determine what kind of geometry they should be
	int geom=0;
	
	// only continue if we have hydrogens
	int numH=h.size();
	int numNH=nh.size();

	if(numH==1 && numNH==1) geom=1;
	if(numH==1 && numNH==2) geom=2;
	if(numH==2 && numNH==1) geom=3;
	if(numH==3 && numNH==1) geom=4;
	// crystallographic water
	if(numH==2 && numNH==0) geom=5;
	// nh4+
	if(numH==4 && numNH==0) geom=6;
        if(numH==1 && numNH==3) geom=7;
	if(numH==2 && numNH==2) geom=8;
	
	if(numH && !geom){
	  ostringstream os;
	  os << "Unexpected geometry for hydrogen bound to atom: " 
	     << m+1 << ":" << a+1 << endl;
	  throw(gromos::Exception("progch", os.str()));
	}
	// we have to have a geometry (this means that there are hydrogens)
	// and a should not be a hydrogen itself. (in the case of H2O we have
	// e.g. H-H bonds. These are only treated via the oxygen.
	if(geom && !sys.mol(m).topology().atom(a).isH()){
	  int r=generate_coordinate(&sys, &gff, m, a, h, nh, geom, eps);
	  replaced+=r;
	  kept+=(numH-r);
	}
      }
    }

    OutG96S oc(cout);
    oc.select("ALL");
    ostringstream os;
    os << "progch found " << replaced + kept << " hydrogen atoms in " 
       << args["coord"] << endl;
    os << kept << " were within " << eps*100 <<"% of minimum energy bond "
       << "length" << endl;
    os << replaced << " were assigned new coordinates based on geometry";
    
    
    oc.writeTitle(os.str());
    
    oc << sys;
    oc.close();
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
    
  }
  return 0;
}

int generate_coordinate(System *sys, GromosForceField *gff, int m, int a, 
			 vector<int> h, vector<int> nh, int geom, double eps)
{
  int count=0;
  switch(geom){
    case(1):
      {
	
	double bond=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	Vec v0=sys->mol(m).pos(a)      - sys->mol(m).pos(h[0]);
	if(fabs(v0.abs() - bond)/bond > eps){
	  
	  double angle=M_PI/180.0*find_angle(sys, gff, m, 
					   Angle(nh[0], a, h[0]), 109.5);
	  int fourth=find_dihedral(sys, m, nh[0], a, h);
	
	  Vec v1=sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0]);
	  Vec v2=sys->mol(m).pos(nh[0])  - sys->mol(m).pos(a);
	  Vec v4=v1.cross(v2);
	  Vec v5=v2.cross(v4).normalize();
	  Vec v6=bond*cos(angle) *v2.normalize() - bond*sin(angle)*v5;
	  
	  sys->mol(m).pos(h[0])=sys->mol(m).pos(a) +v6;
	  count++;
	}
	break;
      }
      
    case(2):
      {
	double bond=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	Vec v0=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	if(fabs(v0.abs() - bond)/bond > eps){
	  
	  Vec v1=sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a);
	  Vec v2=sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a);
	  Vec v3=(v1+v2).normalize();
	  sys->mol(m).pos(h[0])=sys->mol(m).pos(a) - bond*v3;
	  count++;
	}
	break;
      }
    case(3):
      {
	double bond1=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	double bond2=find_bond(sys, gff, m, Bond(a, h[1]), 0.1);
	Vec v01=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	Vec v02=sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
	if(fabs(v01.abs()-bond1)/bond1 > eps ||
	   fabs(v02.abs()-bond2)/bond2 > eps){
	  
	  double angle1=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(nh[0], a, h[0]), 120);
	  double angle2=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(nh[0], a, h[1]), 120);
	  double angle3=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(h[0], a, h[1]), 120);
	  int fourth = find_dihedral(sys, m, nh[0], a, h);
	  
	  Vec v1 = sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0]);
	  Vec v2 = (sys->mol(m).pos(nh[0])  - sys->mol(m).pos(a)).normalize();
	  Vec v4 = v1.cross(v2).normalize();
	  Vec v5 = v2.cross(v4).normalize();
	  Vec v6 = bond1 * cos(angle1) * v2 
	    - bond1 * sin(angle1) * v5;
	  double A=bond2*cos(angle2);
	  double B=( bond1*bond2*cos(angle3) - A*v2.dot(v6) ) / v5.dot(v6);
	  double C=sqrt(bond2*bond2 - A*A - B*B);
	  Vec v7 = A * v2  + B * v5 + C * v4;
	  
	  if(fabs(v01.abs()-bond1)/bond1 > eps) {
	    sys->mol(m).pos(h[0])=sys->mol(m).pos(a) + v6;
	    count++;
	  }
	  if(fabs(v02.abs()-bond2)/bond2 > eps){
	    sys->mol(m).pos(h[1])=sys->mol(m).pos(a) + v7;
	    count++;
	  }
	}
	break;
      }
    case(4):
      {
	// very similar to the non-planar type of above
	double bond1=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	double bond2=find_bond(sys, gff, m, Bond(a, h[1]), 0.1);
	double bond3=find_bond(sys, gff, m, Bond(a, h[2]), 0.1);
	Vec v01=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	Vec v02=sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
	Vec v03=sys->mol(m).pos(a) - sys->mol(m).pos(h[2]);

	if(fabs(v01.abs()-bond1)/bond1 > eps || 
	   fabs(v02.abs()-bond2)/bond2 > eps ||
	   fabs(v03.abs()-bond3)/bond3 > eps){
	  
	  double angle1=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(nh[0], a, h[0]), 109.5);
	  double angle2=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(nh[0], a, h[1]), 109.5);
	  double angle3=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(nh[0], a, h[2]), 109.5);
	  double angle4=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(h[0], a, h[1]), 109.5);
	  double angle5=M_PI/180.0*find_angle(sys, gff, m, 
					      Angle(h[0], a, h[2]), 109.5);
	  int fourth=find_dihedral(sys, m, nh[0], a, h);
	  
	  Vec v1=sys->mol(m).pos(fourth) - sys->mol(m).pos(nh[0]);
	  Vec v2=(sys->mol(m).pos(nh[0])  - sys->mol(m).pos(a)).normalize();
	  Vec v3=sys->mol(m).pos(a)      - sys->mol(m).pos(h[0]);
	  Vec v4=v1.cross(v2).normalize();
	  Vec v5=v2.cross(v4).normalize();
	  Vec v6=bond1*cos(angle1) * v2 - bond1 * sin(angle1)*v5;
	  
	  double A=bond2*cos(angle2);
	  double B=( bond1*bond2*cos(angle4) - A*v2.dot(v6) ) / v5.dot(v6);
	  double C=sqrt(bond2*bond2 - A*A - B*B);
	  Vec v7 = A * v2  + B * v5 + C * v4;
	  A=bond3*cos(angle3);
	  B=(bond1*bond3*cos(angle5) - A*v2.dot(v6) ) / v5.dot(v6);
	  C=sqrt(bond2*bond2 - A*A - B*B);
	  Vec v8 = A*v2 + B * v5 - C * v4;
	  if(fabs(v01.abs()-bond1)/bond1 > eps) {
	    sys->mol(m).pos(h[0])=sys->mol(m).pos(a) + v6;
	    count++;
	  }
	  if(fabs(v02.abs()-bond2)/bond2 > eps) {
	    sys->mol(m).pos(h[1])=sys->mol(m).pos(a) + v7;
	    count++;
	  }
	  if(fabs(v03.abs()-bond3)/bond3 > eps) {
	    sys->mol(m).pos(h[2])=sys->mol(m).pos(a) + v8;
	    count++;
	  }
	}
	break;
      }
    case(5):
      {
	// likely to be a water molecule. Here we have to come up with some
        // random orientation.
	double bond1=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	double bond2=find_bond(sys, gff, m, Bond(a, h[1]), 0.1);
	Vec v01=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	Vec v02=sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
	if(fabs(v01.abs()-bond1)/bond1 > eps || 
	   fabs(v02.abs()-bond2)/bond2 > eps){
	  // first generate a standard molecule. If it really is a water 
	  // molecule, there is probably no angle defined, but putting them 
	  // at 109.5 degrees is not so bad.
	  double angle=M_PI/180.0*find_angle(sys, gff, m, 
					     Angle(h[0], a, h[1]), 109.5);
	  Vec v1(0.0,0.0,bond1);
	  Vec v2(0.0,bond2*sin(angle), bond2*cos(angle));

	  //get three random numbers for the angles
	  //calculate sin and cos of these angles
	  
	  Vec angle_cos, angle_sin;
	  for(int i=0; i<3; i++){
	    double ang=2.0*M_PI/RAND_MAX*double(rand());
	    angle_cos[i]=cos(ang);
	    angle_sin[i]=sin(ang);
	  }
	  
	  // prepare a matrix to perform three random rotations  
	  // The product of three rotation matrices about three axes
	  /*
	   * (  1.0   0.0   0.0)   ( cosy   0.0   siny)   ( cosx   sinx   0.0)
	   * (  0.0  cosz  sinz) X (  0.0   1.0    0.0) X (-sinx   cosx   0.0)
	   * (  0.0 -sinz  cosz)   (-siny   0.0   cosy)   (  0.0    0.0   1.0)
	   */
	  gmath::Matrix rot(3,3);
	  rot(0,0)= angle_cos[0] * angle_cos[1];
	  rot(1,0)= -angle_sin[0] * angle_cos[2]
	    - angle_cos[0] * angle_sin[1] * angle_sin[2];
	  rot(2,0)= angle_sin[0] * angle_sin[2]
	    - angle_cos[0] * angle_sin[1] * angle_cos[2];
	  rot(0,1)= angle_sin[0] * angle_cos[1];
	  rot(1,1)= angle_cos[0] * angle_cos[2]
	    - angle_sin[0] * angle_sin[1] * angle_sin[2];
	  rot(2,1)= -angle_cos[0] * angle_sin[2]
	    - angle_sin[0] * angle_sin[1] * angle_cos[2];
	  rot(0,2)=angle_sin[1];
	  rot(1,2)=angle_cos[1] * angle_sin[2];
	  rot(2,2)=angle_cos[1] * angle_cos[2];
	
	  // rotate the hydrogens and put the coordinates
	  if(fabs(v01.abs()-bond1)/bond1 > eps) {
	    sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + rot*v1;
	    count++;
	  }
	  if(fabs(v02.abs()-bond2)/bond2 > eps) {
	    sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + rot*v2;
	    count++;
	  }
	}
	break;
      }
    case(6):
      {
	// nh4+, simliar to case 5
	double bond1=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	double bond2=find_bond(sys, gff, m, Bond(a, h[1]), 0.1);
	double bond3=find_bond(sys, gff, m, Bond(a, h[2]), 0.1);
	double bond4=find_bond(sys, gff, m, Bond(a, h[3]), 0.1);
	Vec v01=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	Vec v02=sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
	Vec v03=sys->mol(m).pos(a) - sys->mol(m).pos(h[2]);
	Vec v04=sys->mol(m).pos(a) - sys->mol(m).pos(h[3]);
	if(fabs(v01.abs()-bond1)/bond1 > eps || 
	   fabs(v02.abs()-bond2)/bond2 > eps ||
	   fabs(v03.abs()-bond3)/bond3 > eps ||
	   fabs(v04.abs()-bond4)/bond4 > eps){

	  // no angle search, theta is 109.5
	  // phi is 0, 120, 240
	  double angle=M_PI/180.0*109.5;
	  Vec v1(0.0, 0.0, bond1);
	  Vec v2(bond2*sin(angle)*cos(0.0),
		 bond2*sin(angle)*sin(0.0),
		 bond2*cos(angle));
	  Vec v3(bond2*sin(angle)*cos(M_PI/180.0*120.0),
		 bond2*sin(angle)*sin(M_PI/180.0*120.0),
		 bond2*cos(angle));
	  Vec v4(bond2*sin(angle)*cos(M_PI/180.0*240.0),
		 bond2*sin(angle)*sin(M_PI/180.0*240.0),
		 bond2*cos(angle));

	  //get three random numbers for the angles
	  //calculate sin and cos of these angles
	  
	  Vec angle_cos, angle_sin;
	  for(int i=0; i<3; i++){
	    double ang=2.0*M_PI/RAND_MAX*double(rand());
	    angle_cos[i]=cos(ang);
	    angle_sin[i]=sin(ang);
	  }

	  gmath::Matrix rot(3,3);
	  rot(0,0)= angle_cos[0] * angle_cos[1];
	  rot(1,0)= -angle_sin[0] * angle_cos[2]
	    - angle_cos[0] * angle_sin[1] * angle_sin[2];
	  rot(2,0)= angle_sin[0] * angle_sin[2]
	    - angle_cos[0] * angle_sin[1] * angle_cos[2];
	  rot(0,1)= angle_sin[0] * angle_cos[1];
	  rot(1,1)= angle_cos[0] * angle_cos[2]
	    - angle_sin[0] * angle_sin[1] * angle_sin[2];
	  rot(2,1)= -angle_cos[0] * angle_sin[2]
	    - angle_sin[0] * angle_sin[1] * angle_cos[2];
	  rot(0,2)=angle_sin[1];
	  rot(1,2)=angle_cos[1] * angle_sin[2];
	  rot(2,2)=angle_cos[1] * angle_cos[2];

	  // rotate the hydrogens and put the coordinates
	  if(fabs(v01.abs()-bond1)/bond1 > eps) {
	    sys->mol(m).pos(h[0]) = sys->mol(m).pos(a) + rot*v1;
	    count++;
	  }
	  if(fabs(v02.abs()-bond2)/bond2 > eps) {
	    sys->mol(m).pos(h[1]) = sys->mol(m).pos(a) + rot*v2;
	    count++;
	  }
	  if(fabs(v03.abs()-bond3)/bond3 > eps) {
	    sys->mol(m).pos(h[2]) = sys->mol(m).pos(a) + rot*v3;
	    count++;
	  }
	  if(fabs(v04.abs()-bond4)/bond4 > eps) {
	    sys->mol(m).pos(h[3]) = sys->mol(m).pos(a) + rot*v4;
	    count++;
	  }
	}
	break;
      }
    case(7):
      {
	// charged NH, connected to 3 NH atoms.
	double bond=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	Vec v0=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	if(fabs(v0.abs() - bond)/bond > eps){
	  
	  Vec v1=sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a);
	  Vec v2=sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a);
	  Vec v3=sys->mol(m).pos(nh[2]) - sys->mol(m).pos(a);
	  Vec v4=(v1+v2+v3).normalize();
	  sys->mol(m).pos(h[0])=sys->mol(m).pos(a) - bond*v4;
	  count++;
	}
	break;
      }
    case(8):
      {
	// charged NH2, connected to 2 NH atoms
	double bond1=find_bond(sys, gff, m, Bond(a, h[0]), 0.1);
	double bond2=find_bond(sys, gff, m, Bond(a, h[1]), 0.1);
	Vec v0=sys->mol(m).pos(a) - sys->mol(m).pos(h[0]);
	Vec v1=sys->mol(m).pos(a) - sys->mol(m).pos(h[1]);
	if(fabs(v0.abs() - bond1)/bond1 > eps ||
	   fabs(v1.abs() - bond2)/bond2 > eps){
	  
	  Vec v2=sys->mol(m).pos(nh[0]) - sys->mol(m).pos(a);
	  Vec v3=sys->mol(m).pos(nh[1]) - sys->mol(m).pos(a);
          Vec v4=-(v2+v3).normalize();
	  Vec v5=v2.cross(v3).normalize();
	  double angle=M_PI/180.0*find_angle(sys, gff, m, 
					     Angle(h[0], a, h[1]), 109.5);
	  
	  sys->mol(m).pos(h[0])=sys->mol(m).pos(a) +
	    bond1*sin(0.5*angle)*v5 + 
	    bond1*cos(0.5*angle)*v4;
	  sys->mol(m).pos(h[1])=sys->mol(m).pos(a) - 
	    bond2*sin(0.5*angle)*v5 + 
	    bond2*cos(0.5*angle)*v4;
	  count++;
	}
	break;
      }
 

  }
  return count;
}

double find_bond(System *sys, GromosForceField *gff, int m, Bond b, double guess)
{
  BondIterator bi(sys->mol(m).topology());
  double value=0.0;
  for(;bi; ++bi){
    if(bi()[0]==b[0] && bi()[1]==b[1]){
      value=gff->bondType(bi().type()).b0();
      break;
    }
  }
  if(value!=0.0)
    return value;
  else
    return guess;
}

double find_angle(System *sys, GromosForceField *gff, int m, Angle a, double guess)
{
  AngleIterator ai(sys->mol(m).topology());
  double value=0.0;
  
  for(;ai; ++ai){
    if(ai()[0]==a[0] && ai()[1] == a[1] && ai()[2] == a[2]){
      value=gff->angleType(ai().type()).t0();
      break;
    }
  }
  if(value!=0.0)
    return value;
  else
    return guess;
}

int find_dihedral(System *sys, int m, int i, int j, vector<int> h)
{
  
  if(j<i){
    int t=i;
    i=j;
    j=t;
  }
  
  DihedralIterator di(sys->mol(m).topology());
  int fourth=-1;
  
  for(; di; ++di){
    if(di()[1]==i && di()[2]==j){
      for(unsigned int k=0; k<h.size(); k++){
	if(di()[0]==h[k]) fourth=di()[3];
	if(di()[3]==h[k]) fourth=di()[0];
      }
      if(fourth!=-1) break;
    }
  }
  if(fourth==-1){
    
    // cannot find a dihedral, then just take the first atom bounded to i, whicih is not a hydrogen
    Neighbours n(*sys, m, i);
    for(unsigned int k=0; k<n.size(); k++){
      int hydrogen=0;
      for(unsigned int l=0; l<h.size(); l++)
	if(n[k]==h[l]) hydrogen=1;
      if(!hydrogen) fourth=n[k];
    }
    if(fourth==-1)
      throw(gromos::Exception("find_dihedral", "undefined position"));
  }
  
  return fourth;
  
}

