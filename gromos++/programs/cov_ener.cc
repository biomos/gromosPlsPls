//cov_ener

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

int findBond(Property &pp, System &sys);
int findAngle(Property &pp, System &sys);
int findTorsion(Property &pp, System &sys);
double calcBond(double val, Vec par);
double calcAngle(double val, Vec par);
double calcDihedral(double val, Vec par);
double calcImproper(double val, Vec par);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "prop", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";  
  usage += "\t@prop   <property specifier>\n";
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
  
    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    GromosForceField gff(it.forceField());
    // get atoms into AtomSpecifier
    PropertyContainer props(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("prop");
      Arguments::const_iterator to=args.upper_bound("prop");
      for(; iter!=to; iter++)
	{
	  string spec=iter->second.c_str();
	  props.addSpecifier(spec);
	}    
    }
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
    InG96 ic;

    // get force field parameters corresponding with these properties

    // set up an array to contain the forcefield data
    vector<gmath::Vec> par;

    for(unsigned int i=0;i<props.size();i++){
        Vec temp;
	switch(props[i]->atoms().size()){
	    case 2:{
		int t=findBond(*props[i], sys);
                temp[0]=gff.bondType(t).b0();
                temp[1]=gff.bondType(t).fc();
            }
            break;
            case 3:{
                int t=findAngle(*props[i], sys);
                temp[0]=gff.angleType(t).t0();
                temp[1]=gff.angleType(t).fc();
	    }
            break;
            case 4:{
                int t=findTorsion(*props[i], sys);
                // dirty fix to deal with both impropers and dihedrals
                // if t < 0 this means it is an improper, with types counting
		// -1, -2, ...
		// for the following program, the improper is recognised by 
		// the first parameter-value of -1000 (insensible for torsion)
                if(t>=0){
                  temp[0]=gff.dihedralType(t).pd();
                  temp[1]=gff.dihedralType(t).np();
                  temp[2]=gff.dihedralType(t).fc();
		}
                else{
		  t=-1*(t+1);
                  temp[0]=-1000;
                  temp[1]=gff.improperType(t).q0();
                  temp[2]=gff.improperType(t).fc();
		} 
	    }
	    break;
            default:
                throw gromos::Exception("cov_ener", 
		" I get a property with unallowed number of atoms");
	}
        par.push_back(temp);
    }

    // title

    cout << "#" << endl;
  
    cout << "#" << setw(9) << "time";
    cout << "\t\t" << props.toTitle();
    cout << endl;
    
    // set up array for averaging
    vector<double> cov;
    for(unsigned int i=0;i<props.size(); i++)
      cov.push_back(0.0);
    int num_frames=0;
    
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
   
      // calculate the props
      props.calc();

      // print the energies of the properties      
      cout << setw(10) << time << "\t\t";
      // cout << props;
      for(unsigned int i=0;i<props.size();i++){
          double en=0;
	  switch(props[i]->atoms().size()){
	      case 2:
		  en=calcBond(props[i]->getValue(), par[i]);
		  break;
	      case 3:
		  en=calcAngle(props[i]->getValue(), par[i]);
		  break;
	      case 4:
                  if(par[i][0]==-1000)
		    en=calcImproper(props[i]->getValue(), par[i]);
                  else
		    en=calcDihedral(props[i]->getValue(), par[i]);
		  break;
	  }
          cov[i]+=en;
	  
          cout << en  << "\t\t";
      }
      cout << endl;

      time += dt;
      num_frames++;
      
    }
    ic.close();
    
  }
   if(num_frames>1){
     cout << endl;
     cout << "#" << setw(9) << "ave.";
     for(unsigned int i=0; i< props.size(); i++)
       cout << "\t\t" << cov[i]/num_frames;
     cout << endl;
   }
   
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

int findBond(Property &pp, System &sys){
  int m, a, b, f=0;
  if(pp.mols()[0]==pp.mols()[1])
    m=pp.mols()[0];
  else
    throw gromos::Exception("findBond", 
          " covalent interactions are always within the same molecule");
  if(pp.atoms()[0]<pp.atoms()[1]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
  } 
  else {
    a=pp.atoms()[1];
    b=pp.atoms()[0];
  }
  BondIterator bi(sys.mol(m).topology());
  while(bi&&f==0)
    if(bi()[0]==a&&bi()[1]==b) f=1;
    else ++bi;
    if(bi) return bi().type();
    else
      throw gromos::Exception("findBond",
		    " bond not found: "+pp.toTitle());
}
int findAngle(Property &pp, System &sys){
  int m, a, b, c, f=0;
  if(pp.mols()[0]==pp.mols()[1]&&pp.mols()[0]==pp.mols()[2])
    m=pp.mols()[0];
  else
    throw gromos::Exception("findAngle", 
          " covalent interactions are always within the same molecule");
  if(pp.atoms()[0]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
  } 
  else {
    a=pp.atoms()[2];
    b=pp.atoms()[1];
    c=pp.atoms()[0];
  }
  AngleIterator ai(sys.mol(m).topology());
  while(ai&&f==0)
    if(ai()[0]==a&&ai()[1]==b&&ai()[2]==c) f=1;
    else ++ai;
    if(ai) return ai().type();
    else
      throw gromos::Exception("findAngle",
		    " angle not found: "+pp.toTitle());
}

int findTorsion(Property &pp, System &sys){
  int m, a, b, c, d, f=0;
  if(pp.mols()[0]==pp.mols()[1]&&pp.mols()[0]==pp.mols()[2]&&
     pp.mols()[0]==pp.mols()[3])
    m=pp.mols()[0];
  else
    throw gromos::Exception("findDihedral", 
          " covalent interactions are always within the same molecule");
  if(pp.atoms()[1]<pp.atoms()[2]){
    a=pp.atoms()[0];
    b=pp.atoms()[1];
    c=pp.atoms()[2];
    d=pp.atoms()[3];
  } 
  else {
    a=pp.atoms()[3];
    b=pp.atoms()[2];
    c=pp.atoms()[1];
    d=pp.atoms()[0];
  }
  DihedralIterator di(sys.mol(m).topology());
  while(di&&f==0)
    if(di()[0]==a&&di()[1]==b&&di()[2]==c&&di()[3]==d) f=1;
    else ++di;
  if(di) return di().type();
  else{
// maybe we have an improper
      ImproperIterator ii(sys.mol(m).topology());
      while(ii&&f==0)
	if(ii()[0]==a&&ii()[1]==b&&ii()[2]==c&&ii()[3]==d) f=1;
        else ++ii;
      if(ii) return -1*(ii().type()+1);
      else  
        throw gromos::Exception("findDihedral",
		    " (improperer) dihedral not found: "+pp.toTitle());
  }
}

//define pi as a global variable and put the manipulations of par[0] at the begining

double calcBond(double val, Vec par){
  double diff=val*val-par[0]*par[0];
  return 0.25*par[1]*diff*diff;
}
double calcAngle(double val, Vec par){
    val=val*3.1415927/180.0;
    double t0=par[0]*3.1415927/180.0;
    double diff=cos(val) - cos(t0);
    return 0.5*par[1]*diff*diff;
}
double calcDihedral(double val, Vec par){
    val=val*3.1415927/180.0;
    return par[2]*(1+par[0]*cos(par[1]*val));
}
double calcImproper(double val, Vec par){
    if(val>180) val=val-360;
    double diff=val-par[1];
    return 0.5*par[2]*diff*diff;
}
