// shake_analysis calculates (non-bonded) interaction energies 
//                for atoms on which a shake failure has occured and 
//                tries to get a clue as to the reason

#include <cassert>

#include <gmath/Vec.h>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Property.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

using namespace std;
  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "atoms", "props", "cut", 
                    "eps", "kap", "top", "coord"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@atoms  <atoms for which shake fails>\n";
  usage += "\t@cut    <cut-off distance>\n";
  usage += "\t@eps    <epsilon for reaction field correction>\n";
  usage += "\t@kap    <kappa for reaction field correction>\n";
  usage += "\t@top    <number of non-bonded interactions per atom to print>\n";
  usage += "\t@coord  <coordinate file>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // declare the energy class
  Energy en(sys, gff, *pbc);

  //  set atoms
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  
  // get cutoff distance
  double cut=atof(args["cut"].c_str());
  en.setCutOff(cut);
  

  // get RF variables
  double eps=0.0, kap=0.0;
  if(args.count("eps")>0) eps=atof(args["eps"].c_str());
  if(args.count("kap")>0) kap=atof(args["kap"].c_str());
  en.setRF(eps,kap);

  // get top number
  int top=-1;
  if(args.count("top")>0) top=atoi(args["top"].c_str());
  
  // find the properties to calculate
  PropertyContainer props(sys);
  vector<int> num_prop(4,0);
  vector<int> prop_types;
  
  for(int m=0; m<sys.numMolecules(); m++){
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      if(atoms.findAtom(m,bi()[0])!=-1 ||
	 atoms.findAtom(m,bi()[1])!=-1){
	ostringstream os;
	os << "d%" << m+1 << ":" << bi()[0]+1 << "," << bi()[1]+1;
	props.addSpecifier(os.str());
	num_prop[0]++;
	prop_types.push_back(bi().type());
      }
    }
  }
  for(int m=0; m<sys.numMolecules(); m++){
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      if(atoms.findAtom(m,ai()[0])!=-1 ||
	 atoms.findAtom(m,ai()[1])!=-1 ||
	 atoms.findAtom(m,ai()[2])!=-1){
	ostringstream os;
	os << "a%" << m+1 << ":" << ai()[0]+1 << "," << ai()[1]+1 << ","
	   << ai()[2]+1;
	props.addSpecifier(os.str());
	num_prop[1]++;
	prop_types.push_back(ai().type());
      }
    }
  }
  for(int m=0; m<sys.numMolecules(); m++){
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      if(atoms.findAtom(m,ii()[0])!=-1 ||
	 atoms.findAtom(m,ii()[1])!=-1 ||
	 atoms.findAtom(m,ii()[2])!=-1 ||
	 atoms.findAtom(m,ii()[3])!=-1){
	ostringstream os;
	os << "t%" << m+1 << ":" << ii()[0]+1 << "," << ii()[1]+1 << ","
	   << ii()[2]+1 << "," << ii()[3]+1;
	props.addSpecifier(os.str());
	num_prop[2]++;
	prop_types.push_back(ii().type());
      }
    }
  }
  for(int m=0; m<sys.numMolecules(); m++){
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      if(atoms.findAtom(m,di()[0])!=-1 ||
	 atoms.findAtom(m,di()[1])!=-1 ||
	 atoms.findAtom(m,di()[2])!=-1 ||
	 atoms.findAtom(m,di()[3])!=-1){
	ostringstream os;
	os << "t%" << m+1 << ":" << di()[0]+1 << "," << di()[1]+1 << ","
	   << di()[2]+1 << "," << di()[3]+1;
	props.addSpecifier(os.str());
	num_prop[3]++;
	prop_types.push_back(di().type());
	
      }
    }
  }
  en.setProperties(props);
  
  // define input coordinate
  InG96 ic(args["coord"]);
  ic.select("ALL");
  ic >> sys;

  // we have to gather with any method to get covalent interactions 
  // and charge-groups connected
  pbc->gathergr();

  // calculate the covalent energies
  en.calcCov();

  // print out the covalent energies
  if(num_prop[0]+num_prop[1]+num_prop[2]+num_prop[3]){
    cout << "--------------------------------------------------------------"
	 << endl;
    
    cout << "Covalent interactions involving atoms ";
    for(int i=0; i<atoms.size(); i++){
      if(atoms.mol(i)<0) cout << "s";
      else cout << atoms.mol(i)+1;
      cout << ":" << atoms.atom(i)+1 << " ";
    }
    cout << endl;
  }
  
  int count=0;
  if(num_prop[0]){
    cout << endl << "BONDS :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(10) << "atom-"
	 << setw(12) << "atom-"
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(10) << "numbers"
	 << setw(12) << "names"
	 << setw(13) << "constant"
	 << endl;
  }
  for(int i=0; i<num_prop[0]; i++, count++){
    int type = prop_types[count];
    
    cout << setw(4) << props[count]->mols()[0]+1;
    cout << setw(5) << props[count]->atoms()[0]+1 << "-" 
	 << setw(4) << props[count]->atoms()[1]+1
	 << setw(7) << sys.mol(props[count]->mols()[0]).topology().
                          atom(props[count]->atoms()[0]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[1]).topology().
                          atom(props[count]->atoms()[1]).name();
    cout.precision(3);
    cout.setf(ios::scientific, ios::floatfield);
    
    cout << setw(13) << gff.bondType(type).fc();
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(10) << gff.bondType(type).b0();
    cout.precision(5);
    cout << setw(16) << props[count]->getValue();
    cout.setf(ios::scientific, ios::floatfield);
    cout << setw(16) << en.cov(count) << endl;
  }
  if(num_prop[1]){
    cout << endl << "ANGLES :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(15) << "atom-"
	 << setw(17) << "atom-"
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(15) << "numbers"
	 << setw(17) << "names"
	 << setw(13) << "constant"
	 << endl;
  }
  for(int i=0; i<num_prop[1]; i++, count++){
    int type = prop_types[count];
    cout << setw(4) << props[count]->mols()[0]+1;
    cout << setw(5) << props[count]->atoms()[0]+1 << "-" 
	 << setw(4) << props[count]->atoms()[1]+1 << "-" 
	 << setw(4) << props[count]->atoms()[2]+1
	 << setw(7) << sys.mol(props[count]->mols()[0]).topology().
                            atom(props[count]->atoms()[0]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[1]).topology().
                            atom(props[count]->atoms()[1]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[2]).topology().
                            atom(props[count]->atoms()[2]).name();
    cout.precision(3);
    cout.setf(ios::scientific, ios::floatfield);
    
    cout << setw(13) << gff.angleType(type).fc();
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(10) << gff.angleType(type).t0();
    cout.precision(5);
    cout << setw(16) << props[count]->getValue();
    cout.setf(ios::scientific, ios::floatfield);
    cout << setw(16) << en.cov(count) << endl;
  }
  if(num_prop[2]){
    cout << endl << "IMPROPER DIHEDRALS :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(20) << "atom-"
	 << setw(22) << "atom-"
	 << setw(13) << "force-"
	 << setw(10) << "b0"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(20) << "numbers"
	 << setw(22) << "names"
	 << setw(13) << "constant"
	 << endl;
  }
  for(int i=0; i<num_prop[2]; i++, count++){
    int type = prop_types[count];
    
    cout << setw(4) << props[count]->mols()[0]+1;
    cout << setw(5) << props[count]->atoms()[0]+1 << "-" 
	 << setw(4) << props[count]->atoms()[1]+1 << "-"
	 << setw(4) << props[count]->atoms()[2]+1 << "-" 
	 << setw(4) << props[count]->atoms()[3]+1
	 << setw(7) << sys.mol(props[count]->mols()[0]).topology().
                           atom(props[count]->atoms()[0]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[1]).topology().
                           atom(props[count]->atoms()[1]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[2]).topology().
                           atom(props[count]->atoms()[2]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[3]).topology().
                           atom(props[count]->atoms()[3]).name();
    cout.precision(3);
    cout.setf(ios::scientific, ios::floatfield);
    
    cout << setw(13) << gff.improperType(type).fc();
    cout.setf(ios::fixed, ios::floatfield);
    cout << setw(10) << gff.improperType(type).q0();
    cout.precision(5);
    cout << setw(16) << props[count]->getValue();
    cout.setf(ios::scientific, ios::floatfield);
    cout << setw(16) << en.cov(count) << endl;
  }
  if(num_prop[3]){
    cout << endl << "DIHEDRAL ANGLES :" << endl << endl;
    cout << setw(4) << "mol"
	 << setw(20) << "atom-"
	 << setw(22) << "atom-"
	 << setw(13) << "force-"
	 << setw(6)  << "pd"
	 << setw(4)  << "np"
	 << setw(16) << "b in x"
	 << setw(16) << "energy" << endl;
    cout << setw(4) << "# "
	 << setw(20) << "numbers"
	 << setw(22) << "names"
	 << setw(13) << "constant"
	 << endl;
  }
  for(int i=0; i<num_prop[3]; i++, count++){
    int type = prop_types[count];
    cout << setw(4) << props[count]->mols()[0]+1;
    cout << setw(5) << props[count]->atoms()[0]+1 << "-" 
	 << setw(4) << props[count]->atoms()[1]+1 << "-"
	 << setw(4) << props[count]->atoms()[2]+1 << "-" 
	 << setw(4) << props[count]->atoms()[3]+1
	 << setw(7) << sys.mol(props[count]->mols()[0]).topology().
                           atom(props[count]->atoms()[0]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[1]).topology().
                           atom(props[count]->atoms()[1]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[2]).topology().
                           atom(props[count]->atoms()[2]).name() << "-"
	 << setw(4) << sys.mol(props[count]->mols()[3]).topology().
                           atom(props[count]->atoms()[3]).name();
    cout.precision(3);
    cout.setf(ios::scientific, ios::floatfield);
    
    cout << setw(13) << gff.dihedralType(type).fc();
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(1);
    cout << setw(6) << gff.dihedralType(type).pd();
    cout << setw(4) << gff.dihedralType(type).np();
    
    cout.precision(5);
    cout << setw(16) << props[count]->getValue();
    cout.setf(ios::scientific, ios::floatfield);
    cout << setw(16) << en.cov(count) << endl;
  }
  
  // loop over the relevant atoms
  for(int i=0; i<atoms.size(); i++){
    
    //make a pairlist
    SimplePairlist pl(sys, *pbc, cut);
    pl.setAtom(atoms.mol(i), atoms.atom(i));
    pl.setType("CHARGEGROUP");
    pl.calc();
    pl.removeExclusions();

    int atom_i=pl.size();
    
    // add the atom itself to the pairlist
    pl.addAtom(atoms.mol(i), atoms.atom(i));

    // set these atoms for the energy class
    en.setAtoms(pl);
    
    // calculate all the interactions over the pairlist individually
    vector<double> vdw(atom_i), el(atom_i);
    for(int j=0; j<atom_i; j++){
      en.calcPair(j,atom_i, vdw[j], el[j]);
    }
    
    // filter out the largest absolute interactions
    set<int> inter;
    int max_inter, max_atom;
    if(top==-1 || top > atom_i) max_inter=atom_i;
    else max_inter=top;
    vector<int> order;
    double e_tot, max;
    
    for(int k=0; k<max_inter; k++){
      max=0.0;
      max_atom=0;
      
      for(int l=0; l< atom_i; l++){
	e_tot=fabs(vdw[l]+el[l]);
	if(e_tot>=max && inter.count(l)==0){
	  max=e_tot;
	  max_atom=l;
	}
      }
      inter.insert(max_atom);
      order.push_back(max_atom);
    }
    
    // now write out the max_inter interactions
    cout << endl;
    cout << "--------------------------------------------------------------"
	 << endl;
    if(max_inter != atom_i) cout << "Largest " << max_inter << " of ";
    cout << atom_i;
    cout << " non-bonded interactions with atom ";
    cout.precision(4);
    cout.setf(ios::right, ios::floatfield);
    if(pl.mol(atom_i) < 0) cout << "s";
    else cout << pl.mol(atom_i)+1;
    cout << ":" << pl.atom(atom_i)+1 << " (\"" << pl.name(atom_i)
	 << "\"; IAC: " << pl.iac(atom_i)+1 << "; charge: " 
	 << pl.charge(atom_i) << ")" << endl << endl;
    cout << setw(8) << "atom"
	 << setw(5) << "name"
	 << setw(5) << "IAC"
	 << setw(11) << "charge"
	 << setw(15) << "distance"
	 << setw(15) << "vdw"
	 << setw(15) << "coulomb"
	 << setw(15) << "total" << endl;
    
    for(int k=0; k<max_inter; k++){
      cout.precision(4);
      cout.setf(ios::right, ios::floatfield);
      ostringstream os_k;
      if(pl.mol(order[k]) < 0) os_k << "s";
      else os_k << pl.mol(order[k]) + 1;
      os_k << ":" << pl.atom(order[k]) + 1;
      cout << setw(8) << os_k.str()
	   << setw(5) << pl.name(order[k])
	   << setw(5) << pl.iac(order[k]) + 1
	   << setw(11) << pl.charge(order[k]);
      
      gmath::Vec v=pbc->nearestImage(*pl.coord(atom_i), 
				     *pl.coord(order[k]), sys.box());
      double d=(v - *pl.coord(atom_i)).abs();
      cout.precision(5);
      cout << setw(15) << d;
      cout.setf(ios::scientific, ios::floatfield);
      cout << setw(15) << vdw[order[k]]
	   << setw(15) << el[order[k]]
	   << setw(15) << vdw[order[k]]+el[order[k]];
      
      cout << endl;
    }
    cout << endl;
    
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







