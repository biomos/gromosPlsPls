// tstrip.cc

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>

using namespace gcore;
using namespace gio;
using namespace args;

using namespace std;

class atomno
{
public:
  int gromosNum;
  std::string atomSpec;
  int index;
  atomno(int n, std::string s, int i)
  {
    gromosNum=n;
    atomSpec=s;
    index=i;
  }
};

int operator <(const atomno a, const atomno b)
{
  return a.gromosNum < b.gromosNum;
}

int main(int argc, char **argv){

  char *knowns[] = {"topo", "gromosnum", "atomspec"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@topo      <topology>\n";
  usage += "\t@gromosnum <gromos atom number>\n";
  usage += "\t@atomspec  <atomspecifier>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);


    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());


    std::set<atomno> Atomnos;
    utils::AtomSpecifier as(sys);
    
    
    // read and convert gromos numbers
    Arguments::const_iterator iter=args.lower_bound("gromosnum"),
      to=args.upper_bound("gromosnum");
    for(;iter!=to;++iter){
      int grom = atoi(iter->second.c_str())-1;
      int m=0;
      int a=grom;
      while(a>=sys.mol(m).numAtoms()){
	a-=sys.mol(m).numAtoms();
	m++;
	if(m>=sys.numMolecules()){
	  m=-1;
	  break;
	}
      }
      ostringstream os;
      if(m==-1) os << "s"; else os << m+1; os << ":" << a+1;
      as.addSpecifier(os.str());
      Atomnos.insert(atomno(grom+1,os.str(),as.size()-1));
    }
    
    // and do the inverse
    iter=args.lower_bound("atomspec");
    to=args.upper_bound("atomspec");
    for(;iter!=to; ++iter){
      utils::AtomSpecifier bs(sys, iter->second);
      for(int i=0; i < bs.size(); i++){
	int maxmol=bs.mol(i);
	if(maxmol<0) maxmol=sys.numMolecules();
	int grom=bs.atom(i);
	for(int j=0; j< maxmol; j++) grom+=sys.mol(j).numAtoms();
	ostringstream os;
	if(bs.mol(i)<0) os << "s"; else os << bs.mol(i)+1 ;
	os << ":" << bs.atom(i)+1;
	as.addSpecifier(os.str());
	Atomnos.insert(atomno(grom+1, os.str(), as.size()-1));
      }
    }
    if(Atomnos.size()){
      
      cout << setw(10) << "Atom"
	   << setw(10) << "GROMOS"
	   << setw(10) << "Residue"
	   << setw(10) << "Residue"
	   << setw(10) << "Atom"
	   << setw(12) << "Integer"
	   << setw(10) << "Charge" << endl;
      cout << setw(10) << "Specifier"
	   << setw(10) << "number"
	   << setw(10) << "number"
	   << setw(10) << "name"
	   << setw(10) << "name"
	   << setw(12) << "Atom Code"
	   << endl;
    }
    
    set<atomno>::const_iterator b=Atomnos.begin();
    set<atomno>::const_iterator e=Atomnos.end();
    for(;b!=e; ++b){
      int i=(*b).index;
      //determine the residue name
      int maxmol=as.mol(i);
      int resnum=0;
      std::string resname;
      if(maxmol < 0) {
	maxmol=0;
	resnum=as.atom(i)/sys.sol(0).topology().numAtoms();
	resname="SOLV";
      }
      else{
	
	resnum=sys.mol(as.mol(i)).topology().resNum(as.atom(i));
	resname=sys.mol(as.mol(i)).topology().resName(resnum);
      }
      
      for(int m=0; m<maxmol; m++) resnum+=sys.mol(m).topology().numRes();
      cout << setw(10) << (*b).atomSpec
	   << setw(10) << (*b).gromosNum
	   << setw(10) << resnum+1
	   << setw(10) << resname
	   << setw(10) << as.name(i)
	   << setw(12) << as.iac(i)+1
	   << setw(10) << as.charge(i)
	   << endl;
    }
    
      
	
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
