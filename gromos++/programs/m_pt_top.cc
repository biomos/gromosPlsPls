#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "../src/args/Arguments.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/DihedralType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/ImproperType.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/PtTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InPtTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/OutTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

void printtopo(System &sys, vector<PtTopology> &pt, InTopology &it, vector<int> & comb, string title);
void printpttopo(System &sys, vector<PtTopology> &pt, vector<int> & comb, string title);
void printmpttopo(System &sys, vector<PtTopology> &pt, vector<vector<int> > & comb, string title);

int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "spec", "type"};
  int nknowns = 3;
  
  string usage = argv[0];
  usage += "\n\t@topo        <topology>\n";
  usage += "\t@type        <TOPO, PERTTOPO, or MPERT>\n";
  usage += "\t@spec        <specifications file>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read in topology
    InTopology it(args["topo"]);
    System sys(it.system());


    // first read in the specification data
    ifstream fin(args["spec"].c_str());
    int npt;
    fin >> npt;
    vector<string> ptname(npt);
    vector<string> f_atom(npt);
    for(int i=0; i< npt; i++) fin >> ptname[i] >> f_atom[i];
    
    int numcomb=0;
    fin >> numcomb;
    vector<vector <int> > comb;
    vector<int> cc(npt);
    for(int i=0; i< numcomb; ++i){
      for(int j=0; j<npt; j++){
	fin >> cc[j];
	cc[j]--;
      }
      comb.push_back(cc);
    }
      
    // Read the perturbations
    vector<PtTopology> pt;
    {
      for(int i=0; i< npt; i++){
	int start=0;
       
	utils::AtomSpecifier at(sys, f_atom[i]);
	int m_strt=at.mol(0);
	int a_strt=at.atom(0);
	for(int m=0; m<m_strt; m++)
	  start+=sys.mol(m).topology().numAtoms();
	start+=a_strt;
      
	InPtTopology iptt(ptname[i], start);
	pt.push_back(iptt.ptTopo());

	// check the names
	int m, a;
	for(int j=0; j< pt[i].numAtoms(); j++){
	  pt[i].findAtom(sys, m, a, j);
	  if(pt[i].atomName(j) != sys.mol(m).topology().atom(a).name())
	    throw gromos::Exception("m_pt_top",
				    "Atom names in perturbation topology do "
				    "not match topology\n"+pt[i].atomName(j) 
				    + " and " 
				    + sys.mol(m).topology().atom(a).name());
	}
      }
    }
    
    // create a title
    ostringstream title;
    if(args["type"]=="MPERT") title << "Multiple perturbation t";
    else if(args["type"]=="PERTTOPO") title << "Perturbation t";
    else title << "T";
    title << "opology based on\n";
    title << args["topo"] << " and\n";
    title << args["spec"];

    // what do we want for output?
    if(args["type"]=="TOPO"){
      if(comb.size()!=1)
	throw gromos::Exception("m_pt_top", 
				"writing of a topology is not possible for "
				"more than one perturbation. Give type "
				"MPERT\n");
      else
	printtopo(sys,pt,it,comb[0], title.str());
    }
    else if(args["type"]=="PERTTOPO"){
      if(comb.size()!=1)
	throw gromos::Exception("m_pt_top", 
				"writing of a perturbation topology is not "
				"possible for more than one perturbation. "
				"Give type MPERT\n");
      else
	printpttopo(sys,pt,comb[0], title.str());
    }
    else if(args["type"]=="MPERT"){
      printmpttopo(sys,pt,comb, title.str());
    }
    else
      throw gromos::Exception("pt_top", 
	     " type not recognized, use TOPO or PERTTOPO.\n");
    
    return 0;
  }
  catch(gromos::Exception e){
    cerr << e.what() << endl;
    return 1;
  }
}

void printtopo(System &sys, vector<PtTopology> &pt, InTopology &it, vector<int> &comb, string title)
{

  for(unsigned int i=0; i< comb.size(); ++i){
    pt[i].apply(sys, comb[i]);
  }
  OutTopology ot(cout);
  ot.setTitle(title);
  ot.write(sys,it.forceField());
}


void printpttopo(System &sys, vector<PtTopology> &pt, vector<int> &comb, string title)
{
  int numatoms=0;
  utils::AtomSpecifier pt_atoms(sys);
  int m=0,a=0;
  for(unsigned int i=0; i< pt.size(); ++i) {
    numatoms+=pt[i].numAtoms();
    for(int j=0; j< pt[i].numAtoms(); j++){
      pt[i].findAtom(sys,m,a,j);
      pt_atoms.addAtom(m,a);
    }
    
    pt[i].apply(sys, comb[i]);
  }
  pt_atoms.sort();
  
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";

  cout << "PERTATOM\n";
  cout << numatoms << endl;
  cout << "# JLA RESNR ATNAM     IACB      WMB      CGB ISCLJ  ISCC\n";

  // loop over the perturbed atoms
  for(int i=0; i<pt_atoms.size(); ++i){
    int mol=pt_atoms.mol(i);
    int atom=pt_atoms.atom(i);
    
    // the residue number in gromos96 continues over molecules
    int resoff=0, atomoff=0;
    for(int m=0;m<mol;m++){
      resoff+=sys.mol(m).topology().numRes();
      atomoff+=sys.mol(m).numAtoms();
    }
    
    // print out data
    cout.setf(ios::fixed);
    cout << setw(5) << atomoff + atom+1
         << setw(6) << sys.mol(mol).topology().resNum(atom)+resoff+1 << " "
         << setw(4) << pt_atoms.name(i) <<"\t";
    cout << setw(2) << pt_atoms.iac(i)+1
         << setw(9) << setprecision(4) 
         << pt_atoms.mass(i);
    cout << setw(9) << setprecision(3) << pt_atoms.charge(i);
    cout << "     1     1\n";
  }
  cout << "END\n";

  // and all the stuff that we do not take into account
  cout << "PERTATOMPAIR\n   0\nEND\n";
  cout << "PERTBONDH\n   0\nEND\n";
  cout << "PERTBOND\n   0\nEND\n";
  cout << "PERTBANGLEH\n   0\nEND\n";
  cout << "PERTBANGLE\n   0\nEND\n";
  cout << "PERTIMPDIHEDRALH\n   0\nEND\n";
  cout << "PERTIMPDIHEDRAL\n   0\nEND\n";
  cout << "PERTDIHEDRALH\n   0\nEND\n";
  cout << "PERTDIHEDRAL\n   0\nEND\n";
}

void printmpttopo(System &sys, vector<PtTopology> &pt, vector<vector<int> > & comb, 
		  string title)
{
 
  int numatoms=0;
  utils::AtomSpecifier pt_atoms(sys);
  int m=0,a=0;
  for(unsigned int i=0; i< pt.size(); ++i) {
    numatoms+=pt[i].numAtoms();
    for(int j=0; j< pt[i].numAtoms(); j++){
      pt[i].findAtom(sys,m,a,j);
      pt_atoms.addAtom(m,a);
    }
  }
  pt_atoms.sort();

  // make a new pert
  PtTopology pt_new;
  pt_new.setSize(pt_atoms.size(), comb.size());
  for(unsigned int i=0; i< comb.size(); i++){
    // apply perturbations
    for(unsigned int j=0; j< pt.size(); j++)
      pt[j].apply(sys, comb[i][j]) ;
    // fill pt_new
    for(int j=0; j< pt_atoms.size(); ++j){
      pt_new.setIac(j, i, pt_atoms.iac(j));
      pt_new.setCharge(j, i, pt_atoms.charge(j));
    }
  }
  
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";

  cout << "MPERTATOM\n# NUM_ATOMS NUM_PERTURBATIONS\n";
  cout << setw(6) << numatoms << setw(6) << comb.size() << endl;
  cout << "#  ATOM  PERT[1.." << comb.size() << "]\n"
       << "#         IAC CHARGE IAC CHARGE\n"
       << "          ";
  for(unsigned int i=0; i< comb.size(); ++i){
    for(unsigned int j=0; j< comb[i].size(); ++j){
      if(j!=0) cout << "_";
      cout << comb[i][j]+1;
    }
    cout << "    " ;
    
  }
  cout << "\n";
  
  

  // loop over the perturbed atoms
  for(int i=0; i<pt_atoms.size(); ++i){
    int atom=pt_atoms.atom(i);
    
    // print out data
    cout.setf(ios::fixed);
    cout << setw(5) << atom+1
         << setw(4) << pt_atoms.name(i);
    for(unsigned int j=0; j< comb.size(); j++)
      cout << setw(4) << pt_new.iac(i,j)+1
	   << setw(7) << setprecision(3) << pt_new.charge(i,j);
    cout << "\n";
  }
  cout << "END\n";
}
