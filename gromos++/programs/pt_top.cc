/**
 * @file pt_top.cc
 * Combine topologies and perturbation topologies
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pt_top
 * @section pt_top Combine topologies and perturbation topologies
 * @author @ref co
 * @date 7-6-07
 *
 * Combines topologies with perturbation topologies to produce new topologies 
 * or perturbation topologies. Reads a topology and a perturbation topology to
 * produce a new (perturbation) topology. Only non-bonded interactions (IAC and
 * charge) are implemented. The perturbation topology can contain a PERTATOM,
 * PERTATOM03 or MPERTATOM block (see volume IV). The atom numbers in the
 * perturbation topology do not need to match the numbers in the topology
 * exactly. If the topology and perturbation topology do not match in their 
 * atom numbering, a shift can be applied using the @firstatom option.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pttopo</td><td>&lt;perturbation topology with PERTATOM, PERTATOM03 or MPERTATOM block&gt; </td></tr>
 * <tr><td> \@type</td><td>&lt;output format: TOPO, PERTTOPO, or PERTTOPO03&gt; </td></tr>
 * <tr><td> \@npt</td><td>&lt;sequence number of the perturbation in a MPERTATOM block to apply&gt; </td></tr>
 * <tr><td> \@firstatom</td><td>&lt;@ref AtomSpecifier: first atom to which the perturbation will be applied&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pt_top
    @topo      ex.top
    @pttopo    ex.pttop
    @type      PERTATOM03
    @npt       1
    @firstatom 1:1
 @endverbatim
 *
 * <hr>
 */
#include <cassert>
#include <iostream>
#include <sstream>
#include <iomanip>
#include "../src/args/Arguments.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/OutTopology.h"
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
#include "../src/utils/AtomSpecifier.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace args;

class pert
{
  std::vector<int> d_num;
  std::vector<string> d_names;
  std::vector<string> d_pert;
  std::vector< std::vector <double> > d_charge;
  std::vector< std::vector <int> > d_iac;
  std::vector< std::vector <double> > d_mass;
  
public:
  pert(int a, int p);
  ~pert(){};
  void setIac(int a, int p, int iac);
  void setCharge(int a, int p, double q);
  void setMass(int a, int p, double m);
  void setName(int a, string name);
  void setPertName(int p, string name);
  void setNum(int a, int num);
  std::vector<double> charge(int p){return d_charge[p];}
  std::vector<int> iac(int p){return d_iac[p];}
  std::vector<double> mass(int p){return d_mass[p];}
  int iac(int a, int p){return d_iac[p][a];}
  double charge(int a, int p){return d_charge[p][a];}
  double mass(int a, int p){return d_mass[p][a];}
  string name(int a){return d_names[a];}
  string pertName(int p){return d_pert[p];}
  int num(int a){return d_num[a];}
  int numPt(){return d_iac.size();}
  int numAtoms(){return d_names.size();}
  
};

void findatom(System &sys, pert &pt, int &mol, int &atom, int counter);
void printtopo(System &sys, pert &pt, InTopology &it, int iipt, string title);
void printpttopo(System &sys, pert &pt, int iipt, string title);
void printpttopo03(System &sys, pert &pt, int iipt, string title);

int main(int argc, char *argv[]){

  char *knowns[] = {"topo", "pttopo", "firstatom", "npt", "type"};
  int nknowns = 5;
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pttopo    <perturbation topology with PERTATOM, PERTATOM03 or MPERTATOM block>\n";  
  usage += "\t@type      <output format: TOPO, PERTTOPO, or PERTTOPO03>\n";
  usage += "\t@npt       <sequence number of the perturbation in a MPERTATOM block to apply>\n";
  usage += "\t@firstatom <AtomSpecifier: first atom to which the perturbation will be applied>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read in topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //with which atom do we start?
    int start=0;
    {
      Arguments::const_iterator iter=args.lower_bound("firstatom");
      if(iter!=args.upper_bound("firstatom")){
        utils::AtomSpecifier at(sys, iter->second.c_str());
        int m_strt=at.mol(0);
        int a_strt=at.atom(0);
        for(int m=0; m<m_strt; m++)
          start+=sys.mol(m).topology().numAtoms();
        start+=a_strt;
      }
      else start=-1;
      
    }
    
    // read in multiple-perturbation-topology-file
    Ginstream ipt(args["pttopo"]);

    // define a few variables
    string nm, sdum;
    int spt=0;
    int a, p, k, l, iiac;
    double dq, fdum, dmass;

    // determine what type the perturbation topology is
    vector<string> buffer;
    ipt.getblock(buffer);
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("pt_top","PerturbationTopology file " + 
			      ipt.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    string ptstring;
    gio::concatenate(buffer.begin()+1, buffer.end()-1, ptstring);
    istringstream lineStream(ptstring);

    if(buffer[0]=="MPERTATOM"){ lineStream >> a >> p; spt=0; }
    else if(buffer[0]=="PERTATOM"){ lineStream >> a; p = 1; spt=1; }
    else if(buffer[0]=="PERTATOM03") { lineStream >> a; p = 1;spt=2;}
    else throw gromos::Exception("pt_top", 
       " Missing PERTATOM or MPERTATOM block in perturbation topology file"+buffer[0]); 
    

    // create perturbation class to contain all perturbation data
    pert pt(a,p);
    if(!spt){
      string bla;
      for(int i=0; i< pt.numPt(); i++){
	lineStream >> bla;
	pt.setPertName(i, bla);
      }
    }
  
    // read in the perturbation data for the atoms (only iac and q)
    for(int i=0; i<pt.numAtoms(); i++){
      lineStream >> k;
      if(i==0&&start>=0) start-=k;
      pt.setNum(i,k+start);
      if(spt) lineStream >> l;
      lineStream >> nm;
      pt.setName(i, nm);
      for(int j=0; j< pt.numPt(); j++){
        if(spt==2) lineStream >> fdum >> fdum >> fdum;
  	lineStream >> iiac;
        if(spt) {lineStream >> dmass; pt.setMass(i,j,dmass);}
        lineStream >> dq;
	pt.setIac(i,j,iiac-1);
	pt.setCharge(i,j,dq);
        if(spt) lineStream >> l >> l;
      }
    }

    // which perturbation do we use if we have read in an MPERTTOPO-block
    int iipt=0;
    {
      Arguments::const_iterator iter=args.lower_bound("npt");
      if(iter!=args.upper_bound("npt"))
        iipt=atoi(iter->second.c_str())-1;
      if(iipt>=pt.numPt()) throw gromos::Exception("pt_top", 
	  "Higher perturbation specified than available");
    }

    // create a title
    ostringstream title;
    if(args["type"]=="PERTTOPO") title << "Perturbation t";
    else title << "T";
    title << "opology based on\n";
    title << args["topo"] << " and\n";
    title << args["pttopo"];
    if(!spt || start>=0) title << " (";
    if(!spt) title << "perturbation " << iipt + 1;
    if(!spt && start>=0) title << "; ";
    if(start>=0) {
      int mol, atom;
      findatom(sys, pt, mol, atom, 0);
      title << "shifted to start at atomnumber " << pt.num(0)+1;
      title << " (" << mol+1 << ":" << atom+1 << ")";
    }
    if(!spt || start>=0) title << " )";
    
    // what do we want for output?
    if(args["type"]=="TOPO")
      printtopo(sys,pt,it,iipt, title.str());
    else if(args["type"]=="PERTTOPO")
      printpttopo(sys,pt,iipt, title.str());
    else if(args["type"]=="PERTTOPO03")
      printpttopo03(sys,pt,iipt, title.str());
    
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


pert::pert(int a, int p)
{
  for(int i=0;i<p;i++){
    
    vector<int> iac;
    vector<double> charge;
    vector<double> mass;
    
    for(int j=0;j<a;j++){
      iac.push_back(0);
      charge.push_back(0.0);
      mass.push_back(0.0);
      
    }
    d_iac.push_back(iac);
    d_charge.push_back(charge);
    d_mass.push_back(mass);
    d_pert.push_back("");
  }
  for(int j=0;j<a;j++){
    d_names.push_back(" ");
    d_num.push_back(0);
  }
}


void pert::setIac(int a, int p, int iac)
{
  d_iac[p][a]=iac;
}

void pert::setCharge(int a, int p, double q)
{
  d_charge[p][a]=q;
}
void pert::setMass(int a, int p, double m)
{
  d_mass[p][a]=m;
}

void pert::setName(int a, string name)
{
  d_names[a]=name;
}
void pert::setPertName(int p, string name)
{
  d_pert[p]=name;
}

void pert::setNum(int a, int num)
{
  d_num[a]=num;
}


void findatom(System &sys, pert &pt, int &mol, int &atom, int counter)
{
  if(counter>=pt.numAtoms()){
    mol=-1;
    atom=-1;
  }
  else{
    int at_cnt=0;
    for(mol=0;at_cnt<=pt.num(counter)&&mol<sys.numMolecules();mol++)
      at_cnt+=sys.mol(mol).topology().numAtoms();
    mol--;
    atom=pt.num(counter)-at_cnt+sys.mol(mol).topology().numAtoms();
  }
  
}

void printtopo(System &sys, pert &pt, InTopology &it, int iipt, string title)
{
    // create a system to output
    System syo;

    int counter=0;
    int mol, atom;

    //determine molecule and atom for this counter
    findatom(sys, pt, mol, atom, counter);
    
    // loop over the molecules 
    for(int m=0;m<sys.numMolecules();m++){

      if(m<mol||counter==pt.numAtoms()){
	// this molecule does not have perturbed atoms
	// just take the complete topology
	syo.addMolecule(sys.mol(m).topology());
      }
      else{
        MoleculeTopology mto;
        AtomTopology ato;
	
	// loop over the atoms in this molecule
	for(int aa=0;aa<sys.mol(m).topology().numAtoms();aa++){
	  ato=sys.mol(m).topology().atom(aa);
          if(aa==atom){
            if(pt.name(counter)!=sys.mol(m).topology().atom(aa).name()){
	      ostringstream os;
	      os << "Atom names in (perturbation) topologies do not match\n"
		 << "Topology: " << sys.mol(m).topology().atom(aa).name() 
		 << " (" << m+1 << ":" << aa+1 << ")"
		 << "\tPerturbation topology: " << pt.name(counter)
		 << " (" << counter << ")";
	      
	      
              throw gromos::Exception("pt_top", os.str());
	    }
	    
	    ato.setIac(pt.iac(counter, iipt));
	    ato.setCharge(pt.charge(counter, iipt));
	    ato.setMass(pt.mass(counter, iipt));
	    
  	    // loop to next atom in the perturbation list
            counter++;
	    findatom(sys, pt, mol, atom, counter);
	  }

	  // add the atom to the perturbed molecule
	  mto.addAtom(ato);
          int rsnm=sys.mol(m).topology().resNum(aa);
          mto.setResName(rsnm, sys.mol(m).topology().resName(rsnm));
	  mto.setResNum(aa, rsnm);

	}
	// do the bonds, angles etc in this molecule
        BondIterator bi(sys.mol(m).topology());
	for(;bi;++bi)
	  mto.addBond(bi());
	AngleIterator ai(sys.mol(m).topology());
	for(;ai;++ai)
	  mto.addAngle(ai());
	ImproperIterator ii(sys.mol(m).topology());
	for(;ii;++ii)
	  mto.addImproper(ii());
	DihedralIterator di(sys.mol(m).topology());
	for(;di;++di){
          mto.addDihedral(di());
	}

	// add the molecule to the output system	
	syo.addMolecule(mto);

      }
    }
    // we just take the solvent from the old system
    syo.addSolvent(sys.sol(0));

    // and write everything out
    OutTopology ot(cout);
    ot.setTitle(title);
    ot.write(syo,it.forceField());
}


void printpttopo(System &sys, pert &pt, int iipt, string title)
{
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";

  cout << "PERTATOM\n";
  cout << pt.numAtoms() << endl;
  cout << "# JLA RESNR ATNAM     IACB      WMB      CGB ISCLJ  ISCC\n";
  int mol=0, atom=0;
    
  // loop over the atoms
  for(int counter=0; counter < pt.numAtoms(); counter++){
    findatom(sys, pt, mol, atom, counter);

    // the residue number in gromos96 continues over molecules
    int resoff=0;
    for(int m=0;m<mol;m++)
      resoff+=sys.mol(m).topology().numRes();
      
    // check if atomnames are correct
    if(pt.name(counter)!=sys.mol(mol).topology().atom(atom).name())
      throw gromos::Exception("pt_top", 
              "Atom names in (perturbation) topologies do not match\n");

    // print out data
    cout.setf(ios::fixed);
    cout << setw(5) << pt.num(counter)+1
         << setw(6) << sys.mol(mol).topology().resNum(atom)+resoff+1 << " "
         << setw(4) << sys.mol(mol).topology().atom(atom).name() <<"\t";
    cout << setw(2) << pt.iac(counter, iipt)+1
         << setw(9) << setprecision(4) 
         << pt.mass(counter, iipt);
    cout << setw(9) << setprecision(3) << pt.charge(counter, iipt);
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

void printpttopo03(System &sys, pert &pt, int iipt, string title)
{
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";

  cout << "PERTATOM03\n";
  cout << pt.numAtoms() << endl;
  cout << "# JLA RESNR ATNAM IACA    MASSA  CHARGEA  IACB      WMB      CGB ISCLJ  ISCC\n";
  int mol=0, atom=0;
    
  // loop over the atoms
  for(int counter=0; counter < pt.numAtoms(); counter++){
    findatom(sys, pt, mol, atom, counter);

    // the residue number in gromos96 continues over molecules
    int resoff=0;
    for(int m=0;m<mol;m++)
      resoff+=sys.mol(m).topology().numRes();
      
    // check if atomnames are correct
    if(pt.name(counter)!=sys.mol(mol).topology().atom(atom).name())
      throw gromos::Exception("pt_top", 
              "Atom names in (perturbation) topologies do not match\n");

    // print out data
    cout.setf(ios::fixed);
    cout << setw(5) << pt.num(counter)+1
         << setw(6) << sys.mol(mol).topology().resNum(atom)+resoff+1 << " "
         << setw(4) << sys.mol(mol).topology().atom(atom).name() << " ";
    cout << setw(5) << sys.mol(mol).topology().atom(atom).iac()
	 << setw(9) << setprecision(4) 
	 << sys.mol(mol).topology().atom(atom).mass()
	 << setw(9) << setprecision(3) 
	 << sys.mol(mol).topology().atom(atom).charge() << " ";
    cout << setw(5) << pt.iac(counter, iipt)+1
         << setw(9) << setprecision(4) 
         << pt.mass(counter, iipt);
    cout << setw(9) << setprecision(3) << pt.charge(counter, iipt);
    if(sys.mol(mol).topology().atom(atom).iac()!=pt.iac(counter, iipt))
      cout << "     1";
    else
      cout << "     0";
    if(sys.mol(mol).topology().atom(atom).charge()!=pt.charge(counter, iipt))
      cout << "     1\n";
    else
      cout << "     0\n";
  }
  cout << "END\n";
  cout << "#\n";
  
}



   


