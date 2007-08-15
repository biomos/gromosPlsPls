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
 * produce a new (perturbation) topology. The perturbation topology can contain
 * a PERTATOM, PERTATOM03 or MPERTATOM block (see volume IV). The atom numbers 
 * in the perturbation topology do not need to match the numbers in the topology
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
  std::vector< gcore::Bond > d_bonds;
  std::vector< gcore::Angle > d_angles;
  std::vector< gcore::Improper > d_impropers;
  std::vector< gcore::Dihedral > d_dihedrals;
  std::vector< gcore::AtomPair> d_atompair;
  std::vector< std::vector <int> > d_apdata;
  
public:
  pert(int a, int p);
  ~pert(){};
  pert &operator=(const pert &p);
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
  void addBond(int i, int j, int t);
  void addAngle(int i, int j, int k, int t);
  void addImproper(int i, int j, int k, int l, int t);
  void addDihedral(int i, int j, int k, int l, int t);
  void addAtomPair(int i, int j, int k, int l);
  int numBonds(){return d_bonds.size();}
  int numAngles(){return d_angles.size();}
  int numImpropers(){return d_impropers.size();}
  int numDihedrals(){return d_dihedrals.size();}  
  int numAtomPairs(){return d_atompair.size();}
  gcore::Bond bond(int i){return d_bonds[i];}
  gcore::Angle angle(int i){return d_angles[i];}
  gcore::Improper improper(int i){return d_impropers[i];}
  gcore::Dihedral dihedral(int i){return d_dihedrals[i];}
  gcore::AtomPair atompair(int i){return d_atompair[i];}
  std::vector<int> ap_data(int i){return d_apdata[i];}
  
    
    
};

pert readPERTATOM(System &sys, string line, int &start);
pert readMPERTATOM(System &sys, string line, int &start);
pert readPERTATOM03(System &sys, string line, int &start);
void readPERTBONDSTRETCH(System &sys, pert &pt, string line, int start);
void readPERTBONDANGLE(System &sys, pert &pt, string line, int start);
void readPERTIMPROPERDIH(System &sys, pert &pt, string line, int start);
void readPERTPROPERDIH(System &sys, pert &pt, string line, int start);

// gromos96 perturbation topology blocks --Clara Christ
void readPERTATOMPAIR(System &sys, pert &pt, string line, int start);
void readPERTBOND(GromosForceField &ff, System &sys, pert &pt, string line, int start);
void readPERTBANGLE(GromosForceField &ff, System &sys, pert &pt, string line, int start);
void readPERTIMPDIHEDRAL(GromosForceField &ff, System &sys, pert &pt, string line, int start);
void readPERTDIHEDRAL(GromosForceField &ff, System &sys, pert &pt, string line, int start);


void findatom(System &sys, pert &pt, int &mol, int &atom, int counter);
int findbond(System &sys, pert &pt, int m, gcore::Bond b);
int findangle(System &sys, pert &pt, int m, gcore::Angle b);
int findimproper(System &sys, pert &pt, int m, gcore::Improper b);
int finddihedral(System &sys, pert &pt, int m, gcore::Dihedral b);
int findexclusion(gcore::Exclusion e, int i);

void printtopo(System &sys, pert &pt, InTopology &it, int iipt, string title);
void printpttopo(System &sys, pert &pt, int iipt, string title, GromosForceField &ff);
void printpttopo03(System &sys, pert &pt, int iipt, string title, int format);

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
    // --CLARA
    GromosForceField ff(it.forceField());

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

    // create perturbation class to contain all perturbation data
    pert pt(0,0);
    
    // read the file 
    vector<string> buffer;
    while(ipt.getblock(buffer)){
      if(buffer[buffer.size()-1].find("END")!=0)
	throw gromos::Exception("pt_top","PerturbationTopology file " + 
				ipt.name() +
				" is corrupted. No END in "+buffer[0]+
				" block. Got\n"
				+ buffer[buffer.size()-1]);
      
      string ptstring;
      gio::concatenate(buffer.begin()+1, buffer.end()-1, ptstring);
 
      if(buffer[0]=="MPERTATOM"){ 
	pt = readMPERTATOM(sys, ptstring, start);}
      else if(buffer[0]=="PERTATOM"){ 
	pt = readPERTATOM(sys, ptstring, start);}
      else if(buffer[0]=="PERTATOM03" || buffer[0]=="PERTATOMPARAM") { 
	pt = readPERTATOM03(sys, ptstring, start);}
      else if(buffer[0]=="PERTBOND03" ||
	      buffer[0]=="PERTBONDSTRETCHH" || buffer[0]=="PERTBONDSTRETCH")
	readPERTBONDSTRETCH(sys, pt, ptstring, start);
      else if(buffer[0]=="PERTBANGLE03" ||
	      buffer[0]=="PERTBONDANGLEH" || buffer[0]=="PERTBONDANGLE")
	readPERTBONDANGLE(sys, pt, ptstring, start);
      else if(buffer[0]=="PERTIMPDIHEDRAL03" ||
	      buffer[0]=="PERTIMPROPERDIHH" || buffer[0]=="PERTIMPROPERDIH")
	readPERTIMPROPERDIH(sys, pt, ptstring, start);
      else if(buffer[0]=="PERTDIHEDRAL03" ||
	      buffer[0]=="PERTPROPERDIHH" || buffer[0]=="PERTPROPERDIH")
	readPERTPROPERDIH(sys, pt, ptstring, start);
      
      // gromos96 perturabation topology blocks --Clara      
      else if(buffer[0]=="PERTATOMPAIR")
	readPERTATOMPAIR(sys, pt, ptstring, start);
      else if(buffer[0]=="PERTBONDH" || buffer[0]=="PERTBOND")
	readPERTBOND(ff, sys, pt, ptstring, start);
      else if(buffer[0]=="PERTBANGLEH" || buffer[0]=="PERTBANGLE")
	readPERTBANGLE(ff, sys, pt, ptstring, start);
      else if(buffer[0]=="PERTIMPDIHEDRALH" || buffer[0]=="PERTIMPDIHEDRAL")
	readPERTIMPDIHEDRAL(ff, sys, pt, ptstring, start);
      else if(buffer[0]=="PERTDIHEDRALH" || buffer[0]=="PERTDIHEDRAL")
	readPERTDIHEDRAL(ff, sys, pt, ptstring, start);
      else
	std::cerr << "WARNING\n"
		  << "pt_top does not know how to handl block "
		  << buffer[0]
		  << "\ncontents will be ignored\n";
    }
    if(pt.numPt()==0)
      throw gromos::Exception("pt_top", 
			      " No atomic perturbation data read in!");
    
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
    int spt=0;
    if(pt.numPt()==1) spt=1;
    
    ostringstream title;
    if(args["type"]=="PERTTOPO") title << "Perturbation t";
    else title << "T";
    // note the t/T on the previous lines. Beauty is in the details.
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
      printpttopo(sys,pt,iipt, title.str(), ff);
    else if(args["type"]=="PERTTOPO03")
      printpttopo03(sys,pt,iipt, title.str(), 3);
    else if(args["type"]=="GROMOS05")
      printpttopo03(sys,pt,iipt, title.str(), 5);
    
    else
      throw gromos::Exception("pt_top", 
			      " type not recognized, use TOPO, PERTTOPO, PERTTOPO03 or GROMOS05.\n");
    
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

pert &pert::operator=(const pert::pert &p) 
{
  d_num = p.d_num;
  d_names = p.d_names;
  d_pert = p.d_pert;
  d_charge = p.d_charge;
  d_iac = p.d_iac;
  d_mass = p.d_mass;
  d_bonds = p.d_bonds;
  d_angles = p.d_angles;
  d_impropers = p.d_impropers;
  d_dihedrals = p.d_dihedrals;
  
  return *this;
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
void pert::addBond(int i, int j, int t)
{
  gcore::Bond b(i,j);
  b.setType(t);
  d_bonds.push_back(b);
}
void pert::addAngle(int i, int j, int k, int t)
{
  gcore::Angle a(i,j,k);
  a.setType(t);
  d_angles.push_back(a);
}
void pert::addImproper(int i, int j, int k, int l, int t)
{
  gcore::Improper ii(i,j,k,l);
  ii.setType(t);
  d_impropers.push_back(ii);
}
void pert::addDihedral(int i, int j, int k, int l, int t)
{
  gcore::Dihedral d(i,j,k,l);
  d.setType(t);
  d_dihedrals.push_back(d);
}
void pert::addAtomPair(int i, int j, int k, int l)
{
  gcore::AtomPair ap(i,j);
  std::vector<int> d;
  d.push_back(k);
  d.push_back(l);
  d_atompair.push_back(ap);
  d_apdata.push_back(d);
}

int findbond(System &sys, pert &pt, gcore::Bond b, int m){
  int t=-1;
  int offset=0;
  for(int mi=0; mi<m; mi++) offset+=sys.mol(mi).numAtoms();
  for(int i=0; i< pt.numBonds(); ++i){
    if(pt.bond(i)[0] == b[0] &&
       pt.bond(i)[1] == b[1])
      t=i;
  }
  return t;
}
int findangle(System &sys, pert &pt, gcore::Angle b, int m){
  int t=-1;
  int offset=0;
  for(int mi=0; mi<m; mi++) offset+=sys.mol(mi).numAtoms();
  for(int i=0; i< pt.numAngles(); ++i){
    if(pt.angle(i)[0] == b[0] &&
       pt.angle(i)[1] == b[1] &&
       pt.angle(i)[2] == b[2])
      t=i;
  }
  return t;
}
int findimproper(System &sys, pert &pt, gcore::Improper b, int m){
  int t=-1;
  int offset=0;
  for(int mi=0; mi<m; mi++) offset+=sys.mol(mi).numAtoms();
  for(int i=0; i< pt.numImpropers(); ++i){
    if(pt.improper(i)[0] == b[0] &&
       pt.improper(i)[1] == b[1] &&
       pt.improper(i)[2] == b[2] &&
       pt.improper(i)[3] == b[3])
      t=i;
    
  }
  return t;
}
int finddihedral(System &sys, pert &pt, gcore::Dihedral b, int m){
  int t=-1;
  int offset=0;
  for(int mi=0; mi<m; mi++) offset+=sys.mol(mi).numAtoms();
  for(int i=0; i< pt.numDihedrals(); ++i){
    if(pt.dihedral(i)[0] == b[0] &&
       pt.dihedral(i)[1] == b[1] &&
       pt.dihedral(i)[2] == b[2] &&
       pt.dihedral(i)[3] == b[3])
      t=i;
    
  }
  return t;
}
int findexclusion(gcore::Exclusion e, int i){
  int t=-1;
  for(int j=0; j< e.size(); j++)
    if(e.atom(j)==i) t=j;
  return t;
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

      MoleculeTopology mto;
      AtomTopology ato;
	
      // loop over the atoms in this molecule
      for(int aa=0;aa<sys.mol(m).topology().numAtoms();aa++){
	ato=sys.mol(m).topology().atom(aa);
	if(aa==atom){
	  if(pt.name(counter)!=sys.mol(m).topology().atom(aa).name()){
	    /*
	    ostringstream os;
	    os << "Atom names in (perturbation) topologies do not match\n"
	       << "Topology: " << sys.mol(m).topology().atom(aa).name() 
	       << " (" << m+1 << ":" << aa+1 << ")"
	       << "\tPerturbation topology: " << pt.name(counter)
	       << " (" << counter << ")";
	    
	    
	    throw gromos::Exception("pt_top", os.str());
	    */
	    cerr << "WARNING: Atom names in (perturbation) topologies do not match\n"
		 << "Topology: " << sys.mol(m).topology().atom(aa).name() 
		 << " (" << m+1 << ":" << aa+1 << ")"
		 << "\tPerturbation topology: " << pt.name(counter)
		 << " (" << counter << ")" << endl;
	  }
	  
	  ato.setIac(pt.iac(counter, iipt));
	  ato.setCharge(pt.charge(counter, iipt));
	  ato.setMass(pt.mass(counter, iipt));
	  
	  // loop to next atom in the perturbation list
	  counter++;
	  findatom(sys, pt, mol, atom, counter);
	}
	// is ato part of an perturbed AtomPair
	for(int i=0; i< pt.numAtomPairs(); ++i){
	  // we sorted the atompairs, so we only care if aa is the first of a
	  // pair
	  if(pt.atompair(i)[0]==aa){
	    //check if state A is correct
	    if(pt.ap_data(i)[0]==0){
	      if(findexclusion(ato.exclusion(), pt.atompair(i)[1])==-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 0 was specified "
		     << "for state A,\n"
		     << "         but they are not excluded in the topology.\n";
	      if(findexclusion(ato.exclusion14(), pt.atompair(i)[1])!=-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 0 was specified "
		     << "for state A,\n"
		     << "         but they are 1,4-neighbours in the topology.\n";
	    }
	    else if(pt.ap_data(i)[0]==1){
	      if(findexclusion(ato.exclusion(), pt.atompair(i)[1])!=-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 1 was specified "
		     << "for state A,\n"
		     << "         but they are excluded in the topology.\n";
	      if(findexclusion(ato.exclusion14(), pt.atompair(i)[1])!=-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 1 was specified "
		     << "for state A,\n"
		     << "         but they are 1,4-neighbours in the topology.\n";
	    }
	    else if(pt.ap_data(i)[0]==2){
	      if(findexclusion(ato.exclusion(), pt.atompair(i)[1])!=-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 2 was specified "
		     << "for state A,\n"
		     << "         but they are excluded in the topology.\n";
	      if(findexclusion(ato.exclusion14(), pt.atompair(i)[1])==-1)
		cerr << "WARNING: For PERTATOMPAIR " << pt.atompair(i)[0]+1
		     << " and " << pt.atompair(i)[1]+1 << " a 2 was specified "
		     << "for state A,\n"
		     << "         but they are not 1,4-neighbours in the topology.\n";
	    }
	    // do we need to do anything
	    if(pt.ap_data(i)[0]!=pt.ap_data(i)[1]){
	      Exclusion e=ato.exclusion();
	      Exclusion e14=ato.exclusion14();
	      if(pt.ap_data(i)[0]==0 && pt.ap_data(i)[1]==1){
		e.erase(pt.atompair(i)[1]);
	      }
	      else if(pt.ap_data(i)[0]==2 && pt.ap_data(i)[1]==1){
		e14.erase(pt.atompair(i)[1]);
	      }
	      else if(pt.ap_data(i)[0]==0 && pt.ap_data(i)[1]==2){
		e.erase(pt.atompair(i)[1]);
		e14.insert(pt.atompair(i)[1]);
	      }
	      else if(pt.ap_data(i)[0]==1 && pt.ap_data(i)[1]==2){
		e14.insert(pt.atompair(i)[1]);
	      }
	      else if(pt.ap_data(i)[0]==1 && pt.ap_data(i)[1]==0){
		e.insert(pt.atompair(i)[1]);
	      }
	      else if(pt.ap_data(i)[0]==2 && pt.ap_data(i)[2]==0){
		e14.erase(pt.atompair(i)[1]);
		e.insert(pt.atompair(i)[1]);
	      }
	      ato.setExclusion(e);
	      ato.setExclusion14(e14);
	    }
	  }
	}
	
	// add the atom to the perturbed molecule
	mto.addAtom(ato);
	int rsnm=sys.mol(m).topology().resNum(aa);
	mto.setResName(rsnm, sys.mol(m).topology().resName(rsnm));
	mto.setResNum(aa, rsnm);
	
      }
      // do the bonds, angles etc in this molecule
      // keep a list of the bonds that we found
      set<int> foundbonds;
      BondIterator bi(sys.mol(m).topology());
      for(;bi;++bi){
	gcore::Bond b=bi();
	// is this bond perturbed
	int tb=findbond(sys,pt,b,m);
	//cerr << "bonds, tb: " << tb << endl;
	if(tb>=0){
	  foundbonds.insert(tb);
	  b.setType(pt.bond(tb).type());
	}
	mto.addBond(b);
      }
      for(int i=0; i<pt.numBonds(); ++i)
	if(foundbonds.count(i)==0){
	  ostringstream os;
	  os << "Bond " << pt.bond(i)[0]+1 << " - " << pt.bond(i)[1]+1
	     << " not found in topology.";
	  throw gromos::Exception("pt_top", os.str());
	}
      
      set<int> foundangles;
      AngleIterator ai(sys.mol(m).topology());
      for(;ai;++ai){
	gcore::Angle b=ai();
	// is this angle perturbed
	int tb=findangle(sys, pt, b, m);
	//cerr << "angles, tb: " << tb << endl;
	if(tb>=0){
	  //cerr << "inserting angle of type " << pt.angle(tb).type() << endl;
	  foundangles.insert(tb);
	  b.setType(pt.angle(tb).type());
	}
	//mto.addAngle(ai());--CLARA
	mto.addAngle(b);
      }
      for(int i=0; i<pt.numAngles(); ++i)
	if(foundangles.count(i)==0){
	  ostringstream os;
	  os << "Angle " << pt.angle(i)[0]+1 << " - " << pt.angle(i)[1]+1
	     << " - " << pt.angle(i)[2]+1
	     << " not found in topology.";
	  throw gromos::Exception("pt_top", os.str());
	}
      
      set<int> foundimpropers;
      ImproperIterator ii(sys.mol(m).topology());
      for(;ii;++ii){
	gcore::Improper b=ii();
	// is this improper perturbed
	int tb=findimproper(sys,pt,b, m);
	//cerr << "impropers, tb: " << tb << endl;
	if(tb>=0){
	  foundimpropers.insert(tb);
	  b.setType(pt.improper(tb).type());
	}
	//mto.addImproper(ii()); --CLARA
	mto.addImproper(b);
      }
      for(int i=0; i<pt.numImpropers(); ++i)
	if(foundimpropers.count(i)==0){
	  ostringstream os;
	  os << "Improper " << pt.improper(i)[0]+1 << " - " << pt.improper(i)[1]+1
	     << " - " << pt.improper(i)[2]+1 << " - " << pt.improper(i)[3]+1
	     << " not found in topology.";
	  throw gromos::Exception("pt_top", os.str());
	}
    
      set<int> founddihedrals;
      DihedralIterator di(sys.mol(m).topology());
      for(;di;++di){
	gcore::Dihedral b=di();
	// is this dihedral perturbed
	int tb=finddihedral(sys,pt,b, m);
	//cerr << "dihedrals, tb: " << tb << endl;
	if(tb>=0){
	  founddihedrals.insert(tb);
	  b.setType(pt.dihedral(tb).type());
	}
	//	mto.addDihedral(di()); -- CLARA
	mto.addDihedral(b);
      }
      for(int i=0; i<pt.numDihedrals(); ++i)
	if(founddihedrals.count(i)==0){
	  ostringstream os;
	  os << "Dihedral " << pt.dihedral(i)[0]+1 << " - " << pt.dihedral(i)[1]+1
	     << " - " << pt.dihedral(i)[2]+1 << " - " << pt.dihedral(i)[3]+1
	     << " not found in topology.";
	  throw gromos::Exception("pt_top", os.str());
	}

      // add the molecule to the output system	
      syo.addMolecule(mto);
      
    }
    // we just take the solvent from the old system
    syo.addSolvent(sys.sol(0));

    // and write everything out
    OutTopology ot(cout);
    ot.setTitle(title);
    ot.write(syo,it.forceField());
}


void printpttopo(System &sys, pert &pt, int iipt, string title, GromosForceField &ff)
{
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";

  cout << "PERTATOM\n#\n# NJLA: number of perturbed atoms\n";
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

  // atom pairs
  cout << "PERTATOMPAIR\n"
       << "#\n# NEB: number of perturbed atom pairs\n"
       << setw(5) << pt.numAtomPairs() << endl;
  if(pt.numAtomPairs())
    cout << "# IEB, JEB: atom sequence numbers\n"
	 << "# IETA, IETB: type of non-bonded interactions (0,1,2)\n"
	 << "#             in state A and state B\n"
	 << "#             0: excluded\n"
	 << "#             1: normal\n"
	 << "#             2: 1,4-interaction\n";
  cout << "# IEB  JEB IETA IETB\n";
  for(int i=0; i<pt.numAtomPairs(); ++i){
    cout << setw(5) << pt.atompair(i)[0]+1
	 << setw(5) << pt.atompair(i)[1]+1
	 << setw(5) << pt.ap_data(i)[0]
	 << setw(5) << pt.ap_data(i)[1] << endl;
  }
  cout << "END\n";
  

  // Now do the bonded stuff
  if(pt.numBonds())
    cerr << "Warning: pt_top does not distinguish bonds with or without H-atoms.\n"
	 << "Warning: pt_top cannot determine the bond-sequence number correctly.\n";
  
	 
  cout << "PERTBOND\n#\n# NBONG: number of perturbed bonds\n";
  cout << setw(5) << pt.numBonds();
  cout << "\n# IBG  JBG NCBG            CBA   BA0            CBB   BB0\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundbonds;
  for(int m=0; m< sys.numMolecules();++m){
    int cb=0;
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi, ++cb){
      // is this bond perturbed
      int tb=findbond(sys,pt,bi(),m);
      if(tb>=0){
	cout << setw(5) << pt.bond(tb)[0]+1
	     << setw(5) << pt.bond(tb)[1]+1
	     << setw(5) << cb+1 << ' '
	     << setw(14) << setprecision(1) 
	     << ff.bondType(bi().type()).fc() << ' ' 
	     << setw(5) << setprecision(3)
	     << ff.bondType(bi().type()).b0() << ' '
	     << setw(14) << setprecision(1) 
	     << ff.bondType(pt.bond(tb).type()).fc() << ' '
	     << setw(5) << setprecision(3) 
	     << ff.bondType(pt.bond(tb).type()).b0()
	     << endl;
	
	foundbonds.insert(tb);
      }
    }
  }
  cout << "END\n";
  
  for(int i=0; i<pt.numBonds(); ++i)
    if(foundbonds.count(i)==0){
      ostringstream os;
      os << "Bond " << pt.bond(i)[0]+1 << " - " << pt.bond(i)[1]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }
  if(pt.numAngles())
    cerr << "Warning: pt_top does not distinguish angles with or without H-atoms.\n";
  cout << "PERTBANGLE\n"
       << "# number of perturbed angles\n# NTHEG\n"
       << setw(5) << pt.numAngles() << endl;
  cout << "# ITG  JTG  KTG    CTA   TA0    CTB   TB0\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundangles;
  for(int m=0; m< sys.numMolecules(); ++m){
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      // is this angle perturbed
      int tb=findangle(sys,pt,ai(),m);
      if(tb>=0){
	cout << setw(5) << pt.angle(tb)[0]+1
	     << setw(5) << pt.angle(tb)[1]+1
	     << setw(5) << pt.angle(tb)[2]+1 << ' '
	     << setw(6) << setprecision(1) 
	     << ff.angleType(ai().type()).fc() << ' '
	     << setw(5) << setprecision(1)
	     << ff.angleType(ai().type()).t0() << ' '
	     << setw(6) << setprecision(1) 
	     << ff.angleType(pt.angle(tb).type()).fc() << ' '
	     << setw(5) << setprecision(1)
	     << ff.angleType(pt.angle(tb).type()).t0()
	     << endl;
	
	foundangles.insert(tb);
      }
    }
  }
  cout << "END\n";
  
  for(int i=0; i< pt.numAngles();++i)
    if(foundangles.count(i)==0){
      ostringstream os;
      os << "Angle " << pt.angle(i)[0]+1 << " - " << pt.angle(i)[1]+1
	 << " - " << pt.angle(i)[2]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }

  if(pt.numImpropers())
    cerr << "Warning: pt_top does not distinguish impropers with or without H-atoms.\n";
  cout << "PERTIMPDIHEDRAL\n"
       << "# number of perturbed improper dihedrals\n# NQHIG\n"
       << setw(5) << pt.numImpropers() << endl;
  cout << "# IQG  JQG  KQG  LQG    CQA       QA0    CQV       QB0\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundimpropers;
  for(int m=0; m<sys.numMolecules(); m++){
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      // is this angle perturbed
      int tb=findimproper(sys,pt,ii(), m);
      if(tb>=0){
	cout << setw(5) << pt.improper(tb)[0]+1
	     << setw(5) << pt.improper(tb)[1]+1
	     << setw(5) << pt.improper(tb)[2]+1
	     << setw(5) << pt.improper(tb)[3]+1 << ' '
	     << setw(6) << setprecision(3) 
	     << ff.improperType(ii().type()).fc() << ' '
	     << setw(9) << setprecision(5)
	     << ff.improperType(ii().type()).q0() << ' '
	     << setw(6) << setprecision(3)
	     << ff.improperType(pt.improper(tb).type()).fc() << ' '
	     << setw(9) << setprecision(5)
	     << ff.improperType(pt.improper(tb).type()).q0()
	     << endl;
	
	foundimpropers.insert(tb);
      }
    }
  }
  cout << "END\n";
  for(int i=0; i<pt.numImpropers(); ++i)
    if(foundimpropers.count(i)==0){
      ostringstream os;
      os << "Improper " << pt.improper(i)[0]+1 << " - " << pt.improper(i)[1]+1
	 << " - " << pt.improper(i)[2]+1 << " - " << pt.improper(i)[3]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }
  if(pt.numDihedrals())
    cerr << "Warning: pt_top does not distinguish dihedrals with or without H-atoms.\n";
  
  cout << "PERTDIHEDRAL\n"
       << "# number of perturbed dihedrals\n#NPHIG\n"
       << setw(5) << pt.numDihedrals() << endl;
  cout << "# IBG  JBG  KBG  LBG    CPA   PDA   NPA    CPB   PDB   NPB\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> founddihedrals;
  for(int m=0; m<sys.numMolecules(); m++){
    
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      // is this angle perturbed
      int tb=finddihedral(sys,pt,di(), m);
      if(tb>=0){
	int offset=1;
	findatom(sys,pt,mol,atom,pt.dihedral(tb)[0]);
	for(int m=0;m<mol;++m) offset+=sys.mol(m).numAtoms();
	cout << setw(5) << pt.dihedral(tb)[0]+1
	     << setw(5) << pt.dihedral(tb)[1]+1
	     << setw(5) << pt.dihedral(tb)[2]+1
	     << setw(5) << pt.dihedral(tb)[3]+1 << ' '
	     << setw(6) << setprecision(2)
	     << ff.dihedralType(di().type()).fc() << ' '
	     << setw(5) << setprecision(1)
	     << ff.dihedralType(di().type()).pd() << ' '
	     << setw(5) 
	     << ff.dihedralType(di().type()).np() << ' '
	     << setw(6) << setprecision(2)
	     << ff.dihedralType(pt.dihedral(tb).type()).fc() << ' '
	     << setw(5) << setprecision(1)
	     << ff.dihedralType(pt.dihedral(tb).type()).pd() << ' '
	     << setw(5)
	     << ff.dihedralType(pt.dihedral(tb).type()).np()
	     << endl;
	
	founddihedrals.insert(tb);
      }
    }
  }
  cout << "END\n";

  for(int i=0; i<pt.numDihedrals(); ++i)
    if(founddihedrals.count(i)==0){
      ostringstream os;
      os << "Dihedral " << pt.dihedral(i)[0]+1 << " - " << pt.dihedral(i)[1]+1
	 << " - " << pt.dihedral(i)[2]+1 << " - " << pt.dihedral(i)[3]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }



}

void printpttopo03(System &sys, pert &pt, int iipt, string title, int format)
{
  // not so nicely written
  cout << "TITLE\n";
  cout << title;
  cout << "\nEND\n";
  
  if(format==3){
    cout << "PERTATOM03\n";
  }
  else{
    cout << "PERTATOMPARAM\n";
  }
  cout << "# NJLA\n";
  cout << setw(6) << pt.numAtoms() << endl;
  cout << "# JLA RESNR ATNM  IACA      WMA      CGA  IACB      WMB      CGB ALPHLJ ALPHCRF\n";

  
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
    cout << setw(5) << sys.mol(mol).topology().atom(atom).iac()+1
	 << setw(9) << setprecision(4) 
	 << sys.mol(mol).topology().atom(atom).mass()
	 << setw(9) << setprecision(3) 
	 << sys.mol(mol).topology().atom(atom).charge() << " ";
    cout << setw(5) << pt.iac(counter, iipt)+1
	 << setw(9) << setprecision(4) 
	 << pt.mass(counter, iipt);
    cout << setw(9) << setprecision(3) << pt.charge(counter, iipt);
    if(sys.mol(mol).topology().atom(atom).iac()!=pt.iac(counter, iipt))
      cout << "    1.0";
    else
      cout << "    0.0";
    if(sys.mol(mol).topology().atom(atom).charge()!=pt.charge(counter, iipt))
      cout << "     1.0\n";
    else
      cout << "     0.0\n";
  }
  
  cout << "END\n";
  cout << "#\n";

  // atom pairs
  cout << "PERTATOMPAIR\n"
       << "#\n# NEB: number of perturbed atom pairs\n"
       << setw(5) << pt.numAtomPairs() << endl;
  if(pt.numAtomPairs())
    cout << "# IEB, JEB: atom sequence numbers\n"
	 << "# IETA, IETB: type of non-bonded interactions (0,1,2)\n"
	 << "#             in state A and state B\n"
	 << "#             0: excluded\n"
	 << "#             1: normal\n"
	 << "#             2: 1,4-interaction\n";
  cout << "# IEB  JEB IETA IETB\n";
  for(int i=0; i<pt.numAtomPairs(); ++i){
    cout << setw(5) << pt.atompair(i)[0]+1
	 << setw(5) << pt.atompair(i)[1]+1
	 << setw(5) << pt.ap_data(i)[0]
	 << setw(5) << pt.ap_data(i)[1] << endl;
  }
  cout << "END\n";

  if(format==3)
    cout << "PERTBOND03\n";
  else
    cout << "PERTBONDSTRETCH\n";
  
  cout << "# number of perturbed bonds\n# NBONG\n"
       << setw(5) << pt.numBonds() << endl;
  cout << "# IBG  JBG ICBA ICBB\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundbonds;
  for(int m=0; m< sys.numMolecules();++m){
      
    BondIterator bi(sys.mol(m).topology());
    for(;bi;++bi){
      // is this bond perturbed
      int tb=findbond(sys,pt,bi(),m);
      if(tb>=0){
	cout << setw(5) << pt.bond(tb)[0]+1
	     << setw(5) << pt.bond(tb)[1]+1
	     << setw(5) << bi().type()+1
	     << setw(5) << pt.bond(tb).type()+1
	     << endl;
	
	foundbonds.insert(tb);
      }
    }
  }
  cout << "END\n";

  for(int i=0; i<pt.numBonds(); ++i)
    if(foundbonds.count(i)==0){
      ostringstream os;
      os << "Bond " << pt.bond(i)[0]+1 << " - " << pt.bond(i)[1]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }
  
  if(format==3)
    cout << "PERTBANGLE03\n";
  else
    cout << "PERTBONDANGLE\n";
  
  cout << "# number of perturbed angles\n# NTHEG\n"
       << setw(5) << pt.numAngles() << endl;
  cout << "# IBG  JBG  KBG ICBA ICBB\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundangles;
  for(int m=0; m< sys.numMolecules(); ++m){
    AngleIterator ai(sys.mol(m).topology());
    for(;ai;++ai){
      // is this angle perturbed
      int tb=findangle(sys,pt,ai(),m);
      if(tb>=0){
	cout << setw(5) << pt.angle(tb)[0]+1
	     << setw(5) << pt.angle(tb)[1]+1
	     << setw(5) << pt.angle(tb)[2]+1
	     << setw(5) << ai().type()+1
	     << setw(5) << pt.angle(tb).type()+1
	     << endl;
	
	foundangles.insert(tb);
      }
    }
  }
  cout << "END\n";
  
  for(int i=0; i< pt.numAngles();++i)
    if(foundangles.count(i)==0){
      ostringstream os;
      os << "Angle " << pt.angle(i)[0]+1 << " - " << pt.angle(i)[1]+1
	 << " - " << pt.angle(i)[2]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }
  
  if(format==3)
    cout << "PERTIMPDIHEDRAL03\n";
  else
    cout << "PERTIMPROPERDIH\n";
  
  cout << "# number of perturbed improper dihedrals\n# NQHIH\n"
       << setw(5) << pt.numImpropers() << endl;
  cout << "# IBG  JBG  KBG  LBG ICBA ICBB\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> foundimpropers;
  for(int m=0; m<sys.numMolecules(); m++){
    ImproperIterator ii(sys.mol(m).topology());
    for(;ii;++ii){
      // is this angle perturbed
      int tb=findimproper(sys,pt,ii(), m);
      if(tb>=0){
	cout << setw(5) << pt.improper(tb)[0]+1
	     << setw(5) << pt.improper(tb)[1]+1
	     << setw(5) << pt.improper(tb)[2]+1
	     << setw(5) << pt.improper(tb)[3]+1
	     << setw(5) << ii().type()+1
	     << setw(5) << pt.improper(tb).type()+1
	     << endl;
	
	foundimpropers.insert(tb);
      }
    }
  }
  cout << "END\n";
  for(int i=0; i<pt.numImpropers(); ++i)
    if(foundimpropers.count(i)==0){
      ostringstream os;
      os << "Improper " << pt.improper(i)[0]+1 << " - " << pt.improper(i)[1]+1
	 << " - " << pt.improper(i)[2]+1 << " - " << pt.improper(i)[3]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }

  if(format==3)
    cout << "PERTDIHEDRAL03\n";
  else
    cout << "PERTPROPERDIH\n";
  
  cout << "# number of perturbed dihedrals\n#NPHIG\n"
       << setw(5) << pt.numDihedrals() << endl;
  cout << "# IBG  JBG  KBG  LBG ICBA ICBB\n";
  // this may seem a bit clumsy, but we do this by looping over all bonds
  // keep a list of the bonds that we found
  
  set<int> founddihedrals;
  for(int m=0; m<sys.numMolecules(); m++){
    
    DihedralIterator di(sys.mol(m).topology());
    for(;di;++di){
      // is this angle perturbed
      int tb=finddihedral(sys,pt,di(), m);
      if(tb>=0){
	int offset=1;
	findatom(sys,pt,mol,atom,pt.dihedral(tb)[0]);
	for(int m=0;m<mol;++m) offset+=sys.mol(m).numAtoms();
	cout << setw(5) << pt.dihedral(tb)[0]+1
	     << setw(5) << pt.dihedral(tb)[1]+1
	     << setw(5) << pt.dihedral(tb)[2]+1
	     << setw(5) << pt.dihedral(tb)[3]+1
	     << setw(5) << di().type()+1
	     << setw(5) << pt.dihedral(tb).type()+1
	     << endl;
	
	founddihedrals.insert(tb);
      }
    }
  }
  cout << "END\n";

  for(int i=0; i<pt.numDihedrals(); ++i)
    if(founddihedrals.count(i)==0){
      ostringstream os;
      os << "Dihedral " << pt.dihedral(i)[0]+1 << " - " << pt.dihedral(i)[1]+1
	 << " - " << pt.dihedral(i)[2]+1 << " - " << pt.dihedral(i)[3]+1
	 << " not found in topology.";
      throw gromos::Exception("pt_top", os.str());
    }
}



pert::pert readPERTATOM(System &sys, string line, int &start){
  
  istringstream lineStream(line);
  
  // define a few variables
  string nm;
  
  int a, k, l, iiac;
  double dq, dmass;
  
  lineStream >> a; 
  pert pt(a,1);
  
  // read in the perturbation data for the atoms 
  for(int i=0; i<pt.numAtoms(); i++){
    lineStream >> k;
    if(i==0&&start>=0) start-=k;
    pt.setNum(i,k+start);
    lineStream >> l >> nm;
    pt.setName(i, nm);
    for(int j=0; j< pt.numPt(); j++){
      lineStream >> iiac >> dmass >> dq;
      pt.setIac(i,j,iiac-1);
      pt.setMass(i,j,dmass);
      pt.setCharge(i,j,dq);
      lineStream >> l >> l;
    }
  }
  return pt;
}
 
pert readMPERTATOM(System &sys, string line, int &start){
  
  istringstream lineStream(line);

 
   // define a few variables
   string nm;
   
   int a, p, k, iiac;
   double dq;
   
   lineStream >> a >> p; 
   pert pt(a,p);
   
   for(int i=0; i< pt.numPt(); i++){
     lineStream >> nm;
     pt.setPertName(i, nm);
   }
   
    // read in the perturbation data for the atoms (only iac and q)
   for(int i=0; i<pt.numAtoms(); i++){
     lineStream >> k;
     if(i==0&&start>=0) start-=k;
     pt.setNum(i,k+start);
     lineStream >> nm;
     pt.setName(i, nm);
     for(int j=0; j< pt.numPt(); j++){
       lineStream >> iiac >> dq;
       pt.setIac(i,j,iiac-1);
       pt.setCharge(i,j,dq);
     }
   }

   // set the mass in the case of an MPERTATOM block
   for(int i=0; i<pt.numAtoms(); i++){
     int m, a;
     findatom(sys, pt, m, a, i);
     if(a>=0 && m>=0){
       for(int j=0; j< pt.numPt(); j++){
	 pt.setMass(i,j,sys.mol(m).topology().atom(a).mass());
       }
     }
   }
   return pt;
}
 
pert::pert readPERTATOM03(System &sys, string line, int &start){  
  istringstream lineStream(line);
  
  // define a few variables
  string nm;

  int a, k, l, iiac;
  double dq, fdum, dmass;
  
  lineStream >> a; 
  pert pt(a,1);
  
  // read in the perturbation data for the atoms 
  for(int i=0; i<pt.numAtoms(); i++){
    lineStream >> k;
    if(i==0&&start>=0) start-=k;
    pt.setNum(i,k+start);
    lineStream >> l >> nm;
    pt.setName(i, nm);
    for(int j=0; j< pt.numPt(); j++){
      lineStream >>fdum >> fdum >> fdum;
      lineStream >> iiac >> dmass >> dq;
      pt.setIac(i,j,iiac-1);
      pt.setMass(i,j,dmass);
      pt.setCharge(i,j,dq);
      lineStream >> fdum >> fdum;
    }
  }
  return pt;
}

void readPERTBONDSTRETCH(System &sys, pert &pt, string line, int start){
  istringstream lineStream(line);
  int i, j, ta, tb, nb;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){
    lineStream >> i >> j >> ta >> tb;
    pt.addBond(i+start,j+start,tb-1);
  }
}

void readPERTBONDANGLE(System &sys, pert &pt, string line, int start){
  istringstream lineStream(line);
  int i, j,k, ta, tb, nb;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){
    lineStream >> i >> j >> k>> ta >> tb;
    pt.addAngle(i+start,j+start,k+start,tb-1);
  }
}
void readPERTIMPROPERDIH(System &sys, pert &pt, string line, int start){
  istringstream lineStream(line);
  int i, j,k, l, ta, tb, nb;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){
    lineStream >> i >> j >> k>> l>>  ta >> tb;
    pt.addImproper(i+start,j+start,k+start,l+start, tb-1);
  }
}
void readPERTPROPERDIH(System &sys, pert &pt, string line, int start){
  istringstream lineStream(line);
  int i, j,k, l, ta, tb, nb;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){
    lineStream >> i >> j >> k>> l>>  ta >> tb;
    pt.addDihedral(i+start,j+start,k+start,l+start, tb-1);
  }
}

// gromos96 perturbation topology blocks --Clara Christ
void readPERTATOMPAIR(System &sys, pert &pt, string line, int start){
  istringstream lineStream(line);
  int i, j, k, l;
  
  int nb;
  lineStream >> nb;
  for(int ii=0; ii<nb; ++ii){
    lineStream >> i >> j >> k >> l;
    pt.addAtomPair(i+start, j+start, k, l);
  }
} 
void readPERTBOND(GromosForceField &ff, System &sys, pert &pt, string line, int start){
istringstream lineStream(line);
  int i, j, n, nb;
  double fca, b0a, fcb, b0b;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){ 
    lineStream >> i >> j >> n >> fca >> b0a >> fcb >> b0b;
    bool foundBond=false;
    for(int jj=0; jj< ff.numBondTypes(); jj++){
      if(ff.bondType(jj).fc() == fcb && ff.bondType(jj).b0() == b0b){
	pt.addBond(i+start,j+start,jj);
	//cerr << "Adding bond: " << i+start << ", " << j+start << ": " << jj+1 << endl;
	foundBond=true;
      }
    }
    if (foundBond==false){
      ostringstream os;
      os << "No bond type found in topology for bond "  
	 << i << " - " << j 
	 << " in PERTBOND(H) block";
      throw gromos::Exception("pt_top", os.str());
    }
  }
} 
void readPERTBANGLE(GromosForceField &ff, System &sys, pert &pt, string line, int start){
istringstream lineStream(line);
  int i, j, k, nb;
  double fca, t0a, fcb, t0b;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){ 
    lineStream >> i >> j >> k >> fca >> t0a >> fcb >> t0b;
    bool foundAngle=false;
    for(int jj=0; jj< ff.numAngleTypes(); jj++){
      if(ff.angleType(jj).fc() == fcb && ff.angleType(jj).t0() == t0b){
	pt.addAngle(i+start,j+start,k+start,jj);
	//cerr << "Adding angle: " << i+start << ", " << j+start << ", " << k+start << ":  " << jj+1 << endl;
	foundAngle=true;
      }
    }
    if (foundAngle==false){
      ostringstream os;
      os << "No angle type found in topology for angle "  
	 << i << " - " << j << " - " << k
	 << " in PERTBANGLE(H) block";
      throw gromos::Exception("pt_top", os.str());
    }
  }
} 
void readPERTIMPDIHEDRAL(GromosForceField &ff, System &sys, pert &pt, string line, int start){
istringstream lineStream(line);
  int i, j, k, l, nb;
  double fca, q0a, fcb, q0b;
  lineStream >> nb;
  for(int ii=0; ii< nb; ++ii){ 
    lineStream >> i >> j >> k >> l >> fca >> q0a >> fcb >> q0b;
    bool foundImproper=false;
    for(int jj=0; jj< ff.numImproperTypes(); jj++){
      if(ff.improperType(jj).fc() == fcb && ff.improperType(jj).q0() == q0b){
	pt.addImproper(i+start,j+start,k+start,l+start,jj);
	//cerr << "Adding improper: " << i+start << ", " << j+start << ", " << k+start << ", " << l+start << ": " << jj+1 << endl;
	foundImproper=true;
      }
    }
    if (foundImproper==false){
      ostringstream os;
      os << "No improper type found in topology for improper "  
	 << i << " - " << j << " - " << k << " - " << l 
	 << " in PERTIMPDIHEDRAL(H) block";
      throw gromos::Exception("pt_top", os.str());
    }
  }
}
void readPERTDIHEDRAL(GromosForceField &ff, System &sys, pert &pt, string line, int start){
istringstream lineStream(line);
  int i, j, k, l, nb;
  double fca, pda, npa, fcb, pdb, npb;
  lineStream >> nb;
 for(int ii=0; ii< nb; ++ii){ 
   lineStream >> i >> j >> k >> l >> fca >> pda >> npa >> fcb >> pdb >> npb;
    bool foundDihedral=false;
    for(int jj=0; jj< ff.numDihedralTypes(); jj++){
      if(ff.dihedralType(jj).fc() == fcb && ff.dihedralType(jj).pd() == pdb && ff.dihedralType(jj).np() == npb){
	pt.addDihedral(i+start,j+start,k+start,l+start,jj);
	//cerr << "Adding dihedral: " << i+start << ", " << j+start << ", " << k+start << ", " << l+start << ": " << jj+1 << endl;
	foundDihedral=true;
      }
    }
    if (foundDihedral==false){
      ostringstream os;
      os << "No dihedral type found in topology for dihedral "  
	 << i << " - " << j << " - " << k << " - " << l 
	 << " in PERTDIHEDRAL(H) block";
      throw gromos::Exception("pt_top", os.str());
    }
  }
 
}

   


