//ener calculates (non-bonded) interaction energies for specific atoms

#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include "../src/gmath/Vec.h"
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/SimplePairlist.h"
#include "../src/utils/Energy.h"


using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;

 class pert
{
  std::vector<int> d_num;
  std::vector<string> d_names;
  std::vector<string> d_pert;
  std::vector< std::vector <double> > d_charge;
  std::vector< std::vector <int> > d_iac;
  
public:
  pert(int a, int p);
  ~pert(){};
  void setIac(int a, int p, int iac);
  void setCharge(int a, int p, double q);
  void setName(int a, string name);
  void setPertName(int p, string name);
  
  void setNum(int a, int num);
  std::vector<double> charge(int p){return d_charge[p];}
  std::vector<int> iac(int p){return d_iac[p];}
  int iac(int a, int p){return d_iac[p][a];}
  double charge(int a, int p){return d_charge[p][a];}
  string name(int a){return d_names[a];}
  string pertName(int p){return d_pert[p];}
  int num(int a){return d_num[a];}
  int numPt(){return d_iac.size();}
  int numAtoms(){return d_names.size();}
  
};

void findatom(System &sys, pert &pt, int &mol, int &atom, int counter);
void modtopo(System &sys, pert &pt, int iipt);
 
int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "props" << "time" << "cut" << "eps"
         << "kap" << "soft" << "softpar" << "traj" << "firstatom" << "pttopo";

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <atomspecifier>\n";
  usage += "\t@props <propertyspecifier>\n";
  usage += "\t@time <time> <dt>\n";
  usage += "\t@cut <cut-off distance>\n";
  usage += "\t@eps <epsilon for reaction field correction>\n";
  usage += "\t@kap <kappa for reaction field correction>\n";
  usage += "\t@soft <atom specifier for soft atoms>\n";
  usage += "\t@softpar <lam> <a_lj> <a_c>\n";
  usage += "\t@traj  <trajectory files>\n";
  usage += "\t@firstatom <first atom>\n";
  usage += "\t@pttopo <perturbation topology>\n";
  
  
 
try{
  Arguments args(argc, argv, knowns, usage);

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
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());


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
  double dq, fdum;
  
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
      lineStream >> iiac;
      if(spt) lineStream >> fdum;
      lineStream >> dq;
      pt.setIac(i,j,iiac-1);
      pt.setCharge(i,j,dq);
      if(spt) lineStream >> l >> l;
    }
  }
  

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
  en.setAtoms(atoms);
  
  // set properties
  PropertyContainer props(sys, pbc);
  {
    Arguments::const_iterator iter=args.lower_bound("props");
    Arguments::const_iterator to=args.upper_bound("props");
    for(;iter!=to;iter++){
      string p=iter->second.c_str();
      props.addSpecifier(p);
    }
  }
  en.setProperties(props);

  // set non-bonded parameters
  //   get cut-off distance
  {
    Arguments::const_iterator iter=args.lower_bound("cut");
    if(iter!=args.upper_bound("cut"))
      en.setCutOff(atof(iter->second.c_str()));
  }
  //  get epsilon and kappa
  {
    double eps=0.0, kap=0.0;
    Arguments::const_iterator iter=args.lower_bound("eps");
    if(iter!=args.upper_bound("eps"))
      eps=atof(iter->second.c_str());
    iter=args.lower_bound("kap");
    if(iter!=args.upper_bound("kap"))
      kap=atof(iter->second.c_str());
    en.setRF(eps, kap);
  }
  // get soft atom list
  AtomSpecifier soft(sys);
  {
    int lsoft=0;
    Arguments::const_iterator iter=args.lower_bound("soft");
    Arguments::const_iterator to=args.upper_bound("soft");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      soft.addSpecifier(spec);
      lsoft=1;
    }
    //  get al2
    double lam=0, alj=0, a_c=0;
    iter=args.lower_bound("softpar");
    if(iter!=args.upper_bound("softpar")){
      lam=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("softpar")){
      alj=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("softpar"))
      a_c=atof(iter->second.c_str());
    else if(lsoft)
      throw gromos::Exception("Ener", 
	 "soft atoms indicated, but not all parameters defined.\n");
    
    en.setSoft(soft, lam, alj, a_c);
  }
 
  // define input coordinate
  InG96 ic;
  
  // declare some variables for averaging
  int num_frames=0;
  vector<double> nb(pt.numPt(),0.0);

  // open files for output
  ofstream fout[pt.numPt()];
  cout << "opened " << pt.numPt() << " files" << endl;
  
  for(int i=0; i< pt.numPt(); ++i){
    string name="en_"+pt.pertName(i)+".dat";
    cout << "  " << name << endl;
    
    fout[i].open(name.c_str());
    
    fout[i] << "# Time"
	    << "           vanderwaals"
	    << "         electrostatic"
	    << "            non-bonded"
	    << endl;
    fout[i].precision(10);
    fout[i].setf(ios::right, ios::adjustfield);
  }
  
  // loop over all trajectories
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
    ic.open((iter->second).c_str());
    ic.select("ALL");
    
      // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      // we have to gather with any method to get covalent interactions 
      // and charge-groups connected
      pbc->gathergr();

      // make the pairlists
      en.calcPairlist();
      
      // loop over the different perturbations
      for(int p=0; p<pt.numPt(); ++p){
	// set the parameters
	modtopo(sys, pt, p);
	// calculate the interactions
	en.calcNb_interactions();
	
	nb[p] += en.nb();
	
	fout[p] << setw(6) << time << setw(22) << en.vdw()
		<< setw(22) << en.el() << setw(22) << en.nb() << endl;
      }
      cout << time << endl;
      
      time+=dt;
      num_frames++;
    }
  }
  // print out averages
  if(num_frames>1){
    for(int p=0; p< pt.numPt(); ++p){
      fout[p] << "# ave " << setw(36) << nb[p]/num_frames << endl;
      fout[p].close();
    }
  }
}
 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
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




void modtopo(System &sys, pert &pt, int iipt)
{
  int counter=0;
  int mol, atom;
  
  //determine molecule and atom for this counter
  findatom(sys, pt, mol, atom, counter);
  
  // loop over the molecules 
  for(int m=0;m<sys.numMolecules();m++){
    
    if(!(m<mol || counter==pt.numAtoms())){
      
      // loop over the atoms in this molecule
      for(int aa=0;aa<sys.mol(m).topology().numAtoms();aa++){
	if(aa==atom){
	  if(pt.name(counter)!=sys.mol(m).topology().atom(aa).name()){
	    ostringstream os;
	    os << "Atom names in (perturbation) topologies do not match\n"
	       << "Topology: " << sys.mol(m).topology().atom(aa).name() 
	       << " (" << m+1 << ":" << aa+1 << ")"
	       << "\tPerturbation topology: " << pt.name(counter)
	       << " (" << counter << ")";
	    
	    
	    throw gromos::Exception("m_ener", os.str());
	  }
	  
	  sys.mol(m).topology().atom(aa).setIac(pt.iac(counter, iipt));
	  sys.mol(m).topology().atom(aa).setCharge(pt.charge(counter, iipt));
	  
	  // loop to next atom in the perturbation list
	  counter++;
	  findatom(sys, pt, mol, atom, counter);
	}
      }
      
    }
  }
}
    

pert::pert(int a, int p)
{
  for(int i=0;i<p;i++){
    
    vector<int> iac;
    vector<double> charge;
  
    for(int j=0;j<a;j++){
      iac.push_back(0);
      charge.push_back(0.0);
    }
    d_iac.push_back(iac);
    d_charge.push_back(charge);
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
