// gio_InTopology.cc

#include "InTopology.h"
#include "Ginstream.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/AngleType.h"
#include "../gcore/Angle.h"
#include "../gcore/Constraint.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Improper.h"
#include "../gcore/LJType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"

#include <map>
#include <deque>
#include <set>

using namespace gcore;
using gio::InTopology_i;
using gio::InTopology;

// Implementation class
class InTopology_i{
  friend class InTopology;
  Ginstream d_gin;
  GromosForceField d_gff;
  System d_sys;
  string d_version;
  string d_name;
  void init();
  InTopology_i (const char *name):
    d_gin(name), d_version(), d_name(name){
    init();
  }
  ~InTopology_i(){
    d_gin.close();
  }
};

// Constructors

InTopology::InTopology(string name):
  d_this(new InTopology_i(name.c_str())){
}

InTopology::~InTopology(){
  delete d_this;
}

const string &InTopology::title()const{
  return d_this->d_gin.title();
}

const string &InTopology::version()const{
  return d_this->d_version;
}

const System &InTopology::system()const{
  return d_this->d_sys;
}

const GromosForceField &InTopology::forceField()const{
  return d_this->d_gff;
}

void InTopology_i::init(){

  if(!d_gin)
    throw InTopology::Exception("Could not open topology file "+d_gin.name()+".");

  // Temporary variables used to read in stuff
  deque<string> resNames;
  map<int,int> resMap;
  deque<AtomTopology> soluteAtoms;
  deque<AtomTopology> solventAtoms;
  set<Bond> bonds;
  set<Angle> angles;
  set<Improper> impropers;
  set<Dihedral> dihedrals;
  set<Constraint> constraints;


  // generic variables
  double d[4];
  int i[5], num;
  string s;

  // Topphyscon block:
  if(! d_gin.check("TOPPHYSCON"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo TOPPHYSCON block!");
  
  d_gin >> d[0] >> d[1];
  d_gff.setFpepsi(d[0]);
  d_gff.setHbar(d[1]);

  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nTOPPHYSCON block is not OK!");

  // Version Block 
  if(! d_gin.check("TOPVERSION"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo TOPVERSION block!");
  d_gin >> d_version;
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nTOPVERSION block is not OK!");

  // AtomTypename
  if(! d_gin.check("ATOMTYPENAME"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo ATOMTYPENAME block!");
  d_gin >> num;
  for(int j = 0; j<num;j++){
    d_gin>>s;
    d_gff.addAtomTypeName(s);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nATOMTYPENAME block is not OK!");

  // RESNAME block
  if(! d_gin.check("RESNAME"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo RESNAME block!");
  d_gin >> num;
  for (int j=0;j<num;j++){
    d_gin >> s;
    resNames.push_back(s);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\n RESNAME block is not OK!");
  

  // SOLUTEATOM block
  if(! d_gin.check("SOLUTEATOM"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo SOLUTEATOM block!");
  d_gin >> num;
  for (int j=0;j<num;j++){
    soluteAtoms.push_back(AtomTopology());
    d_gin>>i[0];
    if(i[0]!=j+1)
      throw InTopology::Exception("Atom numbers are not sequential!");
    
    // residue number
    d_gin >> i[0];
    resMap[j]=--i[0];
    
    // Atom Name
    d_gin >> s;
    soluteAtoms[j].setName(s);

    // IAC
    d_gin >> i[0];
    soluteAtoms[j].setIac(--i[0]);
    
    // mass
    d_gin >> d[0];
    soluteAtoms[j].setMass(d[0]);

    // charge
    d_gin >> d[0];
    soluteAtoms[j].setCharge(d[0]);

    // charge group code
    d_gin >> i[0];
    soluteAtoms[j].setChargeGroup(i[0]);

    // Exclusions 1-2 and 1-3
    Exclusion *e;
    e = new Exclusion();
    d_gin>>i[0];
    for(int l=0;l<i[0];l++){
      d_gin >> i[1];
      e->insert(--i[1]);
    }
    soluteAtoms[j].setExclusion(*e);
    // Exclusions 1-4
    delete e;
    e=new Exclusion();
    d_gin>>i[0];
    for(int l=0;l<i[0];l++){
      d_gin >> i[1];
      e->insert(--i[1]);
    }
    soluteAtoms[j].setExclusion14(*e);
    delete e;
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nSOLUTEATOM block is not OK!");

  // BONDTYPE block
  if(! d_gin.check("BONDTYPE"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BONDTYPE block!");
  d_gin >> num;
  for (int j = 0 ; j < num ; j++){
    d_gin >> d[0] >> d[1];
    d_gff.addBondType(BondType(d[0],d[1]));
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBONDTYPE block is not OK!");

  // BONDH, BOND blocks
  if(! d_gin.check("BONDH"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BONDH block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2];
    Bond bond(--i[0],--i[1]);
    bond.setType(--i[2]);
    bonds.insert(bond);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBONDH block is not OK!");

  if(! d_gin.check("BOND"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BOND block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2];
    Bond bond(--i[0],--i[1]);
    bond.setType(--i[2]);
    bonds.insert(bond);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBOND block is not OK!");

  // BONDANGLETYPE block
  if(! d_gin.check("BONDANGLETYPE"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BONDANGLETYPE block!");
  d_gin >> num;
  for (int j = 0 ; j < num ; j++){
    d_gin >> d[0] >> d[1];
    d_gff.addAngleType(AngleType(d[0],d[1]));
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBONDANGLETYPE block is not OK!");

  // BONDANGLEH, BONDANGLE blocks
  if(! d_gin.check("BONDANGLEH"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BONDANGLEH block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3];
    Angle angle(--i[0],--i[1],--i[2]);
    angle.setType(--i[3]);
    angles.insert(angle);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBONDANGLEH block is not OK!");

  if(! d_gin.check("BONDANGLE"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo BONDANGLE block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3];
    Angle angle(--i[0],--i[1],--i[2]);
    angle.setType(--i[3]);
    angles.insert(angle);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nBONDANGLE block is not OK!");

  // IMPDIHEDRALTYPE block
  if(! d_gin.check("IMPDIHEDRALTYPE"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo IMPDIHEDRALTYPE block!");
  d_gin >> num;
  for (int j = 0 ; j < num ; j++){
    d_gin >> d[0] >> d[1];
    d_gff.addImproperType(ImproperType(d[0],d[1]));
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nIMPDIHEDRALTYPE block is not OK!");

  // IMPDIHEDRALH, IMPDIHEDRAL blocks
  if(! d_gin.check("IMPDIHEDRALH"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo IMPDIHEDRALH block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    Improper improper(--i[0],--i[1],--i[2],--i[3]);
    improper.setType(--i[4]);
    impropers.insert(improper);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nIMPDIHEDRALH block is not OK!");

  if(! d_gin.check("IMPDIHEDRAL"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo IMPDIHEDRAL block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    Improper improper(--i[0],--i[1],--i[2],--i[3]);
    improper.setType(--i[4]);
    impropers.insert(improper);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nIMPDIHEDRAL block is not OK!");
  
  // DIHEDRALTYPE block
  if(! d_gin.check("DIHEDRALTYPE"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo DIHEDRALTYPE block!");
  d_gin >> num;
  for (int j = 0 ; j < num ; j++){
    d_gin >> d[0] >> d[1] >> i[0];
    d_gff.addDihedralType(DihedralType(d[0],d[1],i[0]));
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nDIHEDRALTYPE block is not OK!");

  // DIHEDRALH, DIHEDRAL blocks
  if(! d_gin.check("DIHEDRALH"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo DIHEDRALH block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
    dihedral.setType(--i[4]);
    dihedrals.insert(dihedral);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nDIHEDRALH block is not OK!");

  if(! d_gin.check("DIHEDRAL"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo DIHEDRAL block!");

  d_gin >> num;
  for (int j=0; j<num; j++){
    d_gin >> i[0] >> i[1] >> i[2] >> i[3] >> i[4];
    Dihedral dihedral(--i[0],--i[1],--i[2],--i[3]);
    dihedral.setType(--i[4]);
    dihedrals.insert(dihedral);
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nDIHEDRAL block is not OK!");

  if(! d_gin.check("LJPARAMETERS"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo LJPARAMETERS block!");
  d_gin >> num;
  for (int j = 0 ; j < num ; j++){
    d_gin >> i[0] >> i[1] >> d[0] >> d[1] >> d[2] >> d[3];
    d_gff.setLJType(AtomPair(--i[0],--i[1]),LJType(d[0],d[1],d[2],d[3]));
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\n LJPARAMETERS block is not OK!");
  
  // SOLVENTATOM block --mika
  if(! d_gin.check("SOLVENTATOM"))
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo SOLVENTATOM block!");
  d_gin >> num; // num equals number of atoms per solvent molecule
  for (int j = 0 ; j < num ; j++){
    d_gin >> i[0];
    solventAtoms.push_back(AtomTopology());
        if(i[0]!=j+1)
    throw InTopology::Exception("Solvent Atom numbers are not sequential!");    
	// set name
	d_gin >> s;
    solventAtoms[j].setName(s);
    // set IAC
    d_gin >> i[0];
    solventAtoms[j].setIac(--i[0]);
    // set mass
    d_gin >> d[0];
    solventAtoms[j].setMass(d[0]);
    // set charge
    d_gin >> d[0];
    solventAtoms[j].setCharge(d[0]); 
  }
  if(! d_gin.check())
    throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\n SOLVENTATOM block is not OK!");

  // SOLVENTCONSTR block --mika
    if(! d_gin.check("SOLVENTCONSTR"))
      throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\nNo SOLVENTCONTSTR block!");
   d_gin >> num; // num equals number of constraints
   for (int j = 0 ; j < num ; j++){
      d_gin >> i[0] >> i[1] >> d[0];
      // smack in the constraints
      Constraint constr(--i[0],--i[1]);
      constr.setDist(d[0]);
      constraints.insert(constr);
   }
   if(! d_gin.check())
   throw InTopology::Exception("Topology file "+d_gin.name()+" is corrupted:\n SOLVENTCONSTR block is not OK!");

  // Now parse the stuff into Topologies and the System.

  int last1=0, last=0, lastres=0;
  MoleculeTopology *mt;
  while(soluteAtoms.size()){
    mt=new MoleculeTopology();
    // Detect the last atom of the first molecule & add bonds:
    for(set<Bond>::iterator iter=bonds.begin(),
	  to=bonds.end(); (iter!=to) && (*iter)[0] <=last; ++iter){
      Bond bond=*iter;
      if(bond[1]>last)last=bond[1];
      bond[0]-=last1;
      bond[1]-=last1;
      mt->addBond(bond);
      bonds.erase(iter);
    }
    last++;
    // add Atoms
    
    for(int i=0;i<last-last1; ++i){
      // adapt exclusions:
      Exclusion *e;
      e=new Exclusion();
      for (int l=0;l<soluteAtoms[0].exclusion().size();++l)
	e->insert(soluteAtoms[0].exclusion().atom(l)-last1);
      soluteAtoms[0].setExclusion(*e);
      delete e;
      e=new Exclusion();
      for (int l=0;l<soluteAtoms[0].exclusion14().size();++l)
	e->insert(soluteAtoms[0].exclusion14().atom(l)-last1);
      soluteAtoms[0].setExclusion14(*e);
      delete e;
      
      // now add atom
      mt->addAtom(soluteAtoms[0]);
      soluteAtoms.pop_front();
      int resn=resMap[i+last1]-lastres;
      
      mt->setResNum(i,resn);
      mt->setResName(resn,resNames[resn+lastres]);
    }
    lastres+=mt->numRes();
    // add Angles
    for(set<Angle>::iterator iter=angles.begin(),
	  to=angles.end(); iter != to && (*iter)[0]<last; ++iter){
      Angle angle(*iter);
      angle[0]-=last1;
      angle[1]-=last1;
      angle[2]-=last1;
      mt->addAngle(angle);
      angles.erase(iter);
    }
    // add Dihedrals
    for(set<Dihedral>::iterator iter=dihedrals.begin(),
	  to=dihedrals.end(); iter != to && (*iter)[0]<last; ++iter){
      Dihedral dihedral(*iter);
      dihedral[0]-=last1;
      dihedral[1]-=last1;
      dihedral[2]-=last1;
      dihedral[3]-=last1;
      mt->addDihedral(dihedral);
      dihedrals.erase(iter);
    }
    // add Impropers
    for(set<Improper>::iterator iter=impropers.begin(),
	  to=impropers.end(); iter != to && (*iter)[0]<last; ++iter){
      Improper improper(*iter);
      improper[0]-=last1;
      improper[1]-=last1;
      improper[2]-=last1;
      improper[3]-=last1;
      mt->addImproper(improper);
      impropers.erase(iter);
    }
    d_sys.addMolecule(Molecule(*mt));
    delete mt;
    last1=last;
  }
   
  // add the solvent topology
   SolventTopology *st;

   st=new SolventTopology();
   while (solventAtoms.size()){
   st->addAtom(solventAtoms[0]);
   solventAtoms.pop_front();}

   for(set<Constraint>::iterator iter=constraints.begin(),
     to=constraints.end(); (iter!=to) ; ++iter){

     st->addConstraint(*iter);


    }

    //  lastt++
   // add atom 


   d_sys.addSolvent(Solvent(*st));
   delete st;
}



