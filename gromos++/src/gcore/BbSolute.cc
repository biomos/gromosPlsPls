// gcore_BbSolute.cc

#include "BbSolute.h"
#include "AtomTopology.h"
#include "Exclusion.h"
#include "Bond.h"
#include "Angle.h"
#include "Improper.h"
#include "Dihedral.h"
#include <set>
#include <vector>
#include <map>
#include <new>

using namespace std;
using gcore::BbSolute_i;
using gcore::BbSolute;
using gcore::BbBondIt;
using gcore::BbAngleIt;
using gcore::BbDihIt;
using gcore::BbImpIt;
using gcore::BbBondIt_i;
using gcore::BbAngleIt_i;
using gcore::BbDihIt_i;
using gcore::BbImpIt_i;
using gcore::Bond;
using gcore::Angle;
using gcore::Dihedral;
using gcore::Improper;
using gcore::AtomTopology;
using gcore::Exclusion;

class BbSolute_i{

  friend class BbSolute;
  friend class BbBondIt;
  friend class BbAngleIt;
  friend class BbDihIt;
  friend class BbImpIt;

  vector<AtomTopology> d_atoms;
  vector<Exclusion> d_pexcl;
  set<Bond> d_bonds;
  set<Angle> d_angles;
  set<Dihedral> d_dihedrals;
  set<Improper> d_impropers;
  string d_resName;
  BbSolute_i():
    d_atoms(),
    d_pexcl(),
    d_bonds(),
    d_angles(),
    d_dihedrals(),
    d_impropers(),
    d_resName()
  {}
  ~BbSolute_i(){}
};


BbSolute::BbSolute() : 
  d_this(new BbSolute_i())
{}

BbSolute::BbSolute(const BbSolute& mt):
  d_this(new BbSolute_i())
{
  d_this->d_atoms=(mt.d_this->d_atoms);
  d_this->d_pexcl=(mt.d_this->d_pexcl);
  d_this->d_bonds=(mt.d_this->d_bonds);
  d_this->d_angles=(mt.d_this->d_angles);
  d_this->d_dihedrals=(mt.d_this->d_dihedrals);
  d_this->d_impropers=(mt.d_this->d_impropers);
  d_this->d_resName=(mt.d_this->d_resName);
}

BbSolute::~BbSolute(){delete d_this;}

// Methods

BbSolute &BbSolute::operator=(const BbSolute &mt){
  if (this != &mt){
    this->BbSolute::~BbSolute();
    new(this) BbSolute(mt);
  }
  return *this;
}

void BbSolute::addAtom(const AtomTopology &a){
  d_this->d_atoms.push_back(a);
  return;
}

void BbSolute::addPexcl(const Exclusion &a)
{
  d_this->d_pexcl.push_back(a);
  return;
}

void BbSolute::addBond(const Bond &b){
  // add checks if bond there?
  d_this->d_bonds.insert(b);
}

void BbSolute::addAngle(const Angle &a){
  // add checks if angle there?
  d_this->d_angles.insert(a);
}

void BbSolute::addDihedral(const Dihedral &a){
  // add checks if dihedral there?
  d_this->d_dihedrals.insert(a);
}

void BbSolute::addImproper(const Improper &a){
  // add checks if improper there?
  d_this->d_impropers.insert(a);
}

void BbSolute::setResName(const string &s){
  d_this->d_resName=s;
}

int BbSolute::numAtoms()const{return d_this->d_atoms.size();}

int BbSolute::numPexcl()const{return d_this->d_pexcl.size();}

const AtomTopology &BbSolute::atom(int i)const{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}

const Exclusion &BbSolute::pexcl(int i)const
{
  assert(i < int(d_this->d_pexcl.size()));
  return d_this->d_pexcl[i];
}

const string &BbSolute::resName()const{
  return d_this->d_resName;
}

class BbBondIt_i{
  friend class BbBondIt;
  set<Bond>::iterator d_it;
  const BbSolute *d_mt;
  // not implemented
  BbBondIt_i(const BbBondIt_i&);
  BbBondIt_i &operator=(const BbBondIt_i &);
public:
  BbBondIt_i():
    d_it(){d_mt=0;}
};

BbBondIt::BbBondIt(const BbSolute &mt):
  d_this(new BbBondIt_i())
{
  d_this->d_it=mt.d_this->d_bonds.begin();
  d_this->d_mt=&mt;
}

BbBondIt::~BbBondIt(){delete d_this;}

void BbBondIt::operator++(){
  ++(d_this->d_it);
}

const Bond &BbBondIt::operator()()const{
  return *(d_this->d_it);
}

BbBondIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_bonds.end();
}

class BbAngleIt_i{
  friend class BbAngleIt;
  set<Angle>::iterator d_it;
  const BbSolute *d_mt;
  // not implemented
  BbAngleIt_i(const BbAngleIt_i&);
  BbAngleIt_i &operator=(const BbAngleIt_i &);
public:
  BbAngleIt_i():
    d_it(){d_mt=0;}
};

BbAngleIt::BbAngleIt(const BbSolute &mt):
  d_this(new BbAngleIt_i())
{
  d_this->d_it=mt.d_this->d_angles.begin();
  d_this->d_mt=&mt;
}

BbAngleIt::~BbAngleIt(){delete d_this;}

void BbAngleIt::operator++(){
  ++(d_this->d_it);
}

const Angle &BbAngleIt::operator()()const{
  return *(d_this->d_it);
}

BbAngleIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_angles.end();
}

class BbImpIt_i{
  friend class BbImpIt;
  set<Improper>::iterator d_it;
  const BbSolute *d_mt;
  // not implemented
  BbImpIt_i(const BbImpIt_i&);
  BbImpIt_i &operator=(const BbImpIt_i &);
public:
  BbImpIt_i():
    d_it(){d_mt=0;}
};

BbImpIt::BbImpIt(const BbSolute &mt):
  d_this(new BbImpIt_i())
{
  d_this->d_it=mt.d_this->d_impropers.begin();
  d_this->d_mt=&mt;
}

BbImpIt::~BbImpIt(){delete d_this;}

void BbImpIt::operator++(){
  ++(d_this->d_it);
}

const Improper &BbImpIt::operator()()const{
  return *(d_this->d_it);
}

BbImpIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_impropers.end();
}

class BbDihIt_i{
  friend class BbDihIt;
  set<Dihedral>::iterator d_it;
  const BbSolute *d_mt;
  // not implemented
  BbDihIt_i(const BbDihIt_i&);
  BbDihIt_i &operator=(const BbDihIt_i &);
public:
  BbDihIt_i():
    d_it(){d_mt=0;}
};

BbDihIt::BbDihIt(const BbSolute &mt):
  d_this(new BbDihIt_i())
{
  d_this->d_it=mt.d_this->d_dihedrals.begin();
  d_this->d_mt=&mt;
}

BbDihIt::~BbDihIt(){delete d_this;}

void BbDihIt::operator++(){
  ++(d_this->d_it);
}

const Dihedral &BbDihIt::operator()()const{
  return *(d_this->d_it);
}

BbDihIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_dihedrals.end();
}
 






