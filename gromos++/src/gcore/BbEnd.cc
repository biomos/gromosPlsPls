// gcore_BbEnd.cc

#include "BbEnd.h"
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
using gcore::BbEnd_i;
using gcore::BbEnd;
using gcore::BeBondIt;
using gcore::BeAngleIt;
using gcore::BeDihIt;
using gcore::BeImpIt;
using gcore::BeBondIt_i;
using gcore::BeAngleIt_i;
using gcore::BeDihIt_i;
using gcore::BeImpIt_i;
using gcore::Bond;
using gcore::Angle;
using gcore::Dihedral;
using gcore::Improper;
using gcore::AtomTopology;
using gcore::Exclusion;

class BbEnd_i{

  friend class BbEnd;
  friend class BeBondIt;
  friend class BeAngleIt;
  friend class BeDihIt;
  friend class BeImpIt;

  vector<AtomTopology> d_atoms;
  set<Bond> d_bonds;
  set<Angle> d_angles;
  set<Dihedral> d_dihedrals;
  set<Improper> d_impropers;
  int d_rep;
  string d_resName;
  BbEnd_i():
    d_atoms(),
    d_bonds(),
    d_angles(),
    d_dihedrals(),
    d_impropers(),
    d_rep(),
    d_resName()
  {}
  ~BbEnd_i(){}
};


BbEnd::BbEnd() : 
  d_this(new BbEnd_i())
{}

BbEnd::BbEnd(const BbEnd& mt):
  d_this(new BbEnd_i())
{
  d_this->d_atoms=(mt.d_this->d_atoms);
  d_this->d_bonds=(mt.d_this->d_bonds);
  d_this->d_angles=(mt.d_this->d_angles);
  d_this->d_dihedrals=(mt.d_this->d_dihedrals);
  d_this->d_impropers=(mt.d_this->d_impropers);
  d_this->d_rep=(mt.d_this->d_rep);
  d_this->d_resName=(mt.d_this->d_resName);
}

BbEnd::~BbEnd(){delete d_this;}

// Methods

BbEnd &BbEnd::operator=(const BbEnd &mt){
  if (this != &mt){
    this->BbEnd::~BbEnd();
    new(this) BbEnd(mt);
  }
  return *this;
}

void BbEnd::addAtom(const AtomTopology &a){
  d_this->d_atoms.push_back(a);
  return;
}

void BbEnd::addBond(const Bond &b){
  // add checks if bond there?
  d_this->d_bonds.insert(b);
}

void BbEnd::addAngle(const Angle &a){
  // add checks if angle there?
  d_this->d_angles.insert(a);
}

void BbEnd::addDihedral(const Dihedral &a){
  // add checks if dihedral there?
  d_this->d_dihedrals.insert(a);
}

void BbEnd::addImproper(const Improper &a){
  // add checks if improper there?
  d_this->d_impropers.insert(a);
}

void BbEnd::setResName(const string &s){
  d_this->d_resName=s;
}

void BbEnd::setRep(int i)
{
  d_this->d_rep=i;
}


int BbEnd::numAtoms()const{return d_this->d_atoms.size();}

const AtomTopology &BbEnd::atom(int i)const{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}

const string &BbEnd::resName()const{
  return d_this->d_resName;
}

const int BbEnd::rep()const{return d_this->d_rep;}


class BeBondIt_i{
  friend class BeBondIt;
  set<Bond>::iterator d_it;
  const BbEnd *d_mt;
  // not implemented
  BeBondIt_i(const BeBondIt_i&);
  BeBondIt_i &operator=(const BeBondIt_i &);
public:
  BeBondIt_i():
    d_it(){d_mt=0;}
};

BeBondIt::BeBondIt(const BbEnd &mt):
  d_this(new BeBondIt_i())
{
  d_this->d_it=mt.d_this->d_bonds.begin();
  d_this->d_mt=&mt;
}

BeBondIt::~BeBondIt(){delete d_this;}

void BeBondIt::operator++(){
  ++(d_this->d_it);
}

const Bond &BeBondIt::operator()()const{
  return *(d_this->d_it);
}

BeBondIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_bonds.end();
}

class BeAngleIt_i{
  friend class BeAngleIt;
  set<Angle>::iterator d_it;
  const BbEnd *d_mt;
  // not implemented
  BeAngleIt_i(const BeAngleIt_i&);
  BeAngleIt_i &operator=(const BeAngleIt_i &);
public:
  BeAngleIt_i():
    d_it(){d_mt=0;}
};

BeAngleIt::BeAngleIt(const BbEnd &mt):
  d_this(new BeAngleIt_i())
{
  d_this->d_it=mt.d_this->d_angles.begin();
  d_this->d_mt=&mt;
}

BeAngleIt::~BeAngleIt(){delete d_this;}

void BeAngleIt::operator++(){
  ++(d_this->d_it);
}

const Angle &BeAngleIt::operator()()const{
  return *(d_this->d_it);
}

BeAngleIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_angles.end();
}

class BeImpIt_i{
  friend class BeImpIt;
  set<Improper>::iterator d_it;
  const BbEnd *d_mt;
  // not implemented
  BeImpIt_i(const BeImpIt_i&);
  BeImpIt_i &operator=(const BeImpIt_i &);
public:
  BeImpIt_i():
    d_it(){d_mt=0;}
};

BeImpIt::BeImpIt(const BbEnd &mt):
  d_this(new BeImpIt_i())
{
  d_this->d_it=mt.d_this->d_impropers.begin();
  d_this->d_mt=&mt;
}

BeImpIt::~BeImpIt(){delete d_this;}

void BeImpIt::operator++(){
  ++(d_this->d_it);
}

const Improper &BeImpIt::operator()()const{
  return *(d_this->d_it);
}

BeImpIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_impropers.end();
}

class BeDihIt_i{
  friend class BeDihIt;
  set<Dihedral>::iterator d_it;
  const BbEnd *d_mt;
  // not implemented
  BeDihIt_i(const BeDihIt_i&);
  BeDihIt_i &operator=(const BeDihIt_i &);
public:
  BeDihIt_i():
    d_it(){d_mt=0;}
};

BeDihIt::BeDihIt(const BbEnd &mt):
  d_this(new BeDihIt_i())
{
  d_this->d_it=mt.d_this->d_dihedrals.begin();
  d_this->d_mt=&mt;
}

BeDihIt::~BeDihIt(){delete d_this;}

void BeDihIt::operator++(){
  ++(d_this->d_it);
}

const Dihedral &BeDihIt::operator()()const{
  return *(d_this->d_it);
}

BeDihIt::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_dihedrals.end();
}
 






