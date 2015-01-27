// gcore_MoleculeTopology.cc
#include <cassert>
#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <new>
#include <string>
#include "AtomTopology.h"
#include "AtomPair.h"
#include "Bond.h"
#include "Angle.h"
#include "Improper.h"
#include "Dihedral.h"
#include "CrossDihedral.h"
#include "LJException.h"
#include "MoleculeTopology.h"

using namespace std;
using gcore::MoleculeTopology;
using gcore::BondIterator;
using gcore::BondDipoleIterator;
using gcore::AngleIterator;
using gcore::DihedralIterator;
using gcore::CrossDihedralIterator;
using gcore::ImproperIterator;
using gcore::LJExceptionIterator;
using gcore::Bond;
using gcore::Angle;
using gcore::Dihedral;
using gcore::CrossDihedral;
using gcore::Improper;
using gcore::LJException;
using gcore::AtomTopology;

class gcore::MoleculeTopology_i{

  friend class gcore::MoleculeTopology;
  friend class gcore::BondIterator;
  friend class gcore::BondDipoleIterator;
  friend class gcore::AngleIterator;
  friend class gcore::DihedralIterator;
  friend class gcore::CrossDihedralIterator;
  friend class gcore::ImproperIterator;
  friend class gcore::LJExceptionIterator;

  vector<AtomTopology> d_atoms;
  set<Bond> d_bonds;
  set<Bond> d_dipole_bonds;
  set<Angle> d_angles;
  set<Dihedral> d_dihedrals;
  set<CrossDihedral> d_crossdihedrals;
  set<Improper> d_impropers;
  set<LJException> d_ljexceptions;
  vector<string> d_resNames;
  map<int,int> d_resNums;
  MoleculeTopology_i():
    d_atoms(),
    d_bonds(),
    d_dipole_bonds(),
    d_angles(),
    d_dihedrals(),
    d_crossdihedrals(),
    d_impropers(),
    d_ljexceptions(),
    d_resNames(),
    d_resNums()
  {}
  ~MoleculeTopology_i(){}
};


gcore::MoleculeTopology::MoleculeTopology() : 
  d_this(new MoleculeTopology_i())
{}

MoleculeTopology::MoleculeTopology(const MoleculeTopology& mt):
  d_this(new MoleculeTopology_i())
{
  d_this->d_atoms=(mt.d_this->d_atoms);
  d_this->d_bonds=(mt.d_this->d_bonds);
  d_this->d_dipole_bonds=(mt.d_this->d_dipole_bonds);
  d_this->d_angles=(mt.d_this->d_angles);
  d_this->d_dihedrals=(mt.d_this->d_dihedrals);
  d_this->d_crossdihedrals=(mt.d_this->d_crossdihedrals);
  d_this->d_impropers=(mt.d_this->d_impropers);
  d_this->d_ljexceptions=(mt.d_this->d_ljexceptions);
  d_this->d_resNames=(mt.d_this->d_resNames);
  d_this->d_resNums=(mt.d_this->d_resNums);
}

MoleculeTopology::~MoleculeTopology(){delete d_this;}

// Methods

MoleculeTopology &MoleculeTopology::operator=(const MoleculeTopology &mt){
  if (this != &mt){
    this->MoleculeTopology::~MoleculeTopology();
    new(this) MoleculeTopology(mt);
  }
  return *this;
}

void MoleculeTopology::addAtom(const AtomTopology &a){
  d_this->d_atoms.push_back(a);
  return;
}

void MoleculeTopology::addBond(const Bond &b){
  // add checks if bond there?
  d_this->d_bonds.insert(b);
}

void MoleculeTopology::addDipoleBond(const Bond &b){
  // add checks if bond there?
  d_this->d_dipole_bonds.insert(b);
}

void MoleculeTopology::addAngle(const Angle &a){
  // add checks if angle there?
  d_this->d_angles.insert(a);
}

void MoleculeTopology::addDihedral(const Dihedral &a){
  // add checks if dihedral there?
  d_this->d_dihedrals.insert(a);
}

void MoleculeTopology::addCrossDihedral(const CrossDihedral &a){
  // add checks if crossdihedral there?
  d_this->d_crossdihedrals.insert(a);
}

void MoleculeTopology::addImproper(const Improper &a){
  // add checks if improper there?
  d_this->d_impropers.insert(a);
}

void MoleculeTopology::addLJException(const LJException &lj){
  // add checks if LJ exception there?
  d_this->d_ljexceptions.insert(lj);
}

void MoleculeTopology::setResName(int res, const string &s){
  int num=d_this->d_resNames.size();
  if(res < num){
    d_this->d_resNames[res]=s;
  }
  else{
    for(int i=num; i<res; ++i)
      d_this->d_resNames.push_back(string());
    d_this->d_resNames.push_back(s);
  }
}

void MoleculeTopology::setResNum(int atom, int res){
  d_this->d_resNums[atom]=res;
}

void MoleculeTopology::clearH()
{
  for(unsigned int i=0;i < d_this->d_atoms.size(); i++)
    d_this->d_atoms[i].setH(false);
}

void MoleculeTopology::setHmass(double mass)
{
  for(unsigned int i=0; i< d_this->d_atoms.size(); i++)
    if(d_this->d_atoms[i].mass() - mass < 0.00001)
      d_this->d_atoms[i].setH(true);
}

void MoleculeTopology::setHiac(int iac)
{
   for(unsigned int i=0; i< d_this->d_atoms.size(); i++)
    if(d_this->d_atoms[i].iac() == iac)
      d_this->d_atoms[i].setH(true);
}

int MoleculeTopology::numAtoms()const{return d_this->d_atoms.size();}

int MoleculeTopology::numBonds()const{return d_this->d_bonds.size();}

int MoleculeTopology::numDipoleBonds()const{return d_this->d_dipole_bonds.size();}

int MoleculeTopology::numAngles()const{return d_this->d_angles.size();}

int MoleculeTopology::numImpropers()const{return d_this->d_impropers.size();}

int MoleculeTopology::numDihedrals()const{return d_this->d_dihedrals.size();}

int MoleculeTopology::numCrossDihedrals()const{return d_this->d_crossdihedrals.size();}

int MoleculeTopology::numLJExceptions()const{return d_this->d_ljexceptions.size();}

AtomTopology &MoleculeTopology::atom(int i)
{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}

const AtomTopology &MoleculeTopology::atom(int i)const{
  assert(i < int(d_this->d_atoms.size()));
  return d_this->d_atoms[i];
}

int MoleculeTopology::numRes()const{
  return d_this->d_resNames.size();}

int MoleculeTopology::resNum(int i)const{
  assert(i < int(d_this->d_resNums.size()));
  return d_this->d_resNums[i];
}

const string &MoleculeTopology::resName(int i)const{
  assert(i < int(d_this->d_resNames.size()));
  return d_this->d_resNames[i];
}

class gcore::BondIterator_i{
  friend class gcore::BondIterator;
  set<Bond>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  BondIterator_i(const BondIterator_i&);
  BondIterator_i &operator=(const BondIterator_i &);
public:
  BondIterator_i():
    d_it(){d_mt=0;}
};

gcore::BondIterator::BondIterator(const MoleculeTopology &mt):
  d_this(new BondIterator_i())
{
  d_this->d_it=mt.d_this->d_bonds.begin();
  d_this->d_mt=&mt;
}

BondIterator::~BondIterator(){delete d_this;}

void BondIterator::operator++(){
  ++(d_this->d_it);
}

const Bond &BondIterator::operator()()const{
  return *(d_this->d_it);
}

Bond &BondIterator::operator()(){
  return const_cast<Bond&>(*(d_this->d_it));
}

BondIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_bonds.end();
}

bool BondIterator::last() const{
  // these iterator are not random access iterators and so they don't support
  // end()-1. Thus we have to copy the iterator and advance it an check whether
  // it's at the end.
  set<Bond>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_bonds.end();
}

bool BondIterator::first() const{
  return d_this->d_it == d_this->d_mt->d_this->d_bonds.begin();
}

class gcore::AngleIterator_i{
  friend class gcore::AngleIterator;
  set<Angle>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  AngleIterator_i(const AngleIterator_i&);
  AngleIterator_i &operator=(const AngleIterator_i &);
public:
  AngleIterator_i():
    d_it(){d_mt=0;}
};

gcore::AngleIterator::AngleIterator(const MoleculeTopology &mt):
  d_this(new AngleIterator_i())
{
  d_this->d_it=mt.d_this->d_angles.begin();
  d_this->d_mt=&mt;
}

AngleIterator::~AngleIterator(){delete d_this;}

void AngleIterator::operator++(){
  ++(d_this->d_it);
}

const Angle &AngleIterator::operator()()const{
  return *(d_this->d_it);
}

Angle &AngleIterator::operator()(){
  return const_cast<Angle&>(*(d_this->d_it));
}


AngleIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_angles.end();
}

bool AngleIterator::last() const{
  set<Angle>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_angles.end();
}

bool AngleIterator::first() const{
  return d_this->d_it == d_this->d_mt->d_this->d_angles.begin();
}

class gcore::ImproperIterator_i{
  friend class gcore::ImproperIterator;
  set<Improper>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  ImproperIterator_i(const ImproperIterator_i&);
  ImproperIterator_i &operator=(const ImproperIterator_i &);
public:
  ImproperIterator_i():
    d_it(){d_mt=0;}
};

ImproperIterator::ImproperIterator(const MoleculeTopology &mt):
  d_this(new ImproperIterator_i())
{
  d_this->d_it=mt.d_this->d_impropers.begin();
  d_this->d_mt=&mt;
}

ImproperIterator::~ImproperIterator(){delete d_this;}

void ImproperIterator::operator++(){
  ++(d_this->d_it);
}

const Improper &ImproperIterator::operator()()const{
  return *(d_this->d_it);
}

Improper &ImproperIterator::operator()(){
  return const_cast<Improper&>(*(d_this->d_it));
}

ImproperIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_impropers.end();
}

bool ImproperIterator::last()const{
  set<Improper>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_impropers.end();
}

bool ImproperIterator::first()const{
  return d_this->d_it == d_this->d_mt->d_this->d_impropers.begin();
}

class gcore::DihedralIterator_i{
  friend class gcore::DihedralIterator;
  set<Dihedral>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  DihedralIterator_i(const DihedralIterator_i&);
  DihedralIterator_i &operator=(const DihedralIterator_i &);
public:
  DihedralIterator_i():
    d_it(){d_mt=0;}
};

DihedralIterator::DihedralIterator(const MoleculeTopology &mt):
  d_this(new DihedralIterator_i())
{
  d_this->d_it=mt.d_this->d_dihedrals.begin();
  d_this->d_mt=&mt;
}

DihedralIterator::~DihedralIterator(){delete d_this;}

void DihedralIterator::operator++(){
  ++(d_this->d_it);
}

const Dihedral &DihedralIterator::operator()()const{
  return *(d_this->d_it);
}

Dihedral &DihedralIterator::operator()(){
  return const_cast<Dihedral&>(*(d_this->d_it));
}

DihedralIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_dihedrals.end();
}

bool DihedralIterator::last()const{
  set<Dihedral>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_dihedrals.end();
}

bool DihedralIterator::first()const{
  return d_this->d_it == d_this->d_mt->d_this->d_dihedrals.begin();
}


class gcore::CrossDihedralIterator_i{
  friend class gcore::CrossDihedralIterator;
  set<CrossDihedral>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  CrossDihedralIterator_i(const CrossDihedralIterator_i&);
  CrossDihedralIterator_i &operator=(const CrossDihedralIterator_i &);
public:
  CrossDihedralIterator_i():
    d_it(){d_mt=0;}
};

CrossDihedralIterator::CrossDihedralIterator(const MoleculeTopology &mt):
  d_this(new CrossDihedralIterator_i())
{
  d_this->d_it=mt.d_this->d_crossdihedrals.begin();
  d_this->d_mt=&mt;
}

CrossDihedralIterator::~CrossDihedralIterator(){delete d_this;}

void CrossDihedralIterator::operator++(){
  ++(d_this->d_it);
}

const CrossDihedral &CrossDihedralIterator::operator()()const{
  return *(d_this->d_it);
}

CrossDihedral &CrossDihedralIterator::operator()(){
  return const_cast<CrossDihedral&>(*(d_this->d_it));
}


CrossDihedralIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_crossdihedrals.end();
}

bool CrossDihedralIterator::last()const{
  set<CrossDihedral>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_crossdihedrals.end();
}

bool CrossDihedralIterator::first()const{
  return d_this->d_it == d_this->d_mt->d_this->d_crossdihedrals.begin();
}

class gcore::LJExceptionIterator_i{
  friend class gcore::LJExceptionIterator;
  set<LJException>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  LJExceptionIterator_i(const LJExceptionIterator_i&);
  LJExceptionIterator_i &operator=(const LJExceptionIterator_i &);
public:
  LJExceptionIterator_i():
    d_it(){d_mt=0;}
};

gcore::LJExceptionIterator::LJExceptionIterator(const MoleculeTopology &mt):
  d_this(new LJExceptionIterator_i())
{
  d_this->d_it=mt.d_this->d_ljexceptions.begin();
  d_this->d_mt=&mt;
}

LJExceptionIterator::~LJExceptionIterator(){delete d_this;}

void LJExceptionIterator::operator++(){
  ++(d_this->d_it);
}

const LJException &LJExceptionIterator::operator()()const{
  return *(d_this->d_it);
}

LJException &LJExceptionIterator::operator()(){
  return const_cast<LJException&>(*(d_this->d_it));
}

LJExceptionIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_ljexceptions.end();
}

bool LJExceptionIterator::last() const{
  set<LJException>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_ljexceptions.end();
}

bool LJExceptionIterator::first() const{
  return d_this->d_it == d_this->d_mt->d_this->d_ljexceptions.begin();
}


class gcore::BondDipoleIterator_i{
  friend class gcore::BondDipoleIterator;
  set<Bond>::iterator d_it;
  const MoleculeTopology *d_mt;
  // not implemented
  BondDipoleIterator_i(const BondDipoleIterator_i&);
  BondDipoleIterator_i &operator=(const BondDipoleIterator_i &);
public:
  BondDipoleIterator_i():
    d_it(){d_mt=0;}
};

gcore::BondDipoleIterator::BondDipoleIterator(const MoleculeTopology &mt):
  d_this(new BondDipoleIterator_i())
{
  d_this->d_it=mt.d_this->d_dipole_bonds.begin();
  d_this->d_mt=&mt;
}

BondDipoleIterator::~BondDipoleIterator(){delete d_this;}

void BondDipoleIterator::operator++(){
  ++(d_this->d_it);
}

const Bond &BondDipoleIterator::operator()()const{
  return *(d_this->d_it);
}

Bond &BondDipoleIterator::operator()(){
  return const_cast<Bond&>(*(d_this->d_it));
}

BondDipoleIterator::operator bool()const{
  return d_this->d_it != d_this->d_mt->d_this->d_dipole_bonds.end();
}

bool BondDipoleIterator::last() const{
  // these iterator are not random access iterators and so they don't support
  // end()-1. Thus we have to copy the iterator and advance it an check whether
  // it's at the end.
  set<Bond>::iterator it(d_this->d_it);
  ++it;
  return it == d_this->d_mt->d_this->d_dipole_bonds.end();
}

bool BondDipoleIterator::first() const{
  return d_this->d_it == d_this->d_mt->d_this->d_dipole_bonds.begin();
}
