// GromosForceField.cc

#include "GromosForceField.h"
#include <vector>
#include <map>
#include "AtomPair.h"
#include "MassType.h"
#include "BondType.h"
#include "AngleType.h"
#include "DihedralType.h"
#include "ImproperType.h"
#include "LJType.h"
#include "CGType.h"
#include "LJExcType.h"


using namespace std;

namespace gcore{

class GromosForceField_i{
  friend class GromosForceField;
  double d_fpepsi, d_hbar, d_spdl, d_boltz;
  std::string d_ffcode;
  vector<string> d_atomTypeName;
  map<int, MassType> d_massType;
  map<int, BondType> d_bondType;
  map<int, AngleType> d_angleType;
  map<int, DihedralType> d_dihedralType;
  map<int, ImproperType> d_improperType;
  map<AtomPair,LJType> d_ljType;
  map<AtomPair,LJExcType> d_ljexcType;
  map<AtomPair,CGType> d_cgType;
  GromosForceField_i():
    d_fpepsi(0), d_hbar(0), d_spdl(0), d_boltz(0), d_ffcode("_no_FORCEFIELD_block_given_"),
    d_atomTypeName(), d_massType(), d_bondType(), d_angleType(),
    d_dihedralType(), d_improperType(), d_ljType(), d_ljexcType(), d_cgType()
  {}
  GromosForceField_i(const GromosForceField_i &gff):
    d_fpepsi(gff.d_fpepsi), d_hbar(gff.d_hbar), d_spdl(gff.d_spdl),
    d_boltz(gff.d_boltz), d_ffcode(gff.d_ffcode),
    d_atomTypeName(gff.d_atomTypeName), d_massType(gff.d_massType),
    d_bondType(gff.d_bondType), d_angleType(gff.d_angleType),
    d_dihedralType(gff.d_dihedralType), d_improperType(gff.d_improperType),
    d_ljType(gff.d_ljType), d_ljexcType(gff.d_ljexcType), d_cgType(gff.d_cgType)
  {}

  ~GromosForceField_i(){}
};

GromosForceField::GromosForceField(): d_this(new GromosForceField_i()){}

  GromosForceField::GromosForceField ( const GromosForceField &gff):
    d_this(new GromosForceField_i(*gff.d_this)){}


GromosForceField::~GromosForceField(){delete d_this;}

void GromosForceField::setFpepsi(double fpepsi)
{d_this->d_fpepsi=fpepsi;}

void GromosForceField::setHbar(double hbar)
{d_this->d_hbar=hbar;}

void GromosForceField::setSpdl(double spdl)
{d_this->d_spdl=spdl;}

void GromosForceField::setBoltz(double boltz)
{d_this->d_boltz=boltz;}

void GromosForceField::setForceField(const string str)
{d_this->d_ffcode=str;}

void GromosForceField::addAtomTypeName(const string &str)
{d_this->d_atomTypeName.push_back(str);}

void GromosForceField::addMassType(const MassType &b)
{d_this->d_massType[b.n()] = b;}

void GromosForceField::addBondType(const BondType &b)
{d_this->d_bondType[b.code()] = b;}

void GromosForceField::addAngleType(const AngleType &b)
{d_this->d_angleType[b.code()] = b;}

void GromosForceField::addDihedralType(const DihedralType &b)
{d_this->d_dihedralType[b.code()] = b;}

void GromosForceField::addImproperType(const ImproperType &b)
{d_this->d_improperType[b.code()] = b;}

void GromosForceField::setLJType(const AtomPair &p, const LJType &l)
{ d_this->d_ljType[p]=l;}

void GromosForceField::setLJExcType(const AtomPair &p, const LJExcType &l)
{ 
  d_this->d_ljexcType[p]=l;
}

void GromosForceField::setCGType(const AtomPair &p, const CGType &l)
{ d_this->d_cgType[p]=l;}

double GromosForceField::fpepsi()const{
  return d_this->d_fpepsi;}

double GromosForceField::hbar()const{
  return d_this->d_hbar;}

double GromosForceField::spdl() const{
  return d_this->d_spdl;}

double GromosForceField::boltz()const{
  return d_this->d_boltz;}

std::string GromosForceField::ForceField()const{
  return d_this->d_ffcode;}

int GromosForceField::numAtomTypeNames()const
{ return d_this->d_atomTypeName.size();}

int GromosForceField::numMassTypes()const
{return d_this->d_massType.size();}

int GromosForceField::numBondTypes()const
{ return d_this->d_bondType.size();}

int GromosForceField::numAngleTypes()const
{ return d_this->d_angleType.size();}

int GromosForceField::numImproperTypes()const
{ return d_this->d_improperType.size();}

int GromosForceField::numDihedralTypes()const
{ return d_this->d_dihedralType.size();}

int GromosForceField::numLJTypes()const
{ return d_this->d_ljType.size();}

int GromosForceField::numLJExcTypes()const
{ return d_this->d_ljexcType.size();}

int GromosForceField::numCGTypes()const
{ return d_this->d_cgType.size();}

const MassType &GromosForceField::massType(const int i) const
{return d_this->d_massType[i];}

const double GromosForceField::findMass(const int i)const
{
  for(unsigned int k=0; k<d_this->d_massType.size(); k++)
    if(d_this->d_massType[k].n()==i) return d_this->d_massType[k].am();
  return 0.0;
}

const BondType &GromosForceField::bondType(const int i) const
{ return d_this->d_bondType[i];}

const AngleType &GromosForceField::angleType(const int i) const
{ return d_this->d_angleType[i];}

const DihedralType &GromosForceField::dihedralType(const int i) const
{ return d_this->d_dihedralType[i];}

const ImproperType &GromosForceField::improperType(const int i) const
{ return d_this->d_improperType[i];}

const string &GromosForceField::atomTypeName(const int i) const 
{ return d_this->d_atomTypeName[i];}

const LJType &GromosForceField::ljType(const AtomPair &p) const
{return d_this->d_ljType[p]; }

const LJExcType &GromosForceField::ljexcType(const AtomPair &p) const
{return d_this->d_ljexcType[p]; }

map<AtomPair,LJExcType> &GromosForceField::ljexceptions() const
{return d_this->d_ljexcType; }

const CGType &GromosForceField::cgType(const AtomPair &p) const
{return d_this->d_cgType[p];
  }

  int GromosForceField::dummyAtomType() const {
    // first try to find the atom type named DUM
    int dummyAtomType = -1;
    for (int i = 0; i < numAtomTypeNames(); ++i) {
      if (atomTypeName(i) == "DUM")
        dummyAtomType = i;
    }

    if (dummyAtomType == -1)
      return -1;

    // then check whether this is really a dummy
    bool isDummy = true;
    for (int i = 0; i < numAngleTypes(); ++i) {
      const LJType lj = ljType(AtomPair(i, dummyAtomType));
      if (lj.c12() != 0.0 && lj.c6() != 0.0 && lj.cs12() != 0.0 && lj.cs6() != 0.0) {
        isDummy = false;
        break;
      } // if is dummy
    } // for atoms

    if (isDummy)
      return dummyAtomType;

    return -1;
  }

}
