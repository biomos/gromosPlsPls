// GromosForceField.cc

#include "GromosForceField.h"
#include <vector>
#include <map>
#include "AtomPair.h"
#include "BondType.h"
#include "AngleType.h"
#include "DihedralType.h"
#include "ImproperType.h"
#include "LJType.h"


using namespace std;

namespace gcore{

class GromosForceField_i{
  friend class GromosForceField;
  double d_fpepsi, d_hbar;
  vector<string> d_atomTypeName;
  vector<BondType> d_bondType;
  vector<AngleType> d_angleType;
  vector<DihedralType> d_dihedralType;
  vector<ImproperType> d_improperType;
  map<AtomPair,LJType> d_ljType;
  GromosForceField_i():
    d_fpepsi(0), d_hbar(0), d_atomTypeName(),
    d_bondType(), d_angleType(),
    d_dihedralType(), d_improperType(),
    d_ljType()
  {}
  GromosForceField_i(const GromosForceField_i &gff):
    d_fpepsi(gff.d_fpepsi), d_hbar(gff.d_hbar),
    d_atomTypeName(gff.d_atomTypeName),
    d_bondType(gff.d_bondType), d_angleType(gff.d_angleType),
    d_dihedralType(gff.d_dihedralType), d_improperType(gff.d_improperType),
    d_ljType(gff.d_ljType)
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

void GromosForceField::addAtomTypeName(const string &str)
{d_this->d_atomTypeName.push_back(str);}

void GromosForceField::addBondType(const BondType &b)
{d_this->d_bondType.push_back(b);}

void GromosForceField::addAngleType(const AngleType &b)
{d_this->d_angleType.push_back(b);}

void GromosForceField::addDihedralType(const DihedralType &b)
{d_this->d_dihedralType.push_back(b);}

void GromosForceField::addImproperType(const ImproperType &b)
{d_this->d_improperType.push_back(b);}

void GromosForceField::setLJType(const AtomPair &p, const LJType &l)
{ d_this->d_ljType[p]=l;}

double GromosForceField::fpepsi()const{
  return d_this->d_fpepsi;}

double GromosForceField::hbar()const{
  return d_this->d_hbar;}

int GromosForceField::numAtomTypeNames()const
{ return d_this->d_atomTypeName.size();}

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

}








