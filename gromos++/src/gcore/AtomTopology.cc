//gcore_AtomTopology.cc
#include "AtomTopology.h"
#include "Exclusion.h"
#include <new>

using gcore::AtomTopology;
using gcore::AtomTopology_i;
using gcore::Exclusion;
using namespace std;

class AtomTopology_i {
  friend class AtomTopology;
  int d_iac;
  int d_chGrp;
  double d_charge;
  double d_mass;
  string d_name;
  Exclusion d_excl;
  Exclusion d_excl14;
  double d_radius;
  bool d_isH;

 public:
  // Constructors
  AtomTopology_i();
  ~AtomTopology_i(){}
};

AtomTopology_i::AtomTopology_i(): 
  d_iac(-1),
  d_chGrp(-1),
  d_charge(0), 
  d_mass(0), 
  d_name(""), 
  d_excl(),
  d_excl14(),
  d_radius(0),
  d_isH(false)
{}


AtomTopology::AtomTopology(): d_this(new AtomTopology_i()){}


AtomTopology::AtomTopology(const AtomTopology &at){
  d_this=new AtomTopology_i();
  d_this->d_iac=at.d_this->d_iac;
  d_this->d_chGrp=at.d_this->d_chGrp;
  d_this->d_charge=at.d_this->d_charge;
  d_this->d_mass=at.d_this->d_mass;
  d_this->d_name=at.d_this->d_name;
  d_this->d_excl=at.d_this->d_excl;
  d_this->d_excl14=at.d_this->d_excl14;
  d_this->d_radius=at.d_this->d_radius;
  d_this->d_isH=at.d_this->d_isH;
}

AtomTopology::~AtomTopology(){delete d_this;}

AtomTopology &AtomTopology::operator=(const AtomTopology &at){
  if (this != &at){
    this->AtomTopology::~AtomTopology();
    new(this) AtomTopology(at);
  }
  return *this;
} 

void AtomTopology::setIac(int i){d_this->d_iac=i;}
void AtomTopology::setChargeGroup(int i){d_this->d_chGrp=i;}
void AtomTopology::setCharge(double d){d_this->d_charge = d;}
void AtomTopology::setMass(double d){d_this->d_mass = d;}
void AtomTopology::setName(const string& s){d_this->d_name = s;}
void AtomTopology::setExclusion(const Exclusion &e){d_this->d_excl=e;}
void AtomTopology::setExclusion14(const Exclusion &e){d_this->d_excl14=e;}
void AtomTopology::setradius(double d){d_this->d_radius = d;}
void AtomTopology::setH(bool b){d_this->d_isH = b;}

int AtomTopology::iac()const{return d_this->d_iac;}
int AtomTopology::chargeGroup()const{return d_this->d_chGrp;}
double AtomTopology::charge()const{return d_this->d_charge;}
double AtomTopology::mass()const{return d_this->d_mass;}
const string &AtomTopology::name()const{return d_this->d_name;}
const Exclusion &AtomTopology::exclusion()const{return d_this->d_excl;}
const Exclusion &AtomTopology::exclusion14()const{return d_this->d_excl14;}
double AtomTopology::radius()const{return d_this->d_radius;}
const bool AtomTopology::isH()const{return d_this->d_isH;}
