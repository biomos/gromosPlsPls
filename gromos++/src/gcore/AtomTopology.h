// gcore_AtomTopology.h
#ifndef INCLUDED_GCORE_ATOMTOPOLOGY
#define INCLUDED_GCORE_ATOMTOPOLOGY

#ifndef INDLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

class AtomTopology_i; 
class Exclusion; 

class AtomTopology {
  AtomTopology_i *d_this;

 public:
  // Constructors
  AtomTopology();
  AtomTopology(const AtomTopology&);

  // Destructor
  ~AtomTopology();

  // Methods
  AtomTopology &operator=(const AtomTopology&);
  void setIac(int);
  void setChargeGroup(int);
  void setCharge(double);
  void setMass(double);
  void setName(const std::string&);
  void setExclusion(const Exclusion &e);
  void setExclusion14(const Exclusion &e);
  void setradius(double);

  // Accessors
  int iac()const;
  int chargeGroup()const;
  double charge()const;
  double mass()const;
  const std::string &name()const;
  const Exclusion &exclusion()const;
  const Exclusion &exclusion14()const;
  double radius()const;
};

}
#endif
