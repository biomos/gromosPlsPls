// gcore_GromosForceField.h

#ifndef INCLUDED_GROMOSFORCEFIELD
#define INCLUDED_GROMOSFORCEFIELD

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

class GromosForceField_i;
class BondType;
class AngleType;
class DihedralType;
class ImproperType;
class LJType;
class AtomPair;

class GromosForceField{
  GromosForceField_i *d_this;
  // not implemented
  GromosForceField &operator=(const GromosForceField &);

 public:
  GromosForceField();
  GromosForceField(const GromosForceField &);
  ~GromosForceField();
  // Methods
  void setFpepsi(double fpepsi);
  void setHbar(double hbar);
  void addAtomTypeName(const std::string &str);
  void addBondType(const BondType &b);
  void addAngleType(const AngleType &b);
  void addDihedralType(const DihedralType &b);
  void addImproperType(const ImproperType &b);
  void setLJType(const AtomPair &p, const LJType &l);
  // Accessors
  double fpepsi()const;
  double hbar()const;
  int numAtomTypeNames()const;
  const std::string &atomTypeName(int i) const;
  int numBondTypes()const;
  const BondType &bondType(int i) const;
  int numAngleTypes()const;
  const AngleType &angleType(int i) const;
  int numDihedralTypes()const;
  const DihedralType &dihedralType(int i) const;
  int numImproperTypes()const;
  const ImproperType &improperType(int i) const;
  int numLJTypes()const;
  const LJType &ljType(const AtomPair &p) const;
};

}
#endif
