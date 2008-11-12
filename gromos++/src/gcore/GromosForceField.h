// gcore_GromosForceField.h

#ifndef INCLUDED_GROMOSFORCEFIELD
#define INCLUDED_GROMOSFORCEFIELD

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

class GromosForceField_i;
class MassType;
class BondType;
class AngleType;
class DihedralType;
class ImproperType;
class LJType;
class CGType;
class AtomPair;
/**
 * Class GromosForceField
 * Purpose: contains all force field parameters (ifp-file in gromos96)
 *
 * Description:
 * The GromosForceField contains all force field parameters that can be 
 * read from a topology. This roughly corresponds to the ifp-file in 
 * gromos96. Exceptions are the values of Fpepsi 
 * (@f$\frac{1}{4\pi\epsilon_0}@f$)and Hbar (@f$\hbar@f$), which find a place
 * in both the GromosForceField and the BuildingBlock classes. This is due
 * to the gromos96 structure
 *
 * @class GromosForceField
 * @author R. Buergi
 * @author gcore
 * @sa gcore::MassType
 * @sa gcore::BondType
 * @sa gcore::AngleType
 * @sa gcore::DihedralType
 * @sa gcore::ImproperType
 * @sa gcore::LJType
 * @sa gcore::BuildingBlock
 */
class GromosForceField{
  GromosForceField_i *d_this;
  // not implemented
  GromosForceField &operator=(const GromosForceField &);

 public:
  /**
   * GromosForceField constructor
   */
  GromosForceField();
  /**
   * GromosForceField copy constructor
   */
  GromosForceField(const GromosForceField &);
  /**
   * GromosForceField deconstructor
   */
  ~GromosForceField();
  // Methods
  /**
   * Method to set the value of Fpepsi  (@f$\frac{1}{4\pi\epsilon_0}@f$)
   */
  void setFpepsi(double fpepsi);
  /**
   * Method to set the value of hbar (@f$\hbar@f$)
   */
  void setHbar(double hbar);
  /**
   * Method to set the value of kB (@f$k_B@f$)
   */
  void setBoltz(double boltz);
  /**
   * Method to set the force field code
   */
  void setForceField(std::string code);
  /**
   * Method to add an Atom Type
   */
  void addAtomTypeName(const std::string &str);
  /**
   * Method to add a Mass Type
   * @param b MassType to add
   */
  void addMassType(const MassType &b);
  /**
   * Method to add a Bond Type
   * @param b BondType to add
   */
  void addBondType(const BondType &b);
  /**
   * Method to add an Angle Type
   * @param b AngleType to add
   */
  void addAngleType(const AngleType &b);
  /**
   * Method to add a Dihedral Type
   * @param b DihedralType to add
   */
  void addDihedralType(const DihedralType &b);
  /**
   * Method to add an Improper Type
   * @param b ImproperType to add
   */
  void addImproperType(const ImproperType &b);
  /**
   * Method to set a Lennard Jones interaction for a specific atom pair
   * @param p An AtomPair defined by their Integer Atom Codes (iac's)
   * @param l The corresponding LJType containing the VDW parameters for 
   *          this AtomPair
   */
  void setLJType(const AtomPair &p, const LJType &l);  
  /**
   * Method to set a coarse grain Lennard Jones interaction for a specific 
   * atom pair
   * @param p An AtomPair defined by their Integer Atom Codes (iac's)
   * @param l The corresponding CGType containing the coarse grain VDW 
   *          parameters for this AtomPair
   */
  void setCGType(const AtomPair &p, const CGType &l);

  // Accessors
  /**
   * Accessor, returns the value of Fpepsi 
   *           ( = @f$\frac{1}{4\pi\epsilon_0}@f$)
   */
  double fpepsi()const;
  /**
   * Accessor, returns the value of Hbar ( = @f$\hbar@f$)
   */
  double hbar()const;
  /**
   * Accessor, returns the value of kB ( = @f$k_B@f$)
   */
  double boltz()const;
  /**
   * Accessor, returns the force field code
   */
  std::string ForceField()const;
  /**
   * Accessor, returns the number of Atom Type Names
   */
  int numAtomTypeNames()const;
  /**
   * Accessor, returns the name of the i-th atom type
   */
  const std::string &atomTypeName(int i) const;
  /** 
   * Accessor, returns the number of MassTypes
   */
  int numMassTypes()const;
  /**
   * Accessor, returns the i-th MassType
   */
  const MassType &massType(int i)const;
  /**
   * Method, returns the Mass from a MassType
   * @param i The gromos96 MassType
   * @return The corresponding Mass
   */
  const double findMass(int i)const;
  /** 
   * Accessor, returns the number of BondTypes
   */
  int numBondTypes()const;
  /**
   * Accessor, returns the i-th BondType
   */
  const BondType &bondType(int i) const;
  /**
   * Accessor, returns the number of AngleTypes
   */
  int numAngleTypes()const;
  /**
   * Accessor, returns the i-th AngleType
   */
  const AngleType &angleType(int i) const;
  /**
   * Accessor, returns the number of DihedralTypes
   */
  int numDihedralTypes()const;
  /** 
   * Accessor, returns the i-th DihedralType
   */
  const DihedralType &dihedralType(int i) const;
  /**
   * Accessor, returns the number of ImproperTypes
   */
  int numImproperTypes()const;
  /**
   * Accessor, returns the i-th ImproperType
   */
  const ImproperType &improperType(int i) const;
  /**
   * Accessor, returns the number of LJTypes
   */
  int numLJTypes()const;
  /**
   * Accessor, returns the LJType for the specified AtomPair
   */
  const LJType &ljType(const AtomPair &p) const;
  /**
   * Accessor, returns the number of CGTypes
   */
  int numCGTypes()const;
  /**
   * Accessor, returns the LJType for the specified AtomPair
   */
  const CGType &cgType(const AtomPair &p) const;
  
};

}
#endif
