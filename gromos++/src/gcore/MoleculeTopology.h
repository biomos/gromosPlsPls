// gcore_MoleculeTopology.h

#ifndef INCLUDED_GCORE_MOLECULETOPOLOGY
#define INCLUDED_GCORE_MOLECULETOPOLOGY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class MoleculeTopology_i;
  class GromosForceField;
  class AtomTopology;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BondIterator;
  class AngleIterator;
  class ImproperIterator;
  class DihedralIterator;

  class MoleculeTopology{
    MoleculeTopology_i *d_this;
    // This class contains all topological information

    friend class BondIterator;
    friend class AngleIterator;
    friend class ImproperIterator;
    friend class DihedralIterator;
    
    // not implemented
    MoleculeTopology &operator=(const MoleculeTopology &);

  public:
    MoleculeTopology();
    MoleculeTopology(const MoleculeTopology &);
    ~MoleculeTopology();
    
    // Methods
    void addAtom(const AtomTopology &a);
    void addBond(const Bond &b);
    void addAngle(const Angle &b);
    void addDihedral(const Dihedral &b);
    void addImproper(const Improper &b);
    void setResName(int res, const std::string &s);
    void setResNum(int atom, int res);
    
    // Accessors
    int numAtoms()const;
    const AtomTopology& atom(int i) const; 
    int numRes()const;
    int resNum(int atom) const;
    const std::string &resName(int i)const;
    
  }; /* class MoleculeTopology */


  class BondIterator_i;
  class AngleIterator_i;
  class ImproperIterator_i;
  class DihedralIterator_i;

  class BondIterator{
    BondIterator_i *d_this;
    // not implemented
    BondIterator();
    BondIterator(const BondIterator&);
    BondIterator &operator=(const BondIterator &);
  public:
    BondIterator(const MoleculeTopology &mt);
    ~BondIterator();
    void operator++();
    const Bond &operator()()const;
    operator bool()const;
  };

  class AngleIterator{
    AngleIterator_i *d_this;
    // not implemented
    AngleIterator();
    AngleIterator(const AngleIterator&);
    AngleIterator &operator=(const AngleIterator &);
  public:
    AngleIterator(const MoleculeTopology &mt);
    ~AngleIterator();
    void operator++();
    const Angle &operator()()const;
    operator bool()const;
  };

  class ImproperIterator{
    ImproperIterator_i *d_this;
    // not implemented
    ImproperIterator();
    ImproperIterator(const ImproperIterator&);
    ImproperIterator &operator=(const ImproperIterator &);
  public:
    ImproperIterator(const MoleculeTopology &mt);
    ~ImproperIterator();
    void operator++();
    const Improper &operator()()const;
    operator bool()const;
  };

  class DihedralIterator{
    DihedralIterator_i *d_this;
    // not implemented
    DihedralIterator();
    DihedralIterator(const DihedralIterator&);
    DihedralIterator &operator=(const DihedralIterator &);
  public:
    DihedralIterator(const MoleculeTopology &mt);
    ~DihedralIterator();
    void operator++();
    const Dihedral &operator()()const;
    operator bool()const;
  };

} /* Namespace */ 
#endif



