// gcore_BbSolute.h

#ifndef INCLUDED_GCORE_BBSOLUTE
#define INCLUDED_GCORE_BBSOLUTE

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class BbSolute_i;
  class GromosForceField;
  class AtomTopology;
  class Exclusion;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BbBondIt;
  class BbAngleIt;
  class BbImpIt;
  class BbDihIt;

  class BbSolute{
    BbSolute_i *d_this;
    // This class contains all topological information

    friend class BbBondIt;
    friend class BbAngleIt;
    friend class BbImpIt;
    friend class BbDihIt;
    
    // not implemented
    BbSolute &operator=(const BbSolute &);

  public:
    BbSolute();
    BbSolute(const BbSolute &);
    ~BbSolute();
    
    // Methods
    void addAtom(const AtomTopology &a);
    void addPexcl(const Exclusion &a);
    void addBond(const Bond &b);
    void addAngle(const Angle &b);
    void addDihedral(const Dihedral &b);
    void addImproper(const Improper &b);
    void setResName(const std::string &s);
    
    // Accessors
    int numAtoms()const;
    int numPexcl()const;
    const AtomTopology& atom(int i) const; 
    const Exclusion& pexcl(int i) const;
    const std::string &resName()const;
    
  }; /* class BbSolute */

  
  class BbBondIt_i;
  class BbAngleIt_i;
  class BbImpIt_i;
  class BbDihIt_i;

  class BbBondIt{
    BbBondIt_i *d_this;
    // not implemented
    BbBondIt();
    BbBondIt(const BbBondIt&);
    BbBondIt &operator=(const BbBondIt &);
  public:
    BbBondIt(const BbSolute &mt);
    ~BbBondIt();
    void operator++();
    const Bond &operator()()const;
    operator bool()const;
  };

  class BbAngleIt{
    BbAngleIt_i *d_this;
    // not implemented
    BbAngleIt();
    BbAngleIt(const BbAngleIt&);
    BbAngleIt &operator=(const BbAngleIt &);
  public:
    BbAngleIt(const BbSolute &mt);
    ~BbAngleIt();
    void operator++();
    const Angle &operator()()const;
    operator bool()const;
  };

  class BbImpIt{
    BbImpIt_i *d_this;
    // not implemented
    BbImpIt();
    BbImpIt(const BbImpIt&);
    BbImpIt &operator=(const BbImpIt &);
  public:
    BbImpIt(const BbSolute &mt);
    ~BbImpIt();
    void operator++();
    const Improper &operator()()const;
    operator bool()const;
  };

  class BbDihIt{
    BbDihIt_i *d_this;
    // not implemented
    BbDihIt();
    BbDihIt(const BbDihIt&);
    BbDihIt &operator=(const BbDihIt &);
  public:
    BbDihIt(const BbSolute &mt);
    ~BbDihIt();
    void operator++();
    const Dihedral &operator()()const;
    operator bool()const;
  };
  
} /* Namespace */ 
#endif



