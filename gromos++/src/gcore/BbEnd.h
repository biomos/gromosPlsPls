// gcore_BbEnd.h

#ifndef INCLUDED_GCORE_BBEND
#define INCLUDED_GCORE_BBEND

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class BbEnd_i;
  class AtomTopology;
  class Exclusion;
  class Bond;
  class Angle;
  class Dihedral;
  class Improper;
  class BeBondIt;
  class BeAngleIt;
  class BeImpIt;
  class BeDihIt;

  class BbEnd{
    BbEnd_i *d_this;
    // This class contains all topological information

    friend class BeBondIt;
    friend class BeAngleIt;
    friend class BeImpIt;
    friend class BeDihIt;
    
    // not implemented
    BbEnd &operator=(const BbEnd &);

  public:
    BbEnd();
    BbEnd(const BbEnd &);
    ~BbEnd();
    
    // Methods
    void addAtom(const AtomTopology &a);
    void addBond(const Bond &b);
    void addAngle(const Angle &b);
    void addDihedral(const Dihedral &b);
    void addImproper(const Improper &b);
    void setResName(const std::string &s);
    void setRep(const int i);
    
    // Accessors
    int numAtoms()const;
    const AtomTopology& atom(int i) const; 
    const std::string &resName()const;
    const int rep()const;
    
  }; /* class BbEnd */

  
  class BeBondIt_i;
  class BeAngleIt_i;
  class BeImpIt_i;
  class BeDihIt_i;

  class BeBondIt{
    BeBondIt_i *d_this;
    // not implemented
    BeBondIt();
    BeBondIt(const BeBondIt&);
    BeBondIt &operator=(const BeBondIt &);
  public:
    BeBondIt(const BbEnd &mt);
    ~BeBondIt();
    void operator++();
    const Bond &operator()()const;
    operator bool()const;
  };

  class BeAngleIt{
    BeAngleIt_i *d_this;
    // not implemented
    BeAngleIt();
    BeAngleIt(const BeAngleIt&);
    BeAngleIt &operator=(const BeAngleIt &);
  public:
    BeAngleIt(const BbEnd &mt);
    ~BeAngleIt();
    void operator++();
    const Angle &operator()()const;
    operator bool()const;
  };

  class BeImpIt{
    BeImpIt_i *d_this;
    // not implemented
    BeImpIt();
    BeImpIt(const BeImpIt&);
    BeImpIt &operator=(const BeImpIt &);
  public:
    BeImpIt(const BbEnd &mt);
    ~BeImpIt();
    void operator++();
    const Improper &operator()()const;
    operator bool()const;
  };

  class BeDihIt{
    BeDihIt_i *d_this;
    // not implemented
    BeDihIt();
    BeDihIt(const BeDihIt&);
    BeDihIt &operator=(const BeDihIt &);
  public:
    BeDihIt(const BbEnd &mt);
    ~BeDihIt();
    void operator++();
    const Dihedral &operator()()const;
    operator bool()const;
  };
  
} /* Namespace */ 
#endif



