// gcore_SolventTopology.h

#ifndef INCLUDED_GCORE_SOLVENTTOPOLOGY
#define INCLUDED_GCORE_SOLVENTTOPOLOGY

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{

  class SolventTopology_i;
  class GromosForceField;
  class AtomTopology;
  class Constraint;
  class ConstraintIterator;

  class SolventTopology{
    SolventTopology_i *d_this;
    // This class contains all topological information

    friend class ConstraintIterator;
    
    // not implemented
    SolventTopology &operator=(const SolventTopology &);

  public:
    SolventTopology();
    SolventTopology(const SolventTopology &);
    ~SolventTopology();
    
    // Methods
    void addAtom(const AtomTopology &a);
    void addConstraint(const Constraint &b);
    void setSolvName(const std::string &s);
    
    // Accessors
    int numAtoms()const;
    const AtomTopology& atom(int i) const; 

    const std::string &solvName()const;
    
  }; /* class MoleculeTopology */


  class ConstraintIterator_i;

  class ConstraintIterator{
    ConstraintIterator_i *d_this;
    // not implemented
    ConstraintIterator();
    ConstraintIterator(const ConstraintIterator&);
    ConstraintIterator &operator=(const ConstraintIterator &);
  public:
    ConstraintIterator(const SolventTopology &mt);
    ~ConstraintIterator();
    void operator++();
    const Constraint &operator()()const;
    operator bool()const;
  };

} /* Namespace */ 
#endif



