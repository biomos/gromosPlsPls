// fit_Reference.h

#ifndef INCLUDED_FIT_WEIGHTS
#define INCLUDED_FIT_WEIGHTS

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

namespace gcore{
  class System;
}

namespace fit{

  class Reference_i;
  
  class Reference{
    Reference_i *d_this;

    // not implemented
    Reference();
    Reference(const Reference &);
    Reference &operator=(const Reference &);

  public:
    Reference(gcore::System *sys);
    ~Reference();
    
    // methods
    void addClass(int mol, const std::string &name);
    /*
       adds a class of atoms with names name
       of molecule mol
       Attention: equal weights will be given to all
       reference atoms!
    */
    void addAtom(int m, int i);
    /*    
      adds atom i of molecule m 
      and gives all reference atoms equal
      weights
    */
    void rescale();
      // gives equal weight to all non-zero elements

    void setWeight(int m, int i, double w);
    // sets weight w to atom i of molecule m

    void normalise();
      // normalises weights to 1

   void makePosList (const gcore::System &,int molecule, const std::string &atomtype, std::vector<int> &poslist);

    // accessor
    gcore::System &sys();
    const gcore::System &sys()const;
    // return reference to system.
    double weight(int m, int i)const;
    // return weight of atom i of molecule m
  };
}

#endif
