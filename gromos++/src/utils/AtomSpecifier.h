// utils_AtomSpecifier.h

// Class that contains a sequential list of specific atoms

#ifndef INCLUDED_UTILS_ATOMSPECIFIER
#define INCLUDED_UTILS_ATOMSPECIFIER

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class System;
  class Molecule;
  class MoleculeTopology;
}

namespace utils
{
  class AtomSpecifier{
    std::vector<int> d_mol, d_atom;
    gcore::System *d_sys;
   
  public: 
    // Constructor
    AtomSpecifier(gcore::System &sys);
    // Constructor<br>
    // note that the string 's' is usually supplied by the user and is 
    // interpreted as GROMOS96 numbering (starting with 1)
    AtomSpecifier(gcore::System &sys, std::string s);
    // Deconstructor
    ~AtomSpecifier(){}
   
    // Methods: add specifier-string<br>
    // note that the string 's' is usually supplied by the user and is 
    // interpreted as GROMOS96 numbering (starting with 1)
    int addSpecifier(std::string s);
    // Methods: add atom 'a' in molecule 'm'<br>
    // note that this takes 'a' and 'm' in gromos++ numbering
    // (starting from 0)
    int addAtom(int m, int a);
    // Methods: remove atom 'a' in molecule 'm' from the specifier<br>
    // note that this takes 'a' and 'm' in gromos++ numbering 
    // (starting from 0)
    int removeAtom(int m, int a);
    // Methods: add all atoms of a certain type<br>
    int addType(std::string s);
    // Methods: copy AtomSpecifier 'as' to this
    AtomSpecifier &operator=(const AtomSpecifier &as);
    // Methods: add AtomSpecifier 'as' to this
    AtomSpecifier operator+(const AtomSpecifier &as);
    
    // Accessors: the molecule containing atom 'i' of the AtomSpecifier
    int mol(int i);
    // Accessors: a pointer to the molecule vector
    std::vector <int> *mol();
    // Accessors: the atomnumber of atom 'i' of the AtomSpecifier
    int atom(int i);
    // Accessors: a pointer to the atom vector
    std::vector <int> *atom();
    // The number of atoms in the AtomSpecifier
    int size();
    // A pointer to the system of which the atoms are specified
    gcore::System *sys();
    
    
    // Exception
    struct Exception: public gromos::Exception{
      Exception(const string &what): gromos::Exception("AtomSpecifier", what){}
    };
    
};
  //inline functions and methods

  inline int AtomSpecifier::mol(int i)
    {
      return d_mol[i];
    }
  inline std::vector<int> *AtomSpecifier::mol()
    {
      return &d_mol;
    }
  inline int AtomSpecifier::atom(int i)
    {
      return d_atom[i];
    }
  inline std::vector<int> *AtomSpecifier::atom()
    {
      return &d_atom;
    }
  inline int AtomSpecifier::size()
    {
      return d_atom.size();
    }
  
  inline gcore::System *AtomSpecifier::sys()
    {
      return d_sys;
    }
  
}
#endif
