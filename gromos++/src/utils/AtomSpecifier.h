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
  /**
   * Class AtomSpecifier
   * purpose: contains specific atoms of the system, keeping track of 
   * molecule and atom numbers
   *
   * Description:
   * The AtomSpecifier can be used to look over a specific set of atoms,
   * possibly spanning different molecules. A 'specifier' is a string with 
   * following format <mol>[:<atom>[-<atom>]]. For example "1:3-5" means atoms
   * 3,4 and 5 of molecule 1; "2:5" means atom 5 of molecule 2; "3" means all
   * atoms of molecule 3.
   *
   * @class AtomSpecifier
   * @author C. Oostenbrink
   * @ingroup utils
   * @sa utils::PropertySpecifier
   */
  class AtomSpecifier{
    std::vector<int> d_mol, d_atom;
    gcore::System *d_sys;
   
  public: 
    // Constructor
    /**
     * AtomSpecifier Constructor
     * @param sys The AtomSpecifier needs to know about the system. It 
     *            does not know about any atoms yet.
     */
    AtomSpecifier(gcore::System &sys);
    /**
     * AtomSpecifier Constructor
     * @param sys The AtomSpecifier needs to know about the system.
     * @param s   A string of the correct format. Usually this is provided
     *            by the user, so it is assumed to start numbering at 1
     */
    AtomSpecifier(gcore::System &sys, std::string s);
    /**
     * AtomSpecifier Deconstructor
     */
    ~AtomSpecifier(){}
   
    /**
     * Method to add parse a string to the AtomSpecifier.
     * @param s Is assumed to be user-specified, with numbering starting at 1
     */
    int addSpecifier(std::string s);
    /**
     * Method to add a single molecule to the AtomSpecifier
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule.
     */
    int addAtom(int m, int a);
    /**
     * Method to remove an atom from the AtomSpecifier.
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule
     */
    int removeAtom(int m, int a);
    /**
     * Method to add all atoms with a certain name to the AtomSpecifier
     * @param s Atom name that is to be added (e.g. CA)
     */    
    int addType(std::string s);
    /**
     * Member operator = copies one AtomSpecifier into the other
     */
    AtomSpecifier &operator=(const AtomSpecifier &as);
    /**
     * Member operator + adds two AtomSpecifiers. Just appends the two,
     * there are no checks if atoms might be listed twice now.
     */
    AtomSpecifier operator+(const AtomSpecifier &as);

    /**
     * Accessor, returns the molecule number of the i-th atom in the
     * AtomSpecifier
     */    
    int mol(int i);
    /**
     * Accessor, returns a pointer to the vector containing the molecule
     * numbers.
     */
    std::vector <int> *mol();
    /**
     * Accessor, returns the atom number of the i-th atom in the 
     * AtomSpecifier
     */
    int atom(int i);
    /**
     * Accessor, returns a pointer to the vector containing the atom
     * numbers.
     */
    std::vector <int> *atom();
    /**
     * Accessor, returns the number of atoms in the AtomSpecifier
     */
    int size();
    /**
     * Accessor, returns a pointer to the system on which the AtomSpecifier
     * is based
     */
    gcore::System *sys();
    
    
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says AtomSpecifier, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const string &what): gromos::Exception("AtomSpecifier", what){}
    };
  protected:
    //Internal function
    /**
     * Parse the arguments string into the AtomSpecifier
     */
    void parse(std::string s);
    void _parseAtomsHelper(std::string substring, int &mol);
    
    
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
