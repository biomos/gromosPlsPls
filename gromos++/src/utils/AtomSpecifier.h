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

namespace gmath
{
  class Vec;
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
    std::vector<int> d_mol, d_atom, d_solventType;
    gcore::System *d_sys;
    int d_nsm;
    
  public: 
    // Constructors
    /** 
     * AtomSpecifier standard constructor
     */
    AtomSpecifier(){d_nsm=-1;};
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
     * Method to set the system the atom specifier is referring to
     * @param sys the system
     */
    void setSystem(gcore::System &sys);
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
     * Method to add a single molecule to the AtomSpecifier without redundancy checks
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule.
     */
    int addAtomStrict(int m, int a);
    /**
     * Method to remove an atom from the AtomSpecifier.
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule
     */
    int removeAtom(int m, int a);
    /**
     * Method to remove an atom from the AtomSpecifier.
     *
     * @param int i remove the atoms with index 1 in the specifier
     */
    int removeAtom(int i);
    /**
     * Method to find the index of a specific atom in the AtomSpecifier
     *
     * Numbering is assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule
     */
    int findAtom(int m, int a);
    /**
     * Method to add all atoms of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     */
    int addType(int m, std::string s);
    /**
     * Method to add all atoms (in all molecules) with a certain name to the 
     * AtomSpecifier
     * @param s Atom name that is to be added (e.g. CA)
     */ 
    int addType(std::string s);
    /**
     * Method to add solvent atoms of a type
     * @param s Atom name that is to be added
     */
    int addSolventType(std::string s);
    /**
     * Method to sort the atoms ascending order. Some applications might
     * need the atoms to be ordered. This is just a simple bubble sort
     */
    void sort();
    /**
     * Member operator = copies one AtomSpecifier into the other
     */
    AtomSpecifier &operator=(const AtomSpecifier &as);
    /**
     * Member operator + adds two AtomSpecifiers. Atoms that appeared in
     * both are only listed once.
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
     * Function to empty the AtomSpecifier
     */
    void clear();
    /**
     * Accessor, returns a pointer to the system on which the AtomSpecifier
     * is based
     */
    gcore::System *sys();
    /**
     * Accessor, returns a pointer to the coordinates of the i-th 
     * atom in the AtomSpecifier   
     */
    gmath::Vec *coord(int i);
    /**
     * Accesor, returns the atom name of the i-th atom in the AtomSpecifier
     */
    std::string name(int i);
    /**
     * Accessor, returns the Iac of the i-th atom in the AtomSpecifier
     */
    int iac(int i);
    /**
     * Accessor, returns the charge of the i-th atom in the AtomSpecifier
     */
    double charge(int i);
    
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
    /**
     * Interprets atoms from the AtomSpecifier
     */
    void _parseAtomsHelper(std::string substring, int &mol);
    /**
     * Adds an atom to the atom specifier, checks whether it is in already
     */
    void _appendAtom(int m, int a);
    /**
     * Compares the string s to the name of the a-th atom in the m-th 
     * molecule of the system. If they are the same, it returns true. Only
     * the characters before a possible '?' in s are compared.
     */
    bool _checkName(int m, int a, std::string s);
    /**
     * Method that expands specified solvent types to the real number of
     * solvent molecules. Is called from any accessor number of solvent 
     * molecules has changed  
     */
    int _expandSolvent();
    /**
     * Tells you if the number of solvent molecules in the system has changed
     * if so, the accessors should re-expand the Solvent Types.
     */
    bool _expand();
    /**
     * special function for sorting, to replace 'larger than' in comparisons
     * makes sure that solvent atoms will be last
     * @param int i index of that is compared to
     * @param int m molecule number
     * @param int a atom number
     * @return bool returns true if atom i should come after m:a
     *                      false if atom i should come before m:a
     */
    bool _compare(int i, int m, int a);
    
};
  //inline functions and methods


  inline std::vector<int> *AtomSpecifier::mol()
    {
      return &d_mol;
    }

  inline std::vector<int> *AtomSpecifier::atom()
    {
      return &d_atom;
    }
  
  inline gcore::System *AtomSpecifier::sys()
    {
      return d_sys;
    }


}
#endif
