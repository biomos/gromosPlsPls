/**
 * @file AtomSpecifier.h
 * AtomSpecifier methods
 */

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

// minimal complete headers
#include "../gcore/Molecule.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/System.h"

#include "../utils/VirtualAtom.h"

namespace gmath
{
  class Vec;
}

namespace utils
{
  enum spec_type{ spec_solute, spec_solvent, spec_virtual };
  
  class AtomSpecifier;

  /**
   * Class SpecAtom
   * purpose: interface to access a specific atom (in an AtomSpecifier)
   * unified access to solute, solvent and virtual atoms.
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class SpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class SpecAtom
  {
  public:
    SpecAtom(gcore::System &sys, int m, int a) : d_sys(&sys), d_mol(m), d_atom(a) {}
    SpecAtom(SpecAtom const & s) : d_sys(s.d_sys), d_mol(s.d_mol), d_atom(s.d_atom) {}
    
    virtual ~SpecAtom() {};

    virtual spec_type type() const { return spec_solute; }
    
    virtual SpecAtom * clone() { return new SpecAtom(*this); }

    virtual int mol() const { return d_mol; }
    virtual int atom()const { return d_atom; }
    
    virtual gmath::Vec & pos() { return d_sys->mol(d_mol).pos(d_atom); }
    virtual gmath::Vec const & pos() const { return d_sys->mol(d_mol).pos(d_atom); }
    
    virtual std::string name()const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).name();
    }
    
    virtual int iac() const 
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).iac();
    }
    virtual double charge() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).charge();
    }
    virtual double mass() const
    {
      return d_sys->mol(d_mol).topology().atom(d_atom).mass();
    }
    virtual int resnum() const
    {
      return d_sys->mol(d_mol).topology().resNum(d_atom);
    }
    virtual std::string resname() const
    {
      return d_sys->mol(d_mol).topology().resName(resnum());
    }
    
    virtual void setSystem(gcore::System &sys)
    {
      d_sys = &sys;
    }
    
    virtual std::string toString() const;

    virtual utils::AtomSpecifier conf();
    
    virtual int virtualType() { return 0; }
    
    
  protected:
    gcore::System *d_sys;
    int d_mol;
    int d_atom;

    struct Exception: public gromos::Exception{
      Exception(const std::string &what): 
	gromos::Exception("AtomSpecifier", what){}
    };

  };

  /**
   * Class SolventSpecAtom
   * purpose: interface to access a specific solvent atom (in an AtomSpecifier)
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class SolventSpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class SolventSpecAtom : public SpecAtom
  {
  public:
    SolventSpecAtom(gcore::System &sys, int m, int a) : SpecAtom(sys, m, a) {}
    SolventSpecAtom(SolventSpecAtom const &s) : SpecAtom(s) {}
    
    virtual ~SolventSpecAtom() {};
    
    virtual spec_type type() const { return spec_solvent; }
    
    virtual SpecAtom * clone() { return new SolventSpecAtom(*this); }
    
    virtual int mol() const { return d_mol; }
    virtual int atom()const { return d_atom; }
    
    virtual gmath::Vec & pos() {
      if(d_atom > d_sys->sol(0).numPos())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).pos(d_atom); 
    }
    virtual gmath::Vec const & pos() const { 
      if(d_atom > d_sys->sol(0).numPos())
	throw Exception(" solvent coordinate not read");
      return d_sys->sol(0).pos(d_atom); 
    }
    
    virtual std::string name()const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).name();
    }
    
    virtual int iac() const 
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).iac();
    }
    virtual double charge() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).charge();
    }
    virtual double mass() const
    {
      int num=(d_atom % d_sys->sol(0).topology().numAtoms());
      return d_sys->sol(0).topology().atom(num).mass();
    }
    virtual int resnum() const
    {
      return d_atom / d_sys->sol(0).topology().numAtoms();
    }
    virtual std::string resname() const
    {
      return "SLV";
    }
    
  private:

  };

  /**
   * Class VirtualSpecAtom
   * purpose: interface to access a specific atom (in an AtomSpecifier)
   * unified access to solute, solvent and virtual atoms.
   *
   * Description:
   * An AtomSpecifier is a vector of SpecAtom *, through which all relevant
   * information about the atom is accessible.
   *
   * @class VirtualSpecAtom
   * @author M. Christen
   * @ingroup utils
   * @sa utils:AtomSpecifier
   */
  class VirtualSpecAtom : public SpecAtom
  {
  public:
    VirtualSpecAtom(gcore::System &sys, std::string s, VirtualAtom::virtual_type t);
    VirtualSpecAtom(VirtualSpecAtom const &s) : SpecAtom(s), d_va(s.d_va), d_pos(s.d_pos) {}
    
    virtual ~VirtualSpecAtom() {};

    virtual spec_type type() const { return spec_virtual; }

    virtual SpecAtom * clone() { return new VirtualSpecAtom(*this); }
    
    virtual int mol() const { /* throw Exception(" accessing VA mol");*/ return d_mol; }
    virtual int atom()const { /* throw Exception(" accessing VA atom"); */ return d_atom; }
    
    virtual gmath::Vec & pos() {return d_pos = d_va.pos(); }
    virtual gmath::Vec const & pos() const { return d_pos = d_va.pos(); }
    
    virtual std::string name()const
    {
      return "VA";
    }
    
    virtual int iac() const 
    {
      throw Exception(" accessing VA iac");
      // return 0;
    }
    virtual double charge() const
    {
      throw Exception(" accessing VA charge");
      // return 0;
    }
    virtual double mass() const
    {
      throw Exception(" accessing VA mass");
      // return 0.0;
    }
    virtual int resnum() const
    {
      throw Exception(" accessing VA resnum");
      // return 0;
    }
    virtual std::string resname() const
    {
      throw Exception(" accessing VA resname");
      // return "VA";
    }
    
    virtual void setSystem(gcore::System &sys);

    virtual std::string toString() const
    {
      return d_va.toString();
    }
    virtual utils::AtomSpecifier conf();

    virtual int virtualType()
    {
      return d_va.type();
    }
    
  protected:
    VirtualAtom d_va;
    mutable gmath::Vec d_pos;
  };

  /**
   * @class AtomSpecifier
   * @author C. Oostenbrink, M. Christen
   * @ingroup utils
   *
   * Class AtomSpecifier
   * purpose: contains specific atoms of the system, keeping track of 
   * molecule and atom numbers
   *
   * Description:
   * This class defines (and implements) a general form to
   * access atoms in a system.
   *
   * @section AtomSpecifier Atom Specifier
   * The AtomSpecifier can be used to look over a specific set of atoms,
   * possibly spanning different molecules. A 'specifier' is a string with 
   * following format:
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <mol>[-<mol>]:<atom>[-<atom>] @endverbatim
   * </b></span>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <mol>[-<mol>]:res(<resnr>:<atom>[-<atom>]) @endverbatim
   * </b></span>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim va(<type>, <atomspec>) @endverbatim
   * </b></span>
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim <atomspec>[;<atomspec>] @endverbatim
   * </b></span>
   * <br>
   * where
   * - <mol> is a molecule number, or 'a' for all molecules
   * - <atom> is an atom number, or 'a' for all atoms, or an atom name
   * - <resnr> is a residue number or a residue name
   * - <atomspec> is a Atom Specifier
   * - res() starts a residue specification
   * - va() starts a virtual atom
   *
   * For example "1:3-5" means atoms
   * 3,4 and 5 of molecule 1; "2:5" means atom 5 of molecule 2; 
   * "3:a" means all atoms of molecule 3.
   *
   * An atom can also be a <b>virtual atom</b>.
   * <span style="color:darkred;font-size:larger"><b>
   * @verbatim va(<type>,<AtomSpecifier>) @endverbatim
   * </b></span>
   * <br>
   * with type:
   * - Gromos96 virtual atom types
   * - com (centre of mass)
   * - cog (centre of geometry)
   *
   * <b>See also</b> @ref PropertySpecifier "Property Specifier"
   *
   *
   */
  class AtomSpecifier{
    std::vector<int> d_solventType;
    
    std::vector<SpecAtom *> d_specatom;
    
    gcore::System *d_sys;
    int d_nsm;
    
  public: 
    // Constructors
    /** 
     * AtomSpecifier standard constructor
     */
    AtomSpecifier(){d_sys=NULL; d_nsm=-1;};
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
     * copy constructor!
     */
    AtomSpecifier(AtomSpecifier const & as);
    /**
     * AtomSpecifier Destructor
     */
    ~AtomSpecifier();

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
     * Method to add a single molecule to the AtomSpecifier 
     * without redundancy checks
     *
     * Numbering is here assumed to be gromos++ numbering, starting at 0
     * @param m Number of the molecule the atom belongs to
     * @param a Atom number within that molecule.
     */
    int addAtomStrict(int m, int a);
    /**
     * Method to add an atom from a Gromos96 number.
     * Numbering is here assumed to be Gromos96 (starting at 1)
     * @param a Gromos96 atom number
     */
    int addGromosAtom(int a);
    /**
     * Method to add a complete molecule.
     * Numbering is here assumed to be gromos++ (starting at 0)
     * @param m Molecule number.
     */
    int addMolecule(int m);
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
     * Method to add atoms of the specified range
     * of the specified molecule with a certain name
     * to the AtomSpecifier
     * @param m number of the molecule to consider
     * @param s Atom name that is to be added (e.g. CA)
     * @param beg begin of range
     * @param end end of range
     */
    int addType(int m, std::string s, int beg, int end);
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
     * Accessor, returns the atom number of the i-th atom in the 
     * AtomSpecifier
     */
    int atom(int i);
    /**
     * Accessor, returns the gromos atom number of the i-th atom in the
     * AtomSpecifier
     */
    int gromosAtom(int i);
    /**
     * Accessor, returns a pointer to the vector containing the atom
     * numbers.
     */
    std::vector<SpecAtom *> & atom();
    /**
     * const Accessor, returns a pointer to the vector containing the atom
     * numbers.
     */
    std::vector<SpecAtom *> const & atom()const;
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
     * Accessor, returns the coordinates of the i-th
     * atom in the AtomSpecifier   
     */
    gmath::Vec & pos(int i);
    /**
     * const Accessor, returns the coordinates of the i-th
     * atom in the AtomSpecifier   
     */
    // gmath::Vec const & pos(int i)const;
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
     * Accessor, returns the mass of the i-th atom in the AtomSpecifier
     */
    double mass(int i);
    /**
     * Accessor, returns the residue number of atom i in the AtomSpecifier
     */
    int resnum(int i);
    /**
     * Accessor, returns the residue name of the i-th atom
     * in the AtomSpecifier
     */
    std::string resname(int i);
    /**
     * Method, returns a vector of strings that would reproduce the
     * AtomSpecifier
     */
    std::vector<std::string> toString();
    /**
     * Method, returns a string of just the i-th atom in the 
     * AtomSpecifier
     */
    std::string toString(int i);
    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says AtomSpecifier, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): 
	gromos::Exception("AtomSpecifier", what){}
    };
  protected:
    //Internal function
    /**
     * Parse the arguments string into the AtomSpecifier
     */
    void parse(std::string s);
    /**
     * parse a single specifier (no ';')
     */
    void parse_single(std::string s);
    /**
     * parse the molecules, deliver in mol
     */
    void parse_molecule(std::string s, std::vector<int> & mol);
    /**
     * parse a range of molecules
     */
    void parse_mol_range(std::string s, std::vector<int> & mol);
    /**
     * parse the atom part
     * can also include residues
     */
    void parse_atom(int mol, std::string s);
    /**
     * parse an atom from a range of possible atoms
     * the range is either all atoms of the molecule mol
     * or just the atoms inside a residue of the molecule mol
     */
    void parse_atom_range(int mol, int beg, int end, std::string s);
    /**
     * parse virtual atoms.
     * standard virtual types plus
     * centre of geometry (cog) and centre of mass (com)
     */
    void parse_va(std::string s);
    /**
     * parse residues
     */
    void parse_res(int mol, std::string res, std::string atom);
    /**
     * parse residues if specified by type
     */
    void parse_res_type(int mol, std::string res, std::string s);
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
     * Compares the string s to the name of the residue 
     * the  a-th atom in the m-th molecule of the system belong to.
     * If they are the same, it returns true. Only
     * the characters before a possible '?' in s are compared.
     */
    bool _checkResName(int m, int a, std::string s);
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
    /**
     * find matching closing bracket for a given opening bracket
     * at position it in the string s
     */
    std::string::size_type find_matching_bracket(std::string s, 
						 char bra='(',
						 std::string::size_type it=0);
    /**
     * find the atoms belonging to residue res of molecule mol
     */
    void res_range(int mol, int res, int &beg, int &end);
    
};
  //inline functions and methods

  inline std::vector<SpecAtom *> & AtomSpecifier::atom()
  {
    return d_specatom;
  }

  inline std::vector<SpecAtom *> const & AtomSpecifier::atom()const
  {
    return d_specatom;
  }

  inline gcore::System *AtomSpecifier::sys()
  {
    return d_sys;
  }

  inline VirtualSpecAtom::VirtualSpecAtom(gcore::System &sys, 
					  std::string s,
					  VirtualAtom::virtual_type t)
    : SpecAtom(sys, -3, -1),
      d_va(sys, AtomSpecifier(sys, s), t)
  {}

  inline void VirtualSpecAtom::setSystem(gcore::System &sys)
  {
    d_va.setSystem(sys);
    SpecAtom::setSystem(sys);
  }

  inline AtomSpecifier VirtualSpecAtom::conf()
  {
    return d_va.conf();
  }
  inline AtomSpecifier SpecAtom::conf()
  {
    AtomSpecifier as;
    return as;
  }
}



#endif

