// gcore_System.h

#ifndef INCLUDED_GCORE_SYSTEM
#define INCLUDED_GCORE_SYSTEM

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gcore{

  class Solvent;
  class Molecule;
  class Box;
  /**
   * Class System
   * The system class in gromos++ contains everything for the 
   * that you were interested in. Coordinates and topological information 
   * of the Molecules and Solvent.
   *
   * Description:
   * The System class contains all information of the molecules you have
   * or will simulate. Coordinates (including box sizes) and topological 
   * information for your solute and solvent molecules are accessible 
   * from the system.
   *
   * @class System
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::Molecule
   * @sa gcore::Solvent
   * @sa gcore::Box
   */
  class System{
    std::vector<Molecule*> d_mol;
    std::vector<Solvent*> d_sol;
    Box *d_box;
    

  public:
    //Constructors
    /**
     * System Constructor
     */
    System();
    /**
     * System copy constructor
     */
    System(const System &sys);
    /**
     * System deconstructor
     */
    ~System();
    /**
     * Boolean to indicate whether a Box block has been read in.
     */
    bool hasBox;
    /**
     * Boolean to indicate whether a Velocity block has been read in.
     */
    bool hasVel;
    // Methods
    /**
     * Member operator = copies one System into the other
     */
    System &operator=(const System &sys);
    /**
     * Method to add a Molecule to your system
     *
     * After adding a molecule to your system, you can still modify
     * the coordinates, but the topological information is now considered
     * as const.
     * @param mol Molecule to be added
     */
    void addMolecule(const Molecule &mol);
    /**
     * Method to add a Solvent to your system
     *
     * After adding a solvent to your system, you can still modify the 
     * coordinates, but the topological information is now considered as
     * const.<br>
     * Note that even though gromos96 does not support this option, gromos++
     * can in principle handle systems with multiple solvents.
     * @param sol Solvent to be added
     */
    void addSolvent(const Solvent &sol);
    
    // Accessors
    /** 
     * Accessor, returns the i-th Molecule in the System as a const
     */
    const Molecule &mol(int i)const;
    /**
     * Accessor, returns the i-th Molecule in the System
     */
    Molecule &mol(int i);
    /**
     * Accessor, returns the i-th Solvent in the System as a const. Remember
     * that this is the i-th solvent type. All solvent molecules of this
     * type are stored within the Solvent class
     */
    const Solvent &sol(int i)const;
    /**
     * Accessor, returns the i-th Solvent in the System. Remember
     * that this is the i-th solvent type. All solvent molecules of this
     * type are stored within the Solvent class
     */
    Solvent &sol(int i);
    /**
     * Accessor, returns the Box dimensions of the System as a const
     */
    const Box &box()const;
    /** 
     * Accessor, returns the Box dimensions of the System
     */
    Box &box();
    /**
     * Accessor, returns the number of Molecules in the System
     */
    int numMolecules()const;
    /**
     * Accessor, returns the number of Solvents in the System. Note 
     * that this again the number of different Solvent types in the 
     * system. The number of solvent molecules (or rather total number
     * of solvent atoms) can be obtained from the Solvent class.
     */
    int numSolvents()const;
    
};

  inline const Molecule &System::mol(int i)const{
    assert (i < this->numMolecules());
    return *d_mol[i];
  }

  inline Molecule &System::mol(int i){
    assert (i < this->numMolecules());
    return *d_mol[i];
  }

  inline const Solvent &System::sol(int i)const{
    assert (i < this->numSolvents());
    return *d_sol[i];
  }

  inline Solvent &System::sol(int i){
    assert (i < this->numSolvents());
    return *d_sol[i];
  }

  inline const Box &System::box()const{
    return *d_box;
  }

  inline Box &System::box(){
    return *d_box;
  }

  inline int System::numMolecules()const{
    return d_mol.size();
  }

  inline int System::numSolvents()const{
      return d_sol.size();
  }
  
}
#endif
