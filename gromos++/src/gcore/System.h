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

  class System{
    std::vector<Molecule*> d_mol;
    std::vector<Solvent*> d_sol;
    Box *d_box;

  public:
    //Constructors
    System();
    System(const System &sys);
    ~System();

    // Methods
    System &operator=(const System &sys);
    void addMolecule(const Molecule &mol);
    void addSolvent(const Solvent &sol);
    
    // Accessors
    const Molecule &mol(int i)const;
    Molecule &mol(int i);
    const Solvent &sol(int i)const;
    Solvent &sol(int i);
    const Box &box()const;
    Box &box();
    int numMolecules()const;
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
