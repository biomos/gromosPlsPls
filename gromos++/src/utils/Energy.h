// utils_Energy.h

// Energy class: calculates all energies of a system

#ifndef INCLUDED_UTILS_ENERGY
#define INCLUDED_UTILS_ENERGY

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_SET
#include <set>
#define INCLUDED_SET
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif

namespace gcore{
  class GromosForceField;
  class System;
}
namespace bound{
  class Boundary;
}
namespace gmath{
  class Vec;
}

namespace utils
{
  class AtomSpecifier;
  class PropertyContainer;
  class Property;
  class Energy;
  class Energy{
    gcore::System *d_sys;
    gcore::GromosForceField *d_gff;
    bound::Boundary *d_pbc;
    utils::AtomSpecifier *d_as;
    utils::PropertyContainer *d_pc;
    std::vector<double> d_cov, d_vdw_m, d_vdw_s, d_el_m, d_el_s;
    std::vector<gmath::Vec> d_covpar;
    utils::AtomSpecifier *d_soft;
    double d_al2, d_eps, d_kap, d_cut, d_p_vdw, d_p_el;
    std::vector<std::set<int> > d_ex, d_third;
  public: 
    // Constructor
    Energy(gcore::System &sys, gcore::GromosForceField &gff, 
           bound::Boundary &pbc);
    // Deconstructor
    ~Energy(){}
   
    // Methods: set the atoms to consider
    int setAtoms(utils::AtomSpecifier &as);

    // Methods: set the covalent properties to consider
    int setProperties(utils::PropertyContainer &pc);

    // Methods: set an atomspecifier with soft atoms, specify al2
    void setSoft(utils::AtomSpecifier &soft, double al2);

    // Methods: set properties for reaction field
    void setRF(double eps, double kap);

    // Methods: set cut-off
    void setCutOff(double cut);
    
    // Methods: calculate all energies
    void calc();

    // Methods: calculate pairwise interaction between to atoms,
    // numbering corresponds to the atomspecifier
    void calcPair(int i, int j, double &vdw, double &el);
    
    // Accessors: the total interaction energy
    double tot();

    // Accessors: the total non-bonded interaction energy
    double nb();

    // Accessors: the total VdW interaction energy
    double vdw();

    // Accessors: the total El interaction energy
    double el();

    // Accessors: the total covalent interaction energy
    double cov();

    // Accessors: the VdW energy of atom i (corresponding to atomspecifier)
    double vdw(int i);

    // Accessors: the El energy of atom i (corresponding to atomspecifier)
    double el(int i);

    // Accessors: the covalent energy of property i (corresponging to 
    // propertycontainer)
    double cov(int i);

    // Accessors: the VdW energy of atom i with solute
    double vdw_m(int i);

    // Accessors: the VdW energy of atom i with solvent
    double vdw_s(int i);

    // Accessors: the El energy of atom i with solute
    double el_m(int i);

    // Accessors: the El energy of atom i with solvent
    double el_s(int i);
    
    // Exception
    struct Exception: public gromos::Exception{
      Exception(const string &what): gromos::Exception("Energy", what){}
    };

    // define pi
    static const double pi=3.14159265359;

  protected:
    gmath::Vec calcChgrp(int i);
    int findBond(utils::Property &pp);
    int findAngle(utils::Property &pp);
    int findDihedral(utils::Property &pp);
    double calcBond(double val, gmath::Vec par);
    double calcAngle(double val, gmath::Vec par);
    double calcDihedral(double val, gmath::Vec par);
    double calcImproper(double val, gmath::Vec par); 
};
//inline functions and methods

inline double Energy::tot()
{
  return this->nb() + this->cov();
}
inline double Energy::nb()
{
  return this->vdw() + this->el();
}
inline void Energy::setCutOff(double cut)
{
  d_cut=cut;
}
inline void Energy::setSoft(utils::AtomSpecifier &soft, double al2)
{
  d_soft=&soft;
  d_al2=al2;
}
inline void Energy::setRF(double eps, double kap)
{
  d_eps=eps;
  d_kap=kap;
}
}
#endif
