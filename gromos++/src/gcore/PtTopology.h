#ifndef INCLUDED_GCORE_PTTOPOLOGY
#define INCLUDED_GCORE_PTTOPOLOGY

namespace gcore
{
  /**
   * Class PtTopology
   * contains one or several perturbation topologies for 
   * nonbonded interactions
   *
   * @class PtTopology
   * @author C. Oostenbrink, N. Schmid
   * @ingroup gcore
   */
  class PtTopology
  {
    std::vector<int> d_atomnum;
    std::vector<std::string> d_atomnames;
    std::vector<std::string> d_pertnames;
    std::vector< std::vector <int> > d_iac;
    std::vector< std::vector <double> > d_charge;
    std::vector< std::vector <double> > d_mass;
    std::vector< std::vector <double> > d_polarisability;
    std::vector< std::vector <double> > d_dampingLevel;
    std::vector<double> d_alphaLJ;
    std::vector<double> d_alphaCRF;
    bool d_hasPolaristaionParams;
  
  public:
    /**
     * constructor
     */
    PtTopology() : d_hasPolaristaionParams(false) {};
    /**
     * deconstructor
     */
    ~PtTopology(){};
    /**
     * set the dimensions of the perturbation topology
     * @param a number of atoms that are to be perturbed
     * @param p number of perturbation topologies that 
     *          should be stored
     */
    void setSize(int a, int p);
    /**
     * function to set the iac for atom a in perturbation p
     */
    void setIac(int a, int p, int iac);
    /**
     * function to set the mass for atom a in perturbation p
     */
    void setMass(int a, int p, double m);
    /**
     * function to set the charge for atom a in perturbation p
     */
    void setCharge(int a, int p, double q);
    /**
     * function to set the polarisability for atom a in perturbation p
     */
    void setPolarisability(int a, int p, double al);
    /**
     * function to set the damping level for atom a in perturbation p
     */
    void setDampingLevel(int a, int p, double l);
    /**
     * function to set the atom name of atom a
     */
    void setAtomName(int a, std::string name);
    /**
     * function to set the name of perturbation p
     */
    void setPertName(int p, std::string name);
    /**
     * function to set the atom number of atom a
     */
    void setAtomNum(int a, int num);
    /**
     * function to set the softness for LJ for of atom a
     */
    void setAlphaLJ(int a, double al);
    /**
     * function to set the softness for CRF for of atom a
     */
    void setAlphaCRF(int a, double al);
    /**
     * function to set whether the topology contained polarisation
     * parameters
     */
    void setHasPolarisationParameters(bool pol);
    /**
     * function to apply a given perturbation to the system
     */
    void apply(gcore::System &sys, int iipt=0);
    /**
     * function to return the molecule and atom number
     * from the atom number in the perturbation topology
     */
    void findAtom(gcore::System &sys, int &m, int &a, int counter);
    
    /**
     * accessor to a vector of the charges in perturbation p
     */
    std::vector<double> charge(int p=0){return d_charge[p];}
    /**
     * accessor to the charge of atom a in perturbation p
     */
    double charge(int a, int p){return d_charge[p][a];}
    /**
     * accessor to a vector of the iac in perturbation p
     */
    std::vector<int> iac(int p=0){return d_iac[p];}
    /**
     * accessor to the iac of atom a in perturbation p
     */
    int iac(int a, int p){return d_iac[p][a];}
    /**
     * accessor to a vector of the masses in perturbation p
     */
    std::vector<double> mass(int p=0){return d_charge[p];}
    /**
     * accessor to the mass of atom a in perturbation p
     */
    double mass(int a, int p){return d_mass[p][a];}
    /**
     * accessor to a vector of the polaisabilities in perturbation p
     */
    std::vector<double> polarisability(int p=0){return d_polarisability[p];}
    /**
     * accessor to the polarisability of atom a in perturbation p
     */
    double polarisability(int a, int p){return d_polarisability[p][a];}
    /**
     * accessor to a vector of the damping levels in perturbation p
     */
    std::vector<double> dampingLevel(int p=0){return d_dampingLevel[p];}
    /**
     * accessor to the damping level of atom a in perturbation p
     */
    double dampingLevel(int a, int p){return d_dampingLevel[p][a];}
    
    /**
     * accessor to the name of atom a
     */
    std::string atomName(int a){return d_atomnames[a];}
    /**
     * accessor to the name of perturbation p
     */
    std::string pertName(int p){return d_pertnames[p];}
    /**
     * accessor to the atom number of atom a
     */
    int atomNum(int a){return d_atomnum[a];}
    /**
     * accessor to the softness in LJ of atom a
     */
    double alphaLJ(int a){return d_alphaLJ[a];}
    /**
     * accessor to the softness in CRF of atom a
     */
    double alphaCRF(int a){return d_alphaCRF[a];}
    /**
     * accessor to the number of perturbations in the topology
     */
    int numPt(){return d_iac.size();}
    /**
     * accessor to the number of atoms in the perturbation topology
     */
    int numAtoms(){return d_atomnames.size();}
    /**
     * accessor to polarisation
     */
    bool hasPolarisationParameters(){return d_hasPolaristaionParams;}
  
  };
}

#endif
