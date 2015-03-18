#ifndef INCLUDED_INPDB
#define INCLUDED_INPDB

namespace gio{
  
  class InPDB_i;
  /**
   * Class InPDB handles the reading of pdb files for further us of the ATOM or
   * HETATM informations (name, resname, coordinates, ...).
   *
   * @class InPDB
   * @ingroup gio
   * @author A. Eichenberger
   */
  class InPDB{

  private:
    /**
     * internal data object pointer
     */
    InPDB_i *d_this;
   
  public:
    InPDB(const std::string &filename, bool readATOM = true, bool readHETATM = false);
    /**
     * Destructor
     */
    ~InPDB();

    // Methods
    /**
     * select atoms to be read (ATOM, HETATM, ...)
     */
    void select(const std::string &thing);
    /**
     * read the PDB file
     */
    void read();
    /**
     * Accessor: returns the residue sequence read from PDB
     */
    std::vector<std::string> getResSeq();
    /**
     * renumber the residues sequentially
     */
    void renumberRes();
    /**
     * Modifier: renames the residue i in sequence
     */
    //void changeResSeq(unsigned int i, std::string newname);
    /**
     * Accessor: returns the position of PDB atom i
     */
    gmath::Vec getAtomPos(unsigned int i);
    /**
     * Accessor to return the number of atoms
     */
    unsigned int numAtoms(void);
    /**
     * Accessor to return the residue name
     */
    std::string getResName(unsigned int i);
    /**
     * Sets the residue name of atom i to s
     * @param i atom number
     * @param s residue name to be set
     */
    void setResName(int i, std::string s);
    /**
     * Accessor to return the residue number
     */
    unsigned int getResNumber(unsigned int i);
    /**
     * Accessor to return the atom name
     */
    std::string getAtomName(unsigned int i);
    /**
     * Accessor to return the type (ATOM or HETATM)
     */
    std::string getType(unsigned int i);
    /**
     * Accessor to return sequential atom number
     */
    unsigned int getSerial(unsigned int i);
    /**
     * Accessor to occupancy
     */
    double getOcc(unsigned int i);
    /**
     * Accessor to B-factor
     */
    double getB(unsigned int i);
    /**
     * Accessor to return the chain name
     */
    char getChain(unsigned int i);
    /**
     * Accessor to return the insertion code
     */
    char getICode(unsigned int i);
    /**
     * Accessor to return the element name
     */
    std::string getElement(unsigned int i);
    /**
     * Accessor to return the atom charge
     */
    std::string getCharge(unsigned int i);
    /*
     * Indicates whether the current InPDB has information about the atom serial numbers
     */
    bool hasSerials();
    /*
     * Indicates whether the current InPDB has information about the atom names
     */
    bool hasAtomNames();
    /*
     * Indicates whether the current InPDB has information about the residue names
     */
    bool hasResNames();
    /*
     * Indicates whether the current InPDB has information about the chain ID
     */
    bool hasChainIDs();
    /*
     * Indicates whether the current InPDB has information about the residue sequence number
     */
    bool hasResSeqs();
    /*
     * Indicates whether the current InPDB has information about the insertion code
     */
    bool hasiCodes();
    /*
     * Indicates whether the current InPDB has information about the x-coodinate
     */
    bool hasX();
    /*
     * Indicates whether the current InPDB has information about the y-coodinate
     */
    bool hasY();
    /*
     * Indicates whether the current InPDB has information about the Z-coodinate
     */
    bool hasZ();
    /*
     * Indicates whether the current InPDB has information about the occupancy
     */
    bool hasOccupancies();
    /*
     * Indicates whether the current InPDB has information about the temp. factor
     */
    bool hasTempFactors();
    /*
     * Indicates whether the current InPDB has information about the element names
     */
    bool hasElements();
    /*
     * Indicates whether the current InPDB has information about the charges
     */
    bool hasCharges();
    /**
     * Exception to throw error messages
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("InPDB", what_arg){}
    };
  };
}

#endif
