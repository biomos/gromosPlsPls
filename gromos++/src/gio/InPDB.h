namespace gio{
  
  class InPDB_i;
  /**
   * Class InPDB handles the reading of pdb files for further us of the ATOM or
   * HETATOM informations (name, resname, coordinates, ...).
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
    InPDB(const std::string &filename, bool readATOM = true, bool readHETATOM = false);
    /**
     * Destructor
     */
    ~InPDB();

    // Methods
    /**
     * select atoms to be read (ATOM, HETATOM, ...)
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
     * Accessor to return the residue number
     */
    unsigned int getResNumber(unsigned int i);
    /**
     * Accessor to return the atom name
     */
    std::string getAtomName(unsigned int i);
    /**
     * Accessor to return the chain name
     */
    std::string getChain(unsigned int i);
    /**
     * Exception to throw error messages
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("InPDB", what_arg){}
    };
  };
}

//#endif