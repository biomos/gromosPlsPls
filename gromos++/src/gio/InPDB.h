namespace gio{
  
  class InPDB_i;
  /**
   * Class InPDB
   * Some description...
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
    /**
     * Constructor
     */
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
     * Modifier: renames the residue i in sequence
     */
    void changeResSeq(unsigned int i, std::string newname);
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
     * Accessor to return the atom name
     */
    std::string getAtomName(unsigned int i);
    /**
     * Exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("InPDB", what_arg){}
    };
  };
}

//#endif
