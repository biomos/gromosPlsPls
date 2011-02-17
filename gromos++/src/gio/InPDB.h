// gio_InPDB.h

//#ifndef INCLUDED_GIO_INPDB
//#define INCLUCED_GIO_INPDB

//#ifndef INCLUDED_STRING
//#include <string>
//#define INCLUDED_STRING
//endif

//#ifndef INCLUDED_GROMOS_EXCEPTION
//#include "../gromos/Exception.h"
//#endif

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
     * Exception
     */
    struct Exception: public gromos::Exception{
      Exception(const std::string& what_arg) : gromos::Exception("InPDB", what_arg){}
    };
  };
}

//#endif
