namespace

#include "AtomSpecifier.h"
 utils {

  // the implementation class, just to let the compiler know that it exists
  class iRDF;

  /**
   * Class RDF:
   * a class to calculate radial distribution functions of a given system
   */
  class RDF {

  private:
    
    // THE CLASSES DATA (or a pointer to the implementation class)
    // ================
    //
    /**
     * pointer to the data of the PDB class defined by the implementation
     * class iPDB
     */
    iRDF *d_this;

  public:

    // CONSTRUCTORS
    // ============
    //
    /**
     * Constructor
     * @param sys
     * @param ic
     * @param firsttrj
     * @param lasttrj
     */
    RDF(gcore::System *sys,
            args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj);
    /**
     * Copy Constructor
     */
    RDF(RDF &rdf);


    // DESTRUCTORS
    // ===========
    /**
     Destructor
     */
    ~RDF(void);


    // MEMBER FUNCTIONS
    // ================
    //
    /**
     * Sets the atoms to be considered as centre atoms
     */
     int addCenters(std::string s);
     /**
     * Removes all centre atoms
     */
     void clearCenters(void);
     /**
     * Sets the atoms to be considered as with atoms
     */
     int addWiths(std::string s);
    /**
     * Removes all with atoms
     */
     void clearWiths(void);
     /**
      * Calculate the rdf (all, intra- and intermolecular)
      */
     void calculateAll(void);
     /**
      * Calculate the rdf (intermolecular only)
      */
     void calculateInter(void);
     /**
      * Sets all data points of the d_rdf vector to zero
      */
     void clearRDF(void);
     /**
      * Sets the grid number for the rdf calculation to the number grid
      */
     void setGrid(unsigned int grid);
     /**
      * Sets the cutoff for the rdf calculation to cut
      */
     void setCut(double cut);
     /**
      * Prints the contents of the d_rdf vector
      */
     void print(std::ostream &os);

  };

} /* end of namespace utils */