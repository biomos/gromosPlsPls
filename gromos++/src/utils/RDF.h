namespace utils {

  // the implementation class, just to let the compiler know that it exists
  class iRDF;

   /**
   * Class RDF:
   * a class to calculate radial distribution functions of a given system
   * described by a topology file and the corresponding trajectory files.
   *
   * The radial distribution function, g(r), is defined as
   *
   * @f[ g(r) = \frac{N_J(k)}{4\pi r^2 dr \rho_J} @f]
   *
   * which is the probability of finding a particle of type J at distance r
   * from a central particle of type I, relative to the probability for a homogenous
   * distribution of particle of type J around particle of type I.
   *
   * @class RDF
   * @ingroup utils
   * @author A. Eichenberger
   * */
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
     * Standard constructor, don't forget to set the system and trajectories
     * by hand!
     */
    RDF();
    /**
     * Constructor to initialize the class.
     * @param sys The system
     * @param firsttrj An iterator defining the first trajectory file
     * @param lasttrj An iterator defining the last trajector file
     */
    RDF(gcore::System *sys,
            args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj);
    /**
     * Constructor
     */
    RDF(const RDF &rdf);


    // DESTRUCTORS
    // ===========
    /**
     * Destructor
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
      * Adds an atom to the centre atoms
      */
     void addCentersAtom(int m, int a);
     /**
     * Removes all centre atoms
     */
     void clearCenters(void);
     /**
     * Sets the atoms to be considered as with atoms
     */
     int addWiths(std::string s);
     /**
      * Adds an atom to the with atoms
      */
     void addWithAtom(int m, int a);
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
      * Sets a system to the class
      */
     void setSystem(gcore::System *sys);
     /**
     * Set the iterators to the first and last trajectory file to be used
     * for the calculation.
     */
    void setTrajectories(args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj);
     /**
      * Prints the contents of the d_rdf vector
      */
     void print(std::ostream &os);

  };

} /* end of namespace utils */