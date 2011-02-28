namespace utils {

  class iNS;

  /**
   * Class NS (Neutron Scattering) calculates the scattering intensities ...
   *
   * And some more descriptions here please...
   *
   * @class NS
   * @ingroup utils
   * @author A. Eichenberger
   */
  class NS {

  private:

    /**
     * A pointer to the implementation class containing the data
     */
    class iNS *d_this;
    /**
     * The default constructor: should not be used since it is likely
     * to forget the setting of the ponters and iterators by hand later on...
     */
    NS(void);

  public:

    /**
     * Constructor
     */
    NS(gcore::System *sys, args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj);

    /**
     * Destructor to delete the data (implementation class)
     */
    ~NS(void);

    /**
     * Sets the number of grid points of the resulting spectrum.
     * @param grid the number of grid points
     */
    void setGrid(int grid);
    /**
     * Sets the cutoff used in the RDF calculations
     * @param cut the number of grid points
     */
    void setCut(int cut);
    /**
     * Sets the maximum Q-value of the resulting spectrum.
     * @param Qmax the maximum Q-value
     */
    void setQmax(int Qmax);
    /**
     * Sets the atoms to be considered as centre atoms
     */
    int addCenters(std::string s);
    /**
     * Sets the atoms to be considered as with atoms
     */
    int addWiths(std::string s);
    /**
     * Gets all combination of centre to with atom types and resets the lengths
     * of the corresponding vectors to that length
     */
    int getCombinations(void);
    /**
     * Checks the lengths and initialisation state of the members of NS to
     * test if it is ready for calculation of the scattering intensities
     */
    void check(void);
    /**
     * Sets a system variable to all dependant variables
     * @param sys a gromos system (gcore::System)
     */
    void setSystem(gcore::System *sys);
    /**
     * Set the iterators to the first and last trajectory file to be used
     * for the calculation.
     */
    void setTrajectories(args::Arguments::const_iterator firsttrj,
            args::Arguments::const_iterator lasttrj);
    /**
     * Prints the IACCOMBINATIONS block with all the centre-to-with iac
     * combinations to the outstream os.
     */
    void printComb(std::ostream &os);

  };

} /* end of namespace utils */
