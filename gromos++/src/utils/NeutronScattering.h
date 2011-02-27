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

  public:

    /**
     * Constructor
     */
    NS(gcore::System *sys);

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
     * Gets all combination of centre to with atom types
     */
    int getCombinations(void);

    void printComb();

  };

} /* end of namespace utils */
