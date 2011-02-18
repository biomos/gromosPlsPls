namespace utils {


  struct gromosAminoAcid {
    std::string acid;
    std::string base;
    std::multimap<std::string, std::string> Hdonors;
    std::multimap<std::string, std::string> Hacceptors;
    double pKa;
    double pKb;
    double pKc;
  };

  class gromosAminoAcidLibrary {

  private:
    std::string version;
    std::map<std::string, gromosAminoAcid> lib;

  public:
    void load(std::string &fname);
    void loadHardcoded45A4(void);
    void loadHardcoded53A6(void);
    void writeLibrary(std::ostream &os, std::string title = "");

  };

}