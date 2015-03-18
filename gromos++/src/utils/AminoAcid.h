#ifndef INCLUDED_AMINOACID
#define INCLUDED_AMINOACID

namespace utils {

  struct gromosAminoAcid {
    std::string acid;
    std::string base;
    std::map<std::string, std::vector<std::string> > Hdonors;
    std::map<std::string, std::vector<std::string> > Hacceptors;
    double pKa;
    double pKb;
    double pKc;
  };

  class gromosAminoAcidLibrary {

  private:
    std::string version;
    std::map<std::string, gromosAminoAcid> lib;
    //std::map<std::string, std::vector<std::string> > Hdonors;
    //std::map<std::string, std::vector<std::string> > Hacceptors;

  public:
    void load(std::string &fname);
    void loadHardcoded45A4(void);
    void loadHardcoded53A6(void);
    void writeLibrary(std::ostream &os, std::string title = "");
    std::string pdb2acid(std::string PDBname);
    std::string pdb2base(std::string PDBname);
    double pKa(std::string PDBname);
    double pKb(std::string PDBname);
    double pKc(std::string PDBname);
    std::vector<std::string> rHdonors(std::string PDBname, std::string GROMOSname);
    std::vector<std::string> rHacceptors(std::string PDBname, std::string GROMOSname);
    std::map<std::string, gromosAminoAcid> getAminoAcids();
  };

}
#endif
