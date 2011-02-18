#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>

#include "AminoAcid.h"

using namespace std;

namespace utils {

  void gromosAminoAcidLibrary::loadHardcoded45A4(void) {

    // clear the library if not empty
    if(lib.size() > 0) {
      lib.clear();
    }

    string pdbname;
    gromosAminoAcid gaa;

    // set the library version
    version = "45A4";

    pdbname = "ALA";
    gaa.acid = "ALA";
    gaa.base = "ALA";
    gaa.pKa = 2.33;
    gaa.pKb = 9.71;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("ALA", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("ALA", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ARG";
    gaa.acid = "ARG";
    gaa.base = "ARGN";
    gaa.pKa = 2.03;
    gaa.pKb = 9.00;
    gaa.pKc = 12.10;
    gaa.Hdonors.insert(pair<string, string > ("ARG", "N"));
    gaa.Hdonors.insert(pair<string, string > ("ARG", "NE"));
    gaa.Hdonors.insert(pair<string, string > ("ARG", "NH1"));
    gaa.Hdonors.insert(pair<string, string > ("ARG", "NH2"));
    gaa.Hacceptors.insert(pair<string, string > ("ARG", "O"));
    gaa.Hdonors.insert(pair<string, string > ("ARGN", "N"));
    gaa.Hdonors.insert(pair<string, string > ("ARGN", "NE"));
    gaa.Hdonors.insert(pair<string, string > ("ARGN", "NH1"));
    gaa.Hdonors.insert(pair<string, string > ("ARGN", "NH2"));
    gaa.Hacceptors.insert(pair<string, string > ("ARGN", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("ARGN", "NH1")); // check again later
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ASN";
    gaa.acid = "ASN";
    gaa.base = "ASN";
    gaa.pKa = 2.16;
    gaa.pKb = 8.73;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("ASN", "N"));
    gaa.Hdonors.insert(pair<string, string > ("ASN", "ND2"));
    gaa.Hacceptors.insert(pair<string, string > ("ASN", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("ASN", "OD1"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ASP";
    gaa.acid = "ASPH";
    gaa.base = "ASP";
    gaa.pKa = 1.95;
    gaa.pKb = 9.66;
    gaa.pKc = 3.71;
    gaa.Hdonors.insert(pair<string, string > ("ASP", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("ASP", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("ASP", "OD1"));
    gaa.Hacceptors.insert(pair<string, string > ("ASP", "OD2"));
    gaa.Hdonors.insert(pair<string, string > ("ASPH", "N"));
    gaa.Hdonors.insert(pair<string, string > ("ASPH", "OD2"));
    gaa.Hacceptors.insert(pair<string, string > ("ASPH", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("ASPH", "OD2")); // check later again
    gaa.Hacceptors.insert(pair<string, string > ("ASPH", "OD1"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "CYS";
    gaa.acid = "CYSH";
    gaa.base = "CYS";
    gaa.pKa = 1.91;
    gaa.pKb = 10.28;
    gaa.pKc = 8.14;
    gaa.Hdonors.insert(pair<string, string > ("CYS", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("CYS", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("CYS", "SG"));
    gaa.Hdonors.insert(pair<string, string > ("CYSH", "N"));
    gaa.Hdonors.insert(pair<string, string > ("CYSH", "SG"));
    gaa.Hacceptors.insert(pair<string, string > ("CYSH", "O"));
    gaa.Hdonors.insert(pair<string, string > ("CYS1", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("CYS1", "O"));
    gaa.Hdonors.insert(pair<string, string > ("CYS2", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("CYS2", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    pdbname = "GLN";
    gaa.acid = "GLN";
    gaa.base = "GLN";
    gaa.pKa = 2.18;
    gaa.pKb = 9.00;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("GLN", "N"));
    gaa.Hdonors.insert(pair<string, string > ("GLN", "NE2"));
    gaa.Hacceptors.insert(pair<string, string > ("GLN", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("GLN", "OE"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "GLU";
    gaa.acid = "GLUH";
    gaa.base = "GLU";
    gaa.pKa = 2.16;
    gaa.pKb = 9.58;
    gaa.pKc = 4.15;
    gaa.Hdonors.insert(pair<string, string > ("GLU", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("GLU", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("GLU", "OE1"));
    gaa.Hacceptors.insert(pair<string, string > ("GLU", "OE2"));
    gaa.Hdonors.insert(pair<string, string > ("GLUH", "N"));
    gaa.Hdonors.insert(pair<string, string > ("GLUH", "OE2"));
    gaa.Hacceptors.insert(pair<string, string > ("GLUH", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("GLUH", "OE2")); // check later again
    gaa.Hacceptors.insert(pair<string, string > ("GLUH", "OE1"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "GLY";
    gaa.acid = "GLY";
    gaa.base = "GLY";
    gaa.pKa = 2.34;
    gaa.pKb = 9.58;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("GLY", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("GLY", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "HIS";
    gaa.acid = "HISH";
    gaa.base = "HISX";
    gaa.pKa = 1.70;
    gaa.pKb = 9.09;
    gaa.pKc = 6.04;
    gaa.Hdonors.insert(pair<string, string > ("HISH", "N"));
    gaa.Hdonors.insert(pair<string, string > ("HISH", "ND1"));
    gaa.Hdonors.insert(pair<string, string > ("HISH", "NE2"));
    gaa.Hacceptors.insert(pair<string, string > ("HISH", "O"));
    gaa.Hdonors.insert(pair<string, string > ("HISA", "N"));
    gaa.Hdonors.insert(pair<string, string > ("HISA", "ND1"));
    gaa.Hacceptors.insert(pair<string, string > ("HISA", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("HISA", "NE2"));
    gaa.Hdonors.insert(pair<string, string > ("HISB", "N"));
    gaa.Hdonors.insert(pair<string, string > ("HISB", "NE2"));
    gaa.Hacceptors.insert(pair<string, string > ("HISB", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("HISB", "ND1"));
    gaa.Hdonors.insert(pair<string, string > ("HISX", "N")); // come back here and have fun ;-)
    gaa.Hacceptors.insert(pair<string, string > ("HISX", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "ILE";
    gaa.acid = "ILE";
    gaa.base = "ILE";
    gaa.pKa = 2.26;
    gaa.pKb = 9.60;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("ILE", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("ILE", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "LEU";
    gaa.acid = "LEU";
    gaa.base = "LEU";
    gaa.pKa = 2.32;
    gaa.pKb = 9.58;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("LEU", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("LEU", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "LYS";
    gaa.acid = "LYSH";
    gaa.base = "LYS";
    gaa.pKa = 2.15;
    gaa.pKb = 9.16;
    gaa.pKc = 10.67;
    gaa.Hdonors.insert(pair<string, string > ("LYS", "N"));
    gaa.Hdonors.insert(pair<string, string > ("LYS", "NZ"));
    gaa.Hacceptors.insert(pair<string, string > ("LYS", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("LYS", "NZ")); // check later again
    gaa.Hdonors.insert(pair<string, string > ("LYSH", "N"));
    gaa.Hdonors.insert(pair<string, string > ("LYSH", "NZ"));
    gaa.Hacceptors.insert(pair<string, string > ("LYSH", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "MET";
    gaa.acid = "MET";
    gaa.base = "MET";
    gaa.pKa = 2.16;
    gaa.pKb = 9.08;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("MET", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("MET", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("MET", "SD"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    gaa.acid = "PHE";
    gaa.base = "PHE";
    gaa.pKa = 2.18;
    gaa.pKb = 9.09;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("PHE", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("PHE", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "PRO";
    gaa.acid = "PRO";
    gaa.base = "PRO";
    gaa.pKa = 1.95;
    gaa.pKb = 10.47;
    gaa.pKc = -1.0;
    gaa.Hacceptors.insert(pair<string, string > ("PRO", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "SER";
    gaa.acid = "SER";
    gaa.base = "SER";
    gaa.pKa = 2.13;
    gaa.pKb = 9.05;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("SER", "N"));
    gaa.Hdonors.insert(pair<string, string > ("SER", "OG"));
    gaa.Hacceptors.insert(pair<string, string > ("SER", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("SER", "OH")); // think about later
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "THR";
    gaa.acid = "THR";
    gaa.base = "THR";
    gaa.pKa = 2.20;
    gaa.pKb = 8.96;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("THR", "N"));
    gaa.Hdonors.insert(pair<string, string > ("THR", "OG1"));
    gaa.Hacceptors.insert(pair<string, string > ("THR", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("THR", "OG1")); // think about later
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "TRP";
    gaa.acid = "TRP";
    gaa.base = "TRP";
    gaa.pKa = 2.38;
    gaa.pKb = 9.34;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("TRP", "N"));
    gaa.Hdonors.insert(pair<string, string > ("TRP", "NH1"));
    gaa.Hacceptors.insert(pair<string, string > ("TRP", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "TYR";
    gaa.acid = "TYR";
    gaa.base = "TYR";
    gaa.pKa = 2.24;
    gaa.pKb = 9.04;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("TYR", "N"));
    gaa.Hdonors.insert(pair<string, string > ("TYR", "OH"));
    gaa.Hacceptors.insert(pair<string, string > ("TYR", "O"));
    gaa.Hacceptors.insert(pair<string, string > ("TYR", "OH")); // think about later
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();
    pdbname = "VAL";
    gaa.acid = "VAL";
    gaa.base = "VAL";
    gaa.pKa = 2.27;
    gaa.pKb = 9.52;
    gaa.pKc = -1.0;
    gaa.Hdonors.insert(pair<string, string > ("VAL", "N"));
    gaa.Hacceptors.insert(pair<string, string > ("VAL", "O"));
    lib.insert(pair<string, gromosAminoAcid > (pdbname, gaa));
    gaa.Hdonors.clear();
    gaa.Hacceptors.clear();

  }


  void gromosAminoAcidLibrary::writeLibrary(ostream &os, string title) {

    // only wirte a library if there is something to write
    if (lib.size() > 0) {

      os << "TITLE\n";
      os << title << endl;
      os << "END\n";
      os << "VERSION\n";
      os << version << endl;
      os << "END\n";
      for (map<string, gromosAminoAcid>::iterator it = lib.begin();
              it != lib.end(); it++) {
        os << "AMINOACID\n";
        os << "#" << setw(11) << "PDB name" << setw(12) << "acid name" << setw(12)
                << "base name" << endl;
        os << setw(12) << it->first
                << setw(12) << it->second.acid
                << setw(12) << it->second.base << endl;
        os << "#\n";
        os << "#" << setw(11) << "pKa" << setw(12) << "pKb" << setw(12) << "pKc" << endl;
        os << setw(12) << it->second.pKa << setw(12) << it->second.pKb
                << setw(12) << it->second.pKc << endl;
        os << "#" << setw(11) << "residue" << setw(12) << "H-donors" << endl;
        os << "END\n#\n";
      }

    } else {

    }

  }

}
