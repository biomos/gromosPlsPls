/**
 * @file pdb2seq.cc
 * Creates the building block sequence as well as the pdb2g96 library file
 * from a pdb file only
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pdb2seq
 * @section breas Creates the building block sequence as well as the pdb2g96 library file
 * @author @ref ae @ref bh
 * @date 18.11.2010
 *
 * Here goes the documentation...
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@pdb</td><td>&lt;pdb file&gt; </td></tr>
 * <tr><td>\@pH</dt><td>&lt;pH value of the simulation box&gt; </td></tr>
 * <tr><td>\@gff</td><td>&lt;gromos force field version&gt; </td></tr>
 * <tr><td>[\@select</dt><td>&lt;atoms to be read from PDB: \"ATOMS\" (standard), \"HETATOM\' or \"ALL\"&gt;]
 * <tr><td>[\@head</dt><td>&ltbuilding block (sequence) of head group, e.g. NH3+&gt;] </td></tr>
 * * <tr><td>[\@tail</dt><td>&ltbuilding block (sequence) of tail group, e.g. COO-&gt;] </td></tr>
 * </table>
 *
 *
 * Example using a specification file:
 * @verbatim
  cry
    @pdb     protein.pdb
    @gff     53a6
    @pH      7
 @endverbatim
 *
 * <hr>
 */

#include <iostream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>

#include "../src/args/Arguments.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/AminoAcid.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InPDB.h"


// REMEMBER
// ========
//
// - check the constants below
// - adapt the InPDB to read also gromos PDBs


// CONSTANTS
// =========
//
const double SS_cutoff = 2.2; //<<< Look up!!


 //FUNCTION DECLARATIONS
 // =====================
 //
std::vector<std::string> findSS(gio::InPDB &myPDB);
std::vector<std::string> AcidOrBase(std::vector<std::string> seq, double pH,
        utils::gromosAminoAcidLibrary &gaal);
void writeResSeq(std::ostream &os, std::vector<std::string> seq);


using namespace std;
using namespace args;
using namespace gio;
using namespace utils;


int main(int argc, char **argv) {


  // DEFINE THE COMMAND LINE ARGUMENTS AND USAGE STRING
  // ==================================================
  //
  Argument_List knowns;
  knowns << "pdb" << "gff" << "pH" << "aalib" << "select" << "head" << "tail"
          << "develop";
  //
  string usage = "# " + string(argv[0]);
  usage += "\n\t@pdb      <pdb file to be read>\n";
  usage += "\t@pH       <specification file for the symmetry transformations]>\n";
  usage += "\t[@select  <atoms to be read from PDB: \"ATOMS\" (standard), \"HETATOM\' or \"ALL\">]\n";
  usage += "\t[@aalib   <amino acid library file>]\n";
  usage += "\t[@gff     <GROMOS force field version (if no @aalib sepcified): 45A4, 53A6>]\n";
  usage += "\t[@head    [<building block (sequence) of head group, e.g. NH3+>]]\n";
  usage += "\t[@tail    [<building block (sequence) of tail group, e.g. COO->]]\n";

  try {

    // READ THE COMMAND LINE ARGUMENTS
    // ===============================
    //
    Arguments args(argc, argv, knowns, usage);
    //
    // which PDB file to be used (file name)?
    if (args.count("pdb") != 1) {
      throw gromos::Exception(argv[0], "specify exactly one pdb file (@pdb)");
    }
    //
    // the agromos force fiel version to be used
    string gffversion;
    if (args.count("gff") == 1) {
      gffversion = args.find("gff")->second;
    } else {
      throw gromos::Exception(argv[0], "specify exactly one GROMOS force field version (@gff)");
    }
    //
    // simulation intended to run at which pH value?
    double pH;
    {
      stringstream ss, css;
      if (args.count("pH") == 1) {
        ss.str(args.find("pH")->second);
        css.str(ss.str());
        ss >> pH;
        if (ss.fail() || ss.bad()) {
          stringstream msg;
          msg << "could not convert " << css.str() << " to a valid pH value";
          throw gromos::Exception(argv[0], msg.str());
        }
      } else {
        throw gromos::Exception(argv[0], "no or more than one value indicated as pH (@pH)");
      }
    }

    // selection of atoms read from PDB
    string select = "ATOM";
    if(args.count("select") > 0) {
      select = args["select"];
      if(select != "ATOM" && select != "HETATOM" && select != "ALL") {
        stringstream msg;
        msg << select << " is not a proper selection of atoms to be read from pdb"
                " (@select), allowed is \"ATTOM\", \"HETATOM\" or \"ALL\"";
        throw gromos::Exception(argv[0], msg.str());
      }
    }

    // REMOVE THIS LATER
    if(args.count("develop") < 0) {
      throw gromos::Exception("PROGRAM UNDER DEVELOPMENT", "do not use this program yet");
    }

    // READ THE PDB FILE
    // =================
    //
    InPDB ipdb(args["pdb"]);
    ipdb.select(select);
    ipdb.read();

    // BUILD/READ THE LIBRARY FILES
    // ============================
    //
    utils::gromosAminoAcidLibrary gaal;
    gaal.loadHardcoded45A4();
    //gaal.writeLibrary(cout, "  45A4 AMINO ACID LIBRARY");

    // RESIDUE SEQUENCE FROM PDB
    // =========================
    //
    // extract the PDB residue sequence
    vector<string> resSeq = ipdb.getResSeq();
    // check for disulfide briges
    resSeq = findSS(ipdb);
    // adapt the protonation state of the residues
    resSeq = AcidOrBase(resSeq, pH, gaal);

    // write the (transformed) residue sequence
    writeResSeq(cout, resSeq);
    // add head/tail group and do other corrections (if necessary)
    //

    // HISTIDIN SHIT
    // =============
    //
    // decide about the His protonation state


    // PRINT OUT ALL WARNINGS/ERRORS
    // =============================
    //
    //

    // WRITE THE SEQUENCE AND LIBRARIES
    // ================================
    //
    // - residue sequence
    // - pdb2g96 library
    // - "corrected" PDB
    //

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



//FUNCTION DEFINITIONS

/**
 * Check for S-S bridges
 */
std::vector<std::string> findSS(InPDB &myPDB){
  vector<string> sequence = myPDB.getResSeq();
  for (unsigned int i = 0; i < myPDB.numAtoms(); ++i) {
    if (sequence[myPDB.getResNumber(i) - 1] == "CYS") {
      if (myPDB.getAtomName(i) == "SG") {
        gmath::Vec first_cys = myPDB.getAtomPos(i);
        for (unsigned int j = i + 1; j < myPDB.numAtoms(); ++j) {
          if (sequence[myPDB.getResNumber(j) - 1] == "CYS") {
            if (myPDB.getAtomName(j) == "SG") {
              gmath::Vec second_cys = myPDB.getAtomPos(j);
              double dist = (first_cys - second_cys).abs();
              if (dist < SS_cutoff) {
                sequence[myPDB.getResNumber(i) - 1] = "CYS1";
                sequence[myPDB.getResNumber(j) - 1] = "CYS2";
              }
            }
          }
        }
      }
    }
  }
  return sequence;
}

vector<string> AcidOrBase(vector<string> seq, double pH, utils::gromosAminoAcidLibrary &gaal) {
  for(int i = 0; i < seq.size(); ++i) {
    if(seq[i] != "CYS1" && seq[i] != "CYS2") {
      double pKc = gaal.pKc(seq[i]);
      if(pKc > 0 && pH > pKc) {
        seq[i] = gaal.pdb2base(seq[i]);
      } else {
        seq[i] = gaal.pdb2acid(seq[i]);
      }
    }
  }
  return seq;
}

void writeResSeq(std::ostream &os, std::vector<std::string> seq) {
  os << "RESSEQUENCE";
  for (unsigned int i = 0; i < seq.size(); i++) {
    os << setw(6);
    if (i % 10 == 0) {
      os << endl;
    }
    os << seq[i];
  }
  os << "\nEND\n";
}
