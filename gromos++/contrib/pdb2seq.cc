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
 * <tr><td>\@gff</td><td>&lt;gromos force field version&gt; </td></tr>
 * <tr><td>\@pH</dt><td>&lt;pH value of the simulation box&gt; </td></tr>
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

#include "../src/args/Arguments.h"
#include "../src/gromos/Exception.h"
#include "../src/utils/AminoAcid.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InPDB.h"


/*
 * FUNCTION DECLARATIONS
 * =====================
 */

std::vector<std::string> findSS(gio::InPDB myPDB);


using namespace std;
using namespace args;
using namespace gio;
using namespace utils;


int main(int argc, char **argv) {


  // DEFINE THE COMMAND LINE ARGUMENTS AND USAGE STRING
  // ==================================================
  //
  Argument_List knowns;
  knowns << "pdb" << "gff" << "pH" << "select" << "head" << "tail" << "develop";
  //
  string usage = "# " + string(argv[0]);
  usage += "\n\t@pdb      <pdb file>\n";
  usage += "\t@gff      <GROMOS force field version (45a4, 53a6, ...)>\n";
  usage += "\t@pH       <specification file for the symmetry transformations]>\n";
  usage += "\t[@select  <atoms to be read from PDB: \"ATOMS\" (standard), \"HETATOM\' or \"ALL\">]\n";
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

    // RESIDUE SEQUENCE FROM PDB
    // =========================
    //
    vector<string> resSeq = ipdb.getResSeq();
    //
    // check the extracted residue sequence
    //   - disulfide bridges
    //
    //   - add head/tail group and do other corrections (if necessary)
    //

    // DECIDE ABOUT THE PROTONATION STATE OF THE RESIDUES
    // ==================================================
    //
    // 1. Loop over residie sequence:
    //    PDB --> GROMOS residue names (acid/base/XXX)
    //
    // 2. Loop over residue sequence:
    //    Dicide about special cases (e.g. His)
    //

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
std::vector<std::string> findSS(gio::InPDB myPDB){

  std::vector<std::string> sequence = myPDB.getResSeq();

  for (unsigned int i=0; i<myPDB.numAtoms(); ++i){

  }

}