/**
 * @file prep_xray.cc
 * Creates input-file for xray restraining
 */

#include <fstream>


/**
 * @page programs Program Documentation
 *
 * @anchor prep_xray
 * @section prep_xray Creates input-file for xray restraining
 * @author @ref ff @ref ns
 * @date 3-2-2009
 *
 * Program uses topology, cif-file (cristallographic information file), a
 * element-mapping-file, a position-file to fit, spacegroup in hermann-mauguin form,
 * cell-properties in following form: a b c alpha beta gamma, scattering resolution
 * and a standard B-Factor (used for restraining) to generate an Input-File for a
 * xray-restrained simulation.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@cif</td><td>&lt;cristallographic information file&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;gromos atom code to element mapping file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;file with atom-position to fit on&gt; </td></tr>
 * <tr><td> \@fitatom</td><td>&lt;atom specifier to filter fit&gt; </td></tr>
 * <tr><td> \@spgr</td><td>&lt;spacegroup in hermann-maauguin-form&gt; </td></tr>
 * <tr><td> \@cell</td><td>&lt;cell in form: a b c alpha beta gamma&gt; </td></tr>
 * <tr><td> \@reso</td><td>&lt;scattering resolution&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;standard b-factor&gt;</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  prep_noe
    @topo          ex.top
    @cif           ex.cif
    @pos           ex.g96
    @fitatom       1:CA
    @map           ex.
    @spgr          P 21 21 21
    @cell          50.0 51.0 52.0 90.0 90.0 90.0
    @reso          1.5
    @bfactor       1.0
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>

#include <args/Arguments.h>
#include <gio/Ginstream.h>
#include <gio/InG96.h>
#include <gcore/System.h>
#include <gio/InTopology.h>
#include <gio/StringTokenizer.h>
#include <bound/Boundary.h>
#include <args/BoundaryParser.h>
#include <utils/VirtualAtom.h>
#include <utils/Neighbours.h>
#include <ios>
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/utils/AtomSpecifier.h"

#include "../config.h"

#define HAVE_CLIPPER
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#endif

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;

// Struct for cif-data

struct cifdatatype {
  int h;
  int k;
  int l;
  float sf;
  float stddev_sf;
};

/* function that splits up a string to an array
 * (trim's spaces between)
 */
vector<string> explode(const string & tosplit) {
  std::istringstream iss(tosplit);
  vector<string> split;
  string temp;
  while (!iss.eof()) {
    if (iss.fail())
      throw gromos::Exception("getCifData", "bad line in cif-file"
            " block.");
    iss >> temp;
    split.push_back(temp);
  }
  return split;
}

/* function to read in experimental structure factors out of
 * a .cif-file.*/
vector<cifdatatype> getCifData(const string & filename) {
  string line;
  ifstream file(filename.c_str());
  if (file.fail()) {
    throw gromos::Exception("getCifData", "could not open '"+filename+"'");
  }
  vector<cifdatatype> cifdata;
  cifdatatype tcif;
  int refln_nr = 0, hnr = 0, knr = 0, lnr = 0, Fnr = 0, Fsnr = 0;
  while (true) {
    if (file.eof()) {
      break;
    }
    getline(file, line);
    if (line.substr(0, 6) == "_refln") {
      refln_nr++;
      // check index of desired data
      if (line == "_refln.index_h") {
        hnr = refln_nr;
      } else if (line == "_refln.index_k") {
        knr = refln_nr;
      } else if (line == "_refln.index_l") {
        lnr = refln_nr;
      } else if (line == "_refln.F_meas_au") {
        Fnr = refln_nr;
      } else if (line == "_refln.F_meas_sigma_au") {
        Fsnr = refln_nr;
      }
    } else {
      // Check if before or after _refln-Block
      if (refln_nr == 0) {
        // do nothing / del lines
      } else {
        // check if data or not
        vector<string> split = explode(line);
        if (split.size() < 4) {
          // no valid line...skip
        } else {
          istringstream(split[hnr - 1]) >> tcif.h;
          istringstream(split[knr - 1]) >> tcif.k;
          istringstream(split[lnr - 1]) >> tcif.l;
          istringstream(split[Fnr - 1]) >> tcif.sf;
          if (Fsnr != 0) {
            istringstream(split[Fsnr - 1]) >> tcif.stddev_sf;
          } else {
            tcif.stddev_sf = 0.0;
          }
          cifdata.push_back(tcif);
        }
      }
    }
  }
  return cifdata;
}
/* function to read in gromos atom code to element mapping and
 * save to a map of int's n strings.*/
map<int, string> getElementMapping(const string & filename) {
  Ginstream file(filename);
  vector<std::string> buffer;
  map<int, string> gactoele;
  file.getblock(buffer);
  if (buffer[0] != "FORCEFIELD")
    throw gromos::Exception("getElementMapping",
          "library file does not contain a FORCEFIELD block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("getElementMapping", "No END in FORCEFIELD"
          " block.");

  file.getblock(buffer);
  if (buffer[0] != "ELEMENTMAPPING")
    throw gromos::Exception("getElementMapping",
          "library file does not contain a ELEMENTMAPPING block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("getElementMapping", "No END in ELEMENTMAPPING"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    int gac;
    string element;
    istringstream ss(buffer[i]);
    ss >> gac >> element;
    if (ss.fail())
      throw gromos::Exception("getElementMapping", "bad line in ELEMENTMAPPING"
            " block.");
    gactoele[gac - 1] = element;
  }
  return gactoele;
}

int main(int argc, char *argv[]) {
#ifdef HAVE_CLIPPER
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@cif            <cristallographic information file>\n";
  usage += "\t@map            <gromos atom code to element mapping file>\n";
  usage += "\t@pos            <gromos atom position file to fit on>\n";
  usage += "\t@fitatom        <atom specifier to filter fit>\n";
  usage += "\t@spgr           <spacegroup in hermann-mauguin form>\n";
  usage += "\t@cell           <cell in form: a b c alpha beta gamma>\n";
  usage += "\t@reso           <scattering resolution>\n";
  usage += "\t@bfactor        <standard B-factor>\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "cif" << "map" << "pos" << "fitatom" << "spgr" << "cell" << "reso" << "bfactor";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());


    //READ IN MANUAL VALUES

    // Get Spacegroup Data
    string spgr;
    {
      Arguments::const_iterator iter = args.lower_bound("spgr");
      Arguments::const_iterator to = args.upper_bound("spgr");
      spgr = iter->second;
      for (++iter; iter != to; ++iter) {
        spgr += " ";
        spgr += iter->second;
      }
    }

    //Read in cell
    vector<float> cell;
    {
      Arguments::const_iterator iter = args.lower_bound("cell");
      Arguments::const_iterator to = args.upper_bound("cell");
      int i = 0;
      for (; iter != to; ++iter, ++i) {
        float tcell;
        std::istringstream is(iter->second);
        if (!(is >> tcell)) {
          throw gromos::Exception(argv[0],
                  "Cell parameters not numeric");
        }
        cell.push_back(tcell);
      }
      if (i != 6) {
        throw gromos::Exception(argv[0],
                "Not enough cell parameters");
      }
    }

    //read in scattering resolution
    float reso;
    {
      std::istringstream is(args["reso"]);
      if (!(is >> reso)) {
        throw gromos::Exception(argv[0],
                "Resolution parameter not numeric");
      }
    }

    //read in standard B-factor
    float bfactor;
    {
      std::istringstream is(args["bfactor"]);
      if (!(is >> bfactor)) {
        throw gromos::Exception(argv[0],
                "B-factor parameter not numeric");
      }
    }

    // Read in reference file (into a fitsystem)
    InG96 ig(args["pos"], 0, 1);
    System fitsys(sys);
    ig >> fitsys;
    AtomSpecifier fitatoms(fitsys);

    //try for fit atoms
    {
      Arguments::const_iterator iter = args.lower_bound("fitatom");
      Arguments::const_iterator to = args.upper_bound("fitatom");

      for (; iter != to; iter++) {
        fitatoms.addSpecifier(iter->second);
      }
    }

    if (fitatoms.size() == 0)
      throw gromos::Exception(argv[0],
            "Fit atom specification results in empty set of atoms");


    // Read in cif-file for generating the structure factor list (using clipper)
    vector<cifdatatype> cifdata = getCifData(args["cif"]);

    // Generate Mapping
    map<int, string> gacmapping = getElementMapping(args["map"]);
    vector<string> element;
    for (int j = 0; j < sys.numMolecules(); ++j) {
      for (int i = 0; i < sys.mol(j).numAtoms(); ++i) {
        element.push_back(gacmapping[sys.mol(j).topology().atom(i).iac()]);
      }
    }
    // Generate Solvent Mapping
    vector<string> solvent;
    for (int j=0; j<sys.numSolvents(); ++j){
      for (int i = 0; i < sys.sol(j).topology().numAtoms(); ++i) {
        solvent.push_back(gacmapping[sys.sol(j).topology().atom(i).iac()]);
      }
    }

    // OUTPUT
    // DUMMY TITLE
    cout << "TITLE\n";
    cout << "Xray specification file\n";
    cout << "END\n";

    // XRAYRESSPEC
    cout << "XRAYRESSPEC\n";
    cout << "#" << setw(4) << "H" << setw(5) << "K" << setw(5) << "L";
    cout << setw(10) << "SF" << setw(12) << "STDDEV_SF\n";
    for (unsigned int i = 0; i < cifdata.size(); ++i) {
      cout << setw(5) << cifdata[i].h;
      cout << setw(5) << cifdata[i].k;
      cout << setw(5) << cifdata[i].l;
      cout << setw(10) << cifdata[i].sf;
      cout << setw(11) << cifdata[i].stddev_sf << "\n";
    }
    cout << "END\n";

    //XRAYELEMENTSPEC
    cout << "XRAYELEMENTSPEC\n";
    cout << "# ELEMENT[1..N]\n";
    for (unsigned int i = 0; i < element.size(); ++i) {
      cout << setw(2) << element[i] << " ";
      if ((i + 1) % 20 == 0) {
        cout << "\n";
      }
    }
    cout << "\nEND\n";

    //XRAYSOLVELEMENTSPEC
    cout << "XRAYSOLVELEMENTSPEC\n";
    cout << "# ELEMENT[1..N]\n";
    for (unsigned int i = 0; i < solvent.size(); ++i) {
      cout << setw(2) << solvent[i] << " ";
      if ((i + 1) % 20 == 0) {
        cout << "\n";
      }
    }
    cout << "\nEND\n";

    //XRAYFITATOMS
    //obsolete
    /*
    cout << "XRAYFITATOMS" << endl;
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(9);
    cout << "# selecting " << fitatoms.size() << " atoms" << endl;
    cout.setf(ios::unitbuf);

    for (int i = 0; i < fitatoms.size(); ++i) {
      cout.setf(ios::right, ios::adjustfield);
      int offset = 1;
      for (int j = 0; j < ((fitatoms.mol(i) >= 0) ? fitatoms.mol(i) : sys.numMolecules()); ++j)
        offset += sys.mol(j).topology().numRes();

      cout << setw(5) << fitatoms.resnum(i) + offset;

      cout.setf(ios::left, ios::adjustfield);
      string res = fitatoms.resname(i);

      if (fitatoms.mol(i) < 0) res = "SOLV";
      cout << ' ' << setw(6) << res << setw(6) << fitatoms.name(i);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << fitatoms.gromosAtom(i) + 1
              << setw(15) << (*fitatoms.coord(i))[0]
              << setw(15) << (*fitatoms.coord(i))[1]
              << setw(15) << (*fitatoms.coord(i))[2] << endl;
    }
    cout << "END" << std::endl;
    */
    // XRAYRESPARA
    cout.precision(4);
    cout << "XRAYRESPARA\n";
    cout << "#" << setw(15) << "SPGR\n";
    cout << setw(15) << spgr << "\n";
    cout << "# " << setw(7) << "CELL" << setw(60) << "RESO" << setw(12) << "BFACTOR\n";
    cout << setw(10) << cell[0] << setw(10) << cell[1] << setw(10) << cell[2];
    cout << setw(10) << cell[3] << setw(10) << cell[4] << setw(10) << cell[5];
    cout << setw(10) << reso << setw(10) << bfactor << "\n";
    cout << "END\n";

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;

#else
      throw gromos::Exception(argv[0], "You have to compile GROMOS++ with CCP4"
              " and clipper libraries in order to use the spacegroup feature"
              " of cry.\nUse --with-ccp4 and --with-clipper for configuration.");
#endif
}