/**
 * @file prep_xray.cc
 * Creates input-file for xray restraining
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor prep_xray
 * @section prep_xray Creates input-file for X-ray restraining
 * @author @ref ff @ref ns
 * @date 3-2-2009
 *
 * Program prep_xray converts a crystallographic information file (CIF) containing
 * reflection data into a GROMOS X-ray restraints specification file. Using a
 * mapping file (\@map) it writes out the element names of the solute and solvent
 * atoms according to their integer atom codes.
 * The spacegroup, cell, resolution and standard B-factor of the atoms are written
 * out the outfile in the correct format.
 * The reflection list can be filtered according to the given resolution range.
 * If no resolution is given, it is determined automatically.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@cif</td><td>&lt;cristallographic information file&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;file with IAC-to-elementname mapping&gt; </td></tr>
 * <tr><td> \@spacegroup</td><td>&lt;spacegroup in Hermann-Maauguin format&gt; </td></tr>
 * <tr><td> \@cell</td><td>&lt;cell in form: a b c alpha beta gamma&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution, from and to&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;standard B-factor (in @f$\mathrm{nm}^2@f$)&gt;</td></tr>
 * </table>
 *
 * For the format of the mapping file see the @ref structure_factor structure_factor
 * program documentation.
 *
 * Example:
 * @verbatim
  prep_noe
    @topo          ex.top
    @cif           ex.cif
    @map           ex.map
    @spacegroup    P 21 21 21
    @cell          5.00 5.10 5.20 90.0 90.0 90.0
    @resolution    0.3 0.15
    @bfactor       0.01
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <streambuf>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>
#include <fstream>
#include <limits>
#include <ios>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gio/InCIF.h"
#include "../src/gio/InIACElementNameMapping.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace std;

// Additional Clipper Headers
#include "../config.h"
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>

int main(int argc, char *argv[]) {
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@cif            <cristallographic information file>\n";
  usage += "\t@map            <gromos atom code to element mapping file>\n";
  usage += "\t@spacegroup     <spacegroup in hermann-mauguin form>\n";
  usage += "\t@cell           <cell in form: a b c alpha beta gamma>\n";
  usage += "\t@resolution     <scattering resolution min max>\n";
  usage += "\t@bfactor        <standard B-factor>\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "cif" << "map" << "spacegroup" << "cell" << "resolution" << "bfactor";

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
      Arguments::const_iterator iter = args.lower_bound("spacegroup");
      Arguments::const_iterator to = args.upper_bound("spacegroup");
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
    // create the clipper cell
    clipper::Cell clipperCell(
            clipper::Cell_descr(cell[0]*10.0f, cell[1]*10.0f, cell[2]*10.0f,
                cell[3], cell[4], cell[5]));

    //read in scattering resolution
    bool have_reso = false;
    double reso_min = 0.0, reso_max = numeric_limits<double>::max();
    {
      if (args.count("resolution") == 2) {
        have_reso = true;
        Arguments::const_iterator iter = args.lower_bound("resolution");
        std::istringstream is(iter->second);
        if (!(is >> reso_min)) {
          throw gromos::Exception(argv[0],
                  "Resolution minimum not numeric");
        }
        ++iter;

        is.clear();
        is.str(iter->second);
        if (!(is >> reso_max)) {
          throw gromos::Exception(argv[0],
                  "Resolution maximum not numeric");
        }
      } // two resolutions
      else if (args.count("resolution") == -1) {
        have_reso = false;
      } else {
        throw gromos::Exception(argv[0],
                  "Please give a resolution range or let the program determine it");
      }
    }

    // swap them if someone is confused by the fact that the low numeric
    // value is actually a high resolution
    if (have_reso && reso_min <= reso_max) {
      const double tmp = reso_min;
      reso_min = reso_max;
      reso_max = tmp;
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

    // Read in cif-file for generating the structure factor list (using clipper)
    InCIF ciffile(args["cif"]);
    vector<CIFData> cifdata = ciffile.getData();

    // Generate Mapping
    InIACElementNameMapping mapfile(args["map"]);
    map<int, string> gacmapping = mapfile.getData();
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

    int filtered = 0;
    for (unsigned int i = 0; i < cifdata.size(); ++i) {
      clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
      double hkl_resolution = sqrt(1.0 / hkl.invresolsq(clipperCell)) / 10.0;
      if (have_reso) {
        if (hkl_resolution < reso_max || hkl_resolution > reso_min) {
          ++filtered;
          continue;
        }
      } else {
        // determine maximal resolution (minimal numerical number)
        reso_max = std::min(reso_max, hkl_resolution);
      }
      cout << setw(5) << cifdata[i].h;
      cout << setw(5) << cifdata[i].k;
      cout << setw(5) << cifdata[i].l;
      cout << setw(10) << cifdata[i].f_obs;
      cout << setw(11) << cifdata[i].stddev_f_obs << "\n";
    }
    cout << "END\n";

    if (filtered)
      cerr << filtered << " reflections have been filtered away." << endl;

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

    // XRAYRESPARA
    cout.precision(4);
    cout << "XRAYRESPARA\n";
    cout << "#" << setw(15) << "SPGR\n";
    cout << setw(15) << spgr << "\n";
    cout << "# " << setw(7) << "CELL" << setw(60) << "RESO" << setw(12) << "BFACTOR\n";
    cout << setw(10) << cell[0] << setw(10) << cell[1] << setw(10) << cell[2];
    cout << setw(10) << cell[3] << setw(10) << cell[4] << setw(10) << cell[5];
    cout << setw(10) << reso_max << setw(10) << bfactor << "\n";
    cout << "END\n";

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

#else
int main(int argc, char **argv) {
  cerr << "You have to compile GROMOS++ with CCP4"
              " and clipper libraries in order to use this program.\n"
              "Use --with-ccp4 and --with-clipper for configuration."
          << endl;
  return 1;
}
#endif
