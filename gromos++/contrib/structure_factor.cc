/**
 * @file structure_factor.cc
 * calculates structure factors
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor structure_factor
 * @section structure_factor calculates structure factors
 * @author @ref ns ff
 * @date 8.4.2009
 *
 * Program structure factor calculate X-ray reflection structure factor amplitudes
 * and phases from a given trajectory. Only the atoms given by the @ref AtomSpecifier
 * \@atomssf are considered for the calculation. The atoms' IAC are mapped to their
 * element names according to the rules given in the \@map file. The atoms' B-factors
 * and occupancies are read from a special file (\@befactor) if requested or default
 * to @f$0.01 \mathrm{nm}^2@f$ and 100%.
 * Structure factors are calculated to the given resolution (\@resultion) while
 * the cell information is calculated from the system's box.
 * Symmetry operations are taken into account by specifing a (\@spacegroup).
 * Make sure you only give asymetric unit when using \@spacegroup.
 *
 * Format of the IAC to element-name mapping file:
 * @verbatim
TITLE
map IAC to element name for xxxx force-field
END
ELEMENTMAPPING
# IAC ELEMENTNAME
1 O
11 N
12 C
26 Fe
END
@endverbatim
 *
 * Format of the B-factor and occupancy file:
 * The B-factors have to be given @f\mathrm{nm}^2@f$.
 * @verbatim
TITLE
B-factors and occupancies for all atoms
END
BFACTOROCCUPANCY
# B-factor Occupancy
0.01  1.0
0.02  0.8
END
   @endverbatim
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomssf</td><td>&lt;@ref AtomSpecifier: atoms to consider for structure_factor&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;file with IAC-to-elementname mapping&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;file with experimental B-factors and occupancies&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution [nm]&gt; </td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup in Hermann-Mauguin format, default: P 1&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  structure_factor
    @topo       ex.top
    @time       0 0.1
    @atomssf    1:CA
    @traj       ex.tr
    @map        ex.map
    @bfactor    ex.bfc
    @resolution 0.1
    @spacegroup P 21 21 21

    @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <memory>

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Time.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Physics.h"

#undef DEBUG

#ifdef NDEBUG
#define DEBUG(level, x)
#else
#define DEBUG_LEVEL 5
#define DEBUG(level, x) if (level <= DEBUG_LEVEL) std::cout << "DEBUG: " << x << std::endl
#endif

// Additional Clipper Headers
#include "../config.h"
#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>
#include <clipper/clipper-contrib.h>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace args;
using namespace std;
using namespace gmath;

// Struct for cif-data

struct cifdatatype {
  int h;
  int k;
  int l;
  float sf;
  float uk;
};

// Struct for BFactors and Occupancy
struct bfoccudatatype{
  float bf;
  float oc;
};

/* function to read in gromos atom code to element mapping and
 * save to a map of int's n strings.*/
std::map<int, std::string> getElementMapping(const std::string & filename) {
  gio::Ginstream file(filename);
  std::vector<std::string> buffer;
  std::map<int, std::string> gactoele;

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
    std::string element;
    std::istringstream ss(buffer[i]);
    ss >> gac >> element;
    if (ss.fail())
      throw gromos::Exception("getElementMapping", "bad line in ELEMENTMAPPING"
            " block.");
    gactoele[gac - 1] = element;
  }
  return gactoele;
}

/* function to read in experimental B-factors out of
 * a .bfc-file (used perl-script to generate this out
 * of a pdb-file.*/
std::vector<bfoccudatatype> getBfactors(const std::string & filename) {
  gio::Ginstream file(filename);
  std::vector<std::string> buffer;
  std::vector<bfoccudatatype> bfoccu;
  bfoccudatatype sdata;
  // Get Block
  file.getblock(buffer);
  if (buffer[0] != "BFACTOROCCUPANCY")
    throw gromos::Exception("getBfactors",
          "library file does not contain a BFACTOROCCUPANCY block!");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("getBfactors", "No END in BFACTOROCCUPANCY"
          " block.");
  // read in the lib
  for (size_t i = 1; i < buffer.size() - 1; i++) {
    // Data-Vars
    std::istringstream ss(buffer[i]);
    ss >> sdata.bf >> sdata.oc;
    if (ss.fail())
      throw gromos::Exception("getBfactors", "bad line in BFACTOROCCUPANCY"
            " block.");
    // Push back to vetor
    bfoccu.push_back(sdata);
  }
  return bfoccu;
}

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "traj" << "map" << "atomssf" << "time" << "bfactor"
          << "resolution" << "spacegroup";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atomssf    <atomspecifier: atoms to consider for structure_factor>\n";
  usage += "\t@traj       <trajectory files>\n";
  usage += "\t@map        <IAC-to-ElementName map-file>\n";
  usage += "\t[@bfactor    <experimental B-factors>]\n";
  usage += "\t@resolution <scattering resolution>\n";
  usage += "\t[@spacegroup <spacegroup in Hermann-Maugin format, default: P 1>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);

  try {
    Arguments args(argc, argv, knowns, usage);

    // Hardcoded B-factor conversion factor.
    const float sqpi2=(gmath::pi*gmath::pi*8.0/3.0);

    // Get Spacegroup Data or default to no symmetry (P 1)
    std::string spgrdata("P 1");
    {
      Arguments::const_iterator iter = args.lower_bound("spacegroup");
      Arguments::const_iterator to = args.upper_bound("spacegroup");
      if (iter != to) {
        spgrdata = iter->second;
        for (++iter; iter != to; ++iter) {
          spgrdata += " ";
          spgrdata += iter->second;
        }
      }
    }
    // initialize the spacegroup
    std::auto_ptr<clipper::Spgr_descr> spgrinit;
    try {
     spgrinit = std::auto_ptr<clipper::Spgr_descr>(new clipper::Spgr_descr(spgrdata, clipper::Spgr_descr::HM));
    } catch(clipper::Message_fatal & msg) {
      throw gromos::Exception(argv[0], "Invalid spacegroup: " + msg.text());
    }
    clipper::CSpacegroup spgr(clipper::String("base spgr"), clipper::Spacegroup(*spgrinit));

    // Get resolution as a double
    double resolution;
    {
      std::istringstream is(args["resolution"]);
      if (!(is >> resolution)) {
        throw gromos::Exception(argv[0],
                "Resolution parameter not numeric");
      }
    }

    // get simulation time
    Time time(args);
    // read topology
    InTopology it(args["topo"]);
    // System for calculation
    System sys(it.system());
    AtomSpecifier calcatoms(sys);
    //get structure_factor atoms
    {
      Arguments::const_iterator iter = args.lower_bound("atomssf");
      Arguments::const_iterator to = args.upper_bound("atomssf");

      for (; iter != to; iter++) {
        string spec = iter->second.c_str();
        calcatoms.addSpecifier(spec);
      }
    }
    if (calcatoms.size() == 0)
      throw gromos::Exception(argv[0], "No structure_factor-atoms specified!");
    
    //Get gac-to-ele mapping
    std::map<int, std::string> gacmapping = getElementMapping(args["map"]);

    // Get experimental Bfactors and occupancy
    std::vector<bfoccudatatype> bfoc;
    bool has_bfactor = false;
    if (args.count("bfactor") == 1) {
      bfoc = getBfactors(args["bfactor"]);
      has_bfactor = true;
    }

    //===========================
    // loop over all trajectories
    InG96 ic;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        if (!sys.hasPos)
          throw gromos::Exception(argv[0],
                "Unable to read POSITION(RED) block from "
                "trajectory file.");
        
        if (!sys.hasBox)
          throw gromos::Exception(argv[0],
                "Cannot calculate structure factors without a box.");

        // create the cell
        clipper::Cell_descr cellinit(
                sys.box().K().abs()*10.0f,
                sys.box().L().abs()*10.0f,
                sys.box().M().abs()*10.0f,
                sys.box().alpha(),
                sys.box().beta(),
                sys.box().gamma());
        clipper::CCell cell(spgr, clipper::String("base cell"), clipper::Cell(cellinit));

        // create the resolutions and corresponding lattice
        clipper::CResolution reso(cell, clipper::String("base reso"), clipper::Resolution(float(resolution * 10.0)));
        clipper::CHKL_info hkls(reso, clipper::String("base hkls"), true);
        clipper::CHKL_data<clipper::data32::F_phi> fphi(hkls);

        // Fill Clipper Atom list
        // we do this insight the loop due to solvent molecules!
        std::vector<clipper::Atom> atomvec;
        for (int i = 0; i < calcatoms.size(); i++) {
          clipper::Atom atm;
          // convert to angstrom
          atm.set_coord_orth(clipper::Coord_orth(
                  calcatoms.pos(i)[0] * 10.0,
                  calcatoms.pos(i)[1] * 10.0,
                  calcatoms.pos(i)[2] * 10.0));

          if (has_bfactor) {
            if (i >= int(bfoc.size())) {
              throw gromos::Exception("structre_factor", "Not enough B-factors given");
            }
            atm.set_occupancy(bfoc[i].oc);
            // convert to Angstrom^2
            atm.set_u_iso(bfoc[i].bf * 100.0 / sqpi2);
          } else {
            atm.set_occupancy(1.0);
            atm.set_u_iso(1.0 / sqpi2);
          }
          atm.set_element(gacmapping[calcatoms.iac(i)]);
          atomvec.push_back(atm);
        }
        clipper::Atom_list atoms(atomvec);

        // Calculate structure factors
        clipper::SFcalc_iso_fft<float> sfc;
        sfc(fphi, atoms);

        cout << "# time: " << time << endl;
        cout << "loop_" << endl
                << "_refln.index_h" << endl
                << "_refln.index_k" << endl
                << "_refln.index_l" << endl
                << "_refln.F_meas_au" << endl
                << "_refln.phase_meas" << endl;

        for (clipper::HKL_info::HKL_reference_index ih = fphi.first_data(); !ih.last(); fphi.next_data(ih)) {
          cout << setw(6) << ih.hkl().h()
                  << setw(6) << ih.hkl().k()
                  << setw(6) << ih.hkl().l()
                  << setw(15) << fphi[ih].f()
                  << setw(15) << fphi[ih].phi() * 180.0 / M_PI << "\n";
        }

      } // while frames in file
    } // for traj
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
