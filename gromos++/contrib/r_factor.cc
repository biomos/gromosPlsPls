/**
 * @file r_factor.cc
 * calculates r factors
 */
/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor r_factor
 * @section r_factor calculates crystallographic R factors
 * @author @ref ns ff
 * @date 8.4.2009
 *
 * Program r factor calculates X-ray reflection structure factor amplitudes
 * and phases from a given trajectory and compares them to experimental values.
 * Only the atoms given by the @ref AtomSpecifier
 * \@atomssf are considered for the calculation. The atoms' IAC are mapped to their
 * element names according to the rules given in the \@map file. The atoms' B-factors
 * and occupancies are read from a special file (\@bfactor) if requested or default
 * to @f$ 0.01 \mathrm{nm}^2 @f$ and 100%.
 * Structure factors are calculated to the given resolution (\@resultion) while
 * the cell information is calculated from the system's box.
 * Symmetry operations are taken into account by specifing a (\@spacegroup).
 * Make sure you only give asymetric unit when using \@spacegroup.
 * The program can write the electron density to special files (FRAME_DENSITY_00001.ccp4), if
 * requested (\@density flag).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time and dt"&gt; </td></tr>
 * <tr><td> \@atomssf</td><td>&lt;@ref AtomSpecifier: atoms to consider for structure_factor&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;file with IAC-to-elementname mapping&gt; </td></tr>
 * <tr><td> \@bfactor</td><td>&lt;file with experimental B-factors and occupancies&gt; </td></tr>
 * <tr><td> \@resolution</td><td>&lt;scattering resolution [nm]&gt; </td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup in Hermann-Mauguin format, default: P 1&gt;]</td></tr>
 * <tr><td> \@cif</td><td>&lt;cristallographic information file&gt; </td></tr>
 * <tr><td>[\@density</td><td>&lt;write density to files&gt;]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;convert length unit to Angstrom&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
 r_factor
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
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/Time.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Physics.h"
#include "../src/gio/InIACElementNameMapping.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/gio/InCIF.h"
#include "../src/bound/Boundary.h"
#include "../src/utils/debug.h"

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
using namespace bound;

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "topo" << "pbc" << "traj" << "map" << "atomssf" << "time" << "bfactor"
          << "resolution" << "spacegroup" << "cif" << "density" << "factor";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@pbc         <boundary type> [<gathermethod>]\n";
  usage += "\t[@time       <time and dt>]\n";
  usage += "\t@atomssf     <atomspecifier: atoms to consider for structure_factor>\n";
  usage += "\t@traj        <trajectory files>\n";
  usage += "\t@map         <IAC-to-ElementName map-file>\n";
  usage += "\t[@bfactor    <experimental B-factors>]\n";
  usage += "\t@resolution  <scattering resolution>\n";
  usage += "\t[@spacegroup <spacegroup in Hermann-Maugin format, default: P 1>]\n";
  usage += "\t@cif         <cristallographic information file>\n";
  usage += "\t[@density    <write electron density to special files>]\n";
  usage += "\t[@factor     <convert length unit to Angstrom. default: 10.0>]\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.precision(8);

  try {
    Arguments args(argc, argv, knowns, usage);
    double factor = 10.0;
    {
      Arguments::const_iterator iter = args.lower_bound("factor");
      Arguments::const_iterator to = args.upper_bound("factor");
      if (iter != to) {
        if (!(istringstream(iter->second) >> factor))
          throw gromos::Exception(argv[0], "factor is not numeric.");
        ++iter;
      }
      if (iter != to)
        throw gromos::Exception(argv[0], "factor takes one argument only.");
    }

    // Hardcoded B-factor conversion factor.
    const double sqpi2=(M_PI*M_PI*8.0);

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
    Boundary *pbc = BoundaryParser::boundary(sys, args);
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
    InIACElementNameMapping mapfile(args["map"]);
    std::map<int, std::string> gacmapping = mapfile.getData();

    // Get experimental Bfactors and occupancy
    std::vector<BFactorOccupancyData> bfoc;
    bool has_bfactor = false;
    if (args.count("bfactor") == 1) {
      InBFactorOccupancy bfac_file(args["bfactor"]);
      bfoc = bfac_file.getData();
      has_bfactor = true;
    }

    // Read in cif-file for generating the structure factor list (using clipper)
    InCIF ciffile(args["cif"]);
    vector<CIFData> cifdata = ciffile.getData();

    bool write_density = false;
    if (args.count("density") >= 0) {
      write_density = true;
    }

    //===========================
    // loop over all trajectories
    InG96 ic;

    cout << "#" << setw(14) << "time" << setw(15) << "R" << setw(15) << "k" << endl;
    unsigned int frame_number = 0;
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {
      ic.open(iter->second);
      ic.select("ALL");

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys >> time;
        ++frame_number;
        if (!sys.hasPos)
          throw gromos::Exception(argv[0],
                "Unable to read POSITION(RED) block from "
                "trajectory file.");
        
        if (!sys.hasBox)
          throw gromos::Exception(argv[0],
                "Cannot calculate structure factors without a box.");

        // get the centre of the box
        const Vec centre = (sys.box().K() * 0.5) + (sys.box().L() * 0.5) +
                (sys.box().M() * 0.5);

        // put atom into positive box
        for(int i = 0; i < calcatoms.size(); ++i) {
          calcatoms.pos(i) = calcatoms.pos(i) - pbc->nearestImage(calcatoms.pos(i), centre, sys.box()) + centre;
        }

        // create the cell
        const double a = sys.box().K().abs();
        const double b = sys.box().L().abs();
        const double c = sys.box().M().abs();
        DEBUG(5, "Cell: " << a << " " << b << " " << c)
        const double alpha = sys.box().alpha();
        const double beta = sys.box().beta();
        const double gamma = sys.box().gamma();


        if (!a || !b || !c || !alpha || !beta || !gamma)
          throw gromos::Exception(argv[0], "Box has zero volume!");

        clipper::Cell_descr cellinit(a * factor, b * factor, c * factor,
                alpha, beta, gamma);
        clipper::CCell cell(spgr, clipper::String("base cell"), clipper::Cell(cellinit));

        // create the resolutions and corresponding lattice
        clipper::CResolution reso(cell, clipper::String("base reso"), clipper::Resolution(resolution * factor));
        clipper::CHKL_info hkls(reso, clipper::String("base hkls"), true);
        clipper::CHKL_data<clipper::data64::F_phi> fphi(hkls);
        clipper::CHKL_data<clipper::data64::F_phi> fphi_print(hkls);

        // Fill Clipper Atom list
        // we do this insight the loop due to solvent molecules!
        std::vector<clipper::Atom> atomvec;
        for (int i = 0; i < calcatoms.size(); i++) {
          clipper::Atom atm;
          // convert to angstrom
          atm.set_coord_orth(clipper::Coord_orth(
                  calcatoms.pos(i)[0] * factor,
                  calcatoms.pos(i)[1] * factor,
                  calcatoms.pos(i)[2] * factor));
          
          if (has_bfactor) {
            const unsigned int atom_index = calcatoms.gromosAtom(i);
            if (atom_index >= bfoc.size()) {
              throw gromos::Exception("structre_factor", "Not enough B-factors given");
            }
            atm.set_occupancy(bfoc[atom_index].occupancy);
            // convert to Angstrom^2
            atm.set_u_iso(bfoc[atom_index].b_factor * factor * factor / sqpi2);
          } else {
            atm.set_occupancy(1.0);
            atm.set_u_iso(1.0 / sqpi2);
          }
          atm.set_element(gacmapping[calcatoms.iac(i)]);
          if (atm.occupancy() != 0.0) {
            atomvec.push_back(atm);
          }
        }
        clipper::Atom_list atoms(atomvec);

        // Calculate structure factors
        clipper::SFcalc_iso_fft<double> sfc;
        sfc(fphi, atoms);

        if (write_density)
          fphi_print = std::complex<double>(0.0,0.0);

        // calculate the scaling constant
        double sum_obs_calc = 0.0;
        double sum_calc_calc = 0.0;
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          const double f_obs = cifdata[i].f_obs;
          const double f_calc = fphi[clipper::HKL(cifdata[i].h, cifdata[i].k, cifdata[i].l)].f();
          sum_obs_calc += f_obs * f_calc;
          sum_calc_calc += f_calc * f_calc;
        }
        const double k = sum_obs_calc / sum_calc_calc;
        // and calculate R
        double sum_dev_f = 0.0;
        double sum_obs = 0.0;
        for(unsigned int i = 0; i < cifdata.size(); ++i) {
          const double f_obs = cifdata[i].f_obs;
          const clipper::HKL hkl(cifdata[i].h, cifdata[i].k, cifdata[i].l);
          const double f_calc = fphi[hkl].f();
          sum_dev_f += fabs(f_obs - k*f_calc);
          sum_obs += f_obs;
          if (write_density)
            fphi_print.set_data(hkl, clipper::data64::F_phi(f_obs / k, fphi[hkl].phi()));
        }

        const double R = sum_dev_f / sum_obs;
        cout << time << setw(15) << R << setw(15) << k << endl;

        // write the electron density
        if (write_density) {
          {
            ostringstream file_name;
            file_name << "FRAME_DENSITY_" << setfill('0') << setw(6) << frame_number << ".ccp4";
            const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
            clipper::Xmap<clipper::ftype64> density(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);
            density.fft_from(fphi_print);

            for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = density.first();
                    !ix.last(); ix.next())
              density[ix] /= density.multiplicity(ix.coord());

            clipper::CCP4MAPfile mapfile;

            mapfile.open_write(file_name.str());
            mapfile.export_xmap(density);
            mapfile.close_write();
          }
          {
            ostringstream file_name;
            file_name << "FRAME_CALCDENSITY_" << setfill('0') << setw(6) << frame_number << ".ccp4";
            const clipper::Grid_sampling grid(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), fphi.base_hkl_info().resolution(), 1.5);
            clipper::Xmap<clipper::ftype64> density(fphi.base_hkl_info().spacegroup(), fphi.base_cell(), grid);
            density.fft_from(fphi);

            for (clipper::Xmap<clipper::ftype32>::Map_reference_index ix = density.first();
                    !ix.last(); ix.next())
              density[ix] /= density.multiplicity(ix.coord());

            clipper::CCP4MAPfile mapfile;

            mapfile.open_write(file_name.str());
            mapfile.export_xmap(density);
            mapfile.close_write();
          }
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
