/**
 * @file rmsdmat.cc
 * Calculates the rmsd-matrix for given structures
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rmsdmat
 * @section rmsdmat Calculates the rmsd-matrix for given structures
 * @author @ref mc @ref co
 * @date 22-8-06
 *
 * Program rmsdmat calculates the root-mean-square deviation
 * between all pairs of structures in a given trajectory file. This matrix of
 * RMSD's can subsequently be used by program @ref cluster to perform a
 * conformational clustering. The matrix can be written out in human readable
 * form, or -to save disk space- in binary format. For efficiency reasons, the
 * RMSD values are written in an integer format. The user can specify the
 * required precision of the RMSD values that are stored. In the case of the
 * binary format, the values are stored as unsigned short int if the precision is
 * less or equal to 4, or otherwise as unsigned int.
 *
 * For an atom-positional RMSD matrix different sets of atoms can be selected to perform a rotational
 * least-squares-fit and to calculate the RMS deviation from. The RMSD matrix
 * can also be calculated from differences in internal coordinates defined by
 * a set of properties (e.g. torsional angles or hydrogen bonds).
 * A selection of
 * structures in the trajectory file to consider can be made using the options
 * skip and stride. Structure pairs may occur for which the least-squares
 * rotational fit fails for numerical reasons. In these cases both structures
 * are fit to the reference structure. If no user specified reference structure
 * is available, the first structure in the trajectory is taken as such.
 * Specifying a reference structure allows the program @ref cluster to perform a
 * forced clustering as well, requiring that the first cluster contains the
 * reference structure, regardless of the cluster size.
 *
 * <B>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary conditions&gt; &lt;gather type&gt; </td></tr>
 * <tr><td> [\@atomsfit</td><td>&lt;@ref AtomSpecifier "atoms" to consider for fit&gt;]</td></tr>
 * <tr><td> [\@atomsrmsd</td><td>&lt;@ref AtomSpecifier "atoms" to consider for rmsd&gt; </td></tr>
 * <tr><td> [\@prop</td><td>&lt;@ref PropertSpecifier "properties" to be used for rmsd computation.&gt;]</td></tr>
 * <tr><td> [\@skip</td><td>&lt;skip frames at beginning&gt;] </td></tr>
 * <tr><td> [\@stride</td><td>&lt;use only every step frame&gt;] </td></tr>
 * <tr><td> [\@human</td><td>(write the matrix in human readable form)] </td></tr>
 * <tr><td> [\@precision</td><td>&lt;number of digits in the matrix (default 4)&gt;] </td></tr>
 * <tr><td> [\@ref</td><td>&lt;reference coordinates&gt;] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rmsdmat
    @topo         ex.top
    @pbc          r
    @atomsfit     1:a
    @atomsrmsd    1:CA
    @skip         5
    @stride       10
    @human
    @precision    4
    @ref          exref.coo
    @traj         ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/System.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/PropertyContainer.h"
#include "../src/utils/Property.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/FastRotationalFit.h"

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace bound;
using namespace gio;
using namespace utils;
using namespace args;
using namespace fit;

double props_rmsd(const PropertyContainer & props, int i, int j) {
  double rmsd = 0.0;
  for (PropertyContainer::const_iterator it = props.begin(), to = props.end();
          it != to; ++it) {
    rmsd += (*it)->nearestImageDistance((*it)->getValue(i), (*it)->getValue(j)).abs2();
  }
  return std::sqrt(rmsd / props.size());
}

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "traj" << "pbc" << "ref" << "atomsrmsd" << "atomsfit"
            << "skip" << "stride" << "human" << "precision" << "prop";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@pbc          <boundary conditions> <gather type>\n";
  usage += "\t[@prop        <PropertySpecifier>]\n";
  usage += "\t[@atomsfit    <atoms to consider for fit>]\n";
  usage += "\t[@atomsrmsd   <atoms to consider for rmsd>\n";
  usage += "\t              (only specify if different from atomsfit)]\n";
  usage += "\t[@skip        <skip frames at beginning>]\n";
  usage += "\t[@stride      <use only every step frame>]\n";
  usage += "\t[@human       (write the matrix in human readable form)]\n";
  usage += "\t[@precision   <number of digits in the matrix (default 4)>]\n";
  usage += "\t[@ref         <reference coordinates>]\n";
  usage += "\t@traj         <trajectory files>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // read the fit atoms
    AtomSpecifier fitatoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsfit"),
              to = args.upper_bound("atomsfit");
      for (; iter != to; ++iter) {
        fitatoms.addSpecifier(iter->second);
      }
    }

    // read the rmsd atoms
    AtomSpecifier rmsdatoms(sys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsrmsd"),
              to = args.upper_bound("atomsrmsd");
      for (; iter != to; ++iter) {
        rmsdatoms.addSpecifier(iter->second);
      }
      if (rmsdatoms.size() == 0)
        rmsdatoms = fitatoms;

    }

    // for which atoms do we want to keep the coordinates
    AtomSpecifier atoms(fitatoms + rmsdatoms);

    // if fitatoms != rmsdatoms keep lists of what to do
    vector<bool> fit_spec, rmsd_spec;
    if (atoms.size() != fitatoms.size() || atoms.size() != rmsdatoms.size()) {
      fit_spec.resize(atoms.size(), true);
      rmsd_spec.resize(atoms.size(), true);
      for (int i = 0; i < atoms.size(); ++i) {
        if (fitatoms.findAtom(atoms.mol(i), atoms.atom(i)) < 0)
          fit_spec[i] = false;
        if (rmsdatoms.findAtom(atoms.mol(i), atoms.atom(i)) < 0)
          rmsd_spec[i] = false;
      }
    }

    FastRotationalFit frf(fit_spec, rmsd_spec);

    int skip = args.getValue<int>("skip", false, 0);
    int stride = args.getValue<int>("stride", false, 1);

    // read the precision
    int ii = args.getValue<int>("precision", false, 4);
    int precision = 1;
    for (int i = 0; i < ii; ++i) {
      precision *= 10;
    }

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, refSys, args);

    // read in properties
    PropertyContainer props(sys, pbc);
    {
      Arguments::const_iterator iter = args.lower_bound("prop");
      Arguments::const_iterator to = args.upper_bound("prop");
      for (int i = 0; iter != to; iter++, ++i) {
        props.addSpecifier(iter->second);
      }
    }

    if ((fitatoms.empty() && props.empty()) ||
            (!fitatoms.empty() && !props.empty()))
      throw gromos::Exception(argv[0],
            "Give either fit atoms or properties.");

    // create the vector to store the trajectory
    vector< vector < Vec > > traj;
    vector< Vec > frame(atoms.size());

    // read reference coordinates...
    InG96 ic;
    if (args.count("ref") > 0) {
      ic.open(args["ref"]);
    } else {
      ic.open(args.lower_bound("traj")->second);
    }
    ic >> sys;
    ic.close();

    (*pbc.*gathmethod)();

    // calculate the centre of geometry of the relevant atoms
    Vec cog;
    for (int i = 0; i < atoms.size(); ++i) {
      cog += *atoms.coord(i);
    }
    cog /= atoms.size();

    // put it in the trajectory
    for (int i = 0; i < atoms.size(); ++i) {
      frame[i] = *atoms.coord(i) - cog;
    }
    traj.push_back(frame);

    int framenum = 0;

    // loop over all trajectories
    for (Arguments::const_iterator iter = args.lower_bound("traj");
            iter != args.upper_bound("traj"); ++iter) {

      ic.open(iter->second);

      // loop over all frames
      while (!ic.eof()) {
        ic >> sys;

        if (!((framenum - skip) % stride)) {

          //pbc call
          (*pbc.*gathmethod)();

          if (props.size()) {
            props.calc();
          } else {
            Vec cog;
            for (int i = 0; i < atoms.size(); ++i) {
              cog += *atoms.coord(i);
            }
            cog /= atoms.size();
            for (int i = 0; i < atoms.size(); ++i) {
              frame[i] = *atoms.coord(i) - cog;
            }
            int err = frf.fit(traj[0], frame);

            if (err) {
              ostringstream os;
              os << "Error while fitting to the reference structure\n"
                      << "Error code " << err << " in frame number " << framenum + 1;

              throw gromos::Exception(argv[0], os.str());
            }

            // store coordinates from sys in traj
            traj.push_back(frame);
          }
        }
        framenum++;
      }
      ic.close();
    }

    // everything is in the thing; create the thingy
    // open a file
    ofstream fout;
    bool human = false;
    if (args.count("human") >= 0) {
      fout.open("RMSDMAT.dat");
      human = true;
    } else {
      fout.open("RMSDMAT.bin", ios::out | ios::binary);
    }

    // make a double loop
    int num = traj.size();
    if (props.size())
      num = props.front()->num();
    Matrix rot(3, 3, 0);
    Matrix unit(3, 3, 0);
    for (size_t i = 0; i < 3; i++) unit(i, i) = 1;

    cout << "Read " << num << " out of " << framenum
            << " structures from trajectory" << endl;

    if (human) { // human format
      fout << "TITLE\n"
              << "\trmsd-matrix for " << num - 1 << " + 1 (ref) = "
              << num << " structures\n"
              << "END\n"
              << "RMSDMAT\n"
              << "# number of frames   skip   stride\n"
              << num << "\t" << skip << "\t" << stride << "\n"
              << "# precision\n"
              << precision << "\n";

      for (int i = 0; i < num; ++i) {
        for (int j = i + 1; j < num; ++j) {
          double rmsd = 0.0;
          if (props.size()) {
            rmsd = props_rmsd(props, i, j);
          } else {
            if (frf.fit(rot, traj[i], traj[j])) {
              cout << "error rotational fit on frames " << i + 1 << " and "
                      << j + 1
                      << "\nfitting to reference structure instead" << endl;
              rot = unit;
            }

            rmsd = frf.rmsd(rot, traj[i], traj[j]);
          }

          rmsd *= precision;
          if (rmsd > std::numeric_limits<unsigned int>::max())
            throw gromos::Exception(argv[0], "RMSD value is too big for a 'int'. Adjust @precision.");

          fout << setw(8) << i
                  << setw(8) << j
                  << setw(8) << unsigned(rmsd)
                  << endl;
        }
      }
      fout << "END\n";
    } else { // binary format
      fout.write((char*) &num, sizeof (int));
      fout.write((char*) &skip, sizeof (int));
      fout.write((char*) &stride, sizeof (int));
      fout.write((char*) &precision, sizeof (int));

      if (precision < 1e5) { // small precision -> short format
        std::cout << "using 'unsigned short' as format" << std::endl;

        typedef unsigned short ushort;
        ushort irmsd;
        for (int i = 0; i < num; ++i) {
          for (int j = i + 1; j < num; ++j) {
            double rmsd = 0.0;
            if (props.size()) {
              rmsd = props_rmsd(props, i, j);
            } else {
              if (frf.fit(rot, traj[i], traj[j])) {
                cout << "error rotational fit on frames " << i + 1 << " and "
                        << j + 1
                        << "\nfitting to reference structure instead" << endl;
                rot = unit;
              }

              rmsd = frf.rmsd(rot, traj[i], traj[j]);
            }

            rmsd *= precision;
            if (rmsd > std::numeric_limits<unsigned short>::max())
              throw gromos::Exception(argv[0], "RMSD value is too big for a 'short'. Adjust @precision.");

            irmsd = ushort(rmsd);
            fout.write((char*) &irmsd, sizeof (ushort));
          }
        }
      } else { // higher precision -> int format
        std::cout << "using 'unsigned int' as format" << std::endl;
        unsigned irmsd;
        for (int i = 0; i < num; ++i) {
          for (int j = i + 1; j < num; ++j) {
            double rmsd = 0.0;
            if (props.size()) {
              rmsd = props_rmsd(props, i, j);
            } else {
              if (frf.fit(rot, traj[i], traj[j])) {
                cout << "error rotational fit on frames " << i + 1 << " and "
                        << j + 1
                        << "\nfitting to reference structure instead" << endl;
                rot = unit;
              }

              rmsd = frf.rmsd(rot, traj[i], traj[j]);
            }

            rmsd *= precision;
            if (rmsd > std::numeric_limits<unsigned int>::max())
              throw gromos::Exception(argv[0], "RMSD value is too big for a 'int'. Adjust @precision.");

            irmsd = unsigned(rmsd);
            fout.write((char*) &irmsd, sizeof (unsigned int));
          }
        }
      } // if precision
    } // if human
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
