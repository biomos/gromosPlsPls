/**
 * @file cry.cc
 * Perform symmetry operations on molecules
 */

/**
 * @page programs Program Documentation
 *
 * @anchor cry
 * @section cry Perform symmetry operations on molecules
 * @author @ref mc @ref co
 * @date 8-6-07
 *
 * When using periodic boundary conditions, the computational box containing
 * the molecular system is treated as being translationally invariant. So,
 * periodic boundary conditions can also be used when simulating a crystal as
 * long as the unit cell, or a number of adjacent unit cells is used as
 * computational box. Unless the asymmetric unit is translationally invariant,
 * it cannot be used a computational box. Since crystallographic coordinates of
 * moleculars are generally only provided for the molecules in one asymmetric
 * unit, the coordinates of the other molecules in the unit cell or cells are
 * to be generated by crystallographic symmetry transformations.
 *
 * The program cry can rotate and translate copies of a system to create a 
 * crystal unit cell. Based on a topology and initial (gathered) coordinates,
 * as well as a specification file or a spacegroup for the symmetry transformations.
 * In addition to the specification file a conversion factor can be given
 * if the translation vector is specified in different units than the coordinates.
 * The symmetry operation can be generated automatically by specifing a spacegroup
 * given in Hall or in Hermann-Mauguin symbol notation. Further, a cell containing
 * the unit cells edge lengths and angles in degrees have to be given.
 * A coordinate file is generated that
 * contains coordinates for as many systems as there are transformations
 * defined in the specification file. The corresponding topology can easily be
 * generated using the program @ref com_top (section V-2.2).
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file for the molecules&gt; </td></tr>
 * <tr><td>[\@spec</td><td>&lt;specification file for the symmetry transformations]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;conversion factor for distances&gt;]</td></tr>
 * <tr><td>[\@spacegroup</td><td>&lt;spacegroup symbol&gt;]</td></tr>
 * <tr><td>[\@cell</td><td>&lt;cell edge lengths [nm] and angles [degree]&gt;]</td></tr>
 * </table>
 *
 *
 * Example using a specification file:
 * @verbatim
  cry
    @topo    ex.top
    @pos     exref.coo
    @spec    cry.spec
    @factor  0.1
 @endverbatim
 * Example using a spacegroup
 * @verbatim
  cry
    @topo       ex.top
    @pos        exref.coo
    @spacegroup P 21 21 21
    @cell       10 10 10 90 90 90
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <map>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Triclinic.h"

#include "../config.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#endif

using namespace std;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace bound;
using namespace args;

void read_spec(std::string name,
        vector<Matrix> & rotation,
        vector<Vec> &translation,
        double factor);

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "pos" << "spec" << "factor" << "spacegroup" << "cell";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@pos         <input coordinate file for the molecules>\n";
  usage += "\t[@spec       <specification file for the symmetry transformations]\n";
  usage += "\t[@factor     <conversion factor for distances>]\n";
  usage += "\t[@spacegroup <spacegroup symbol, Hall or Hermann-Mauguin>]\n";
  usage += "\t[@cell       <unit cell edge lengths [nm] and angles [degree]>]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    args.check("topo", 1);
    InTopology it(args["topo"]);
    System sys(it.system());

    // prepare the rotation matrices and translation vectors
    vector<Matrix> rotation;
    vector<Vec> translation;

    // create a final system to work on here.
    System finalsys;

    // set this to true of you want to put the molecules into the box
    bool put_into_box = false;

    // this program can either use a specification file OR a spacegroup
    if (args.count("spec") > 0 && args.count("spacegroup") > 0) {
      throw gromos::Exception(argv[0], "Either give a specification file or a spacegroup.");
    }

    // read the specification file
    if (args.count("spec") > 0) {
      // read the conversion factor
      double factor = 1.0;
      if (args.count("factor") > 0) {
        if (!(istringstream(args["factor"]) >> factor))
          throw gromos::Exception(argv[0], "The conversion factor needs to be numeric.");
      }
      // check input for consistency
      if (args.count("cell") > 0)
        throw gromos::Exception(argv[0], "No cell needed when using specification file.");
      read_spec(args["spec"], rotation, translation, factor);
    }

    // create the matrices and vectors from the spacegroup
    if (args.count("spacegroup") > 0) {
#ifdef HAVE_CLIPPER
      // check input for consistency
     if (args.count("factor") > 0)
        throw gromos::Exception(argv[0], "No conversion factor needed when using spacegroup.");
     
     // read the cell
      std::vector<double> cell_data;
      {
        Arguments::const_iterator iter = args.lower_bound("cell");
        Arguments::const_iterator to = args.upper_bound("cell");
        int i = 0;
        for (; iter != to; ++iter, ++i) {
          double temp_cell;
          std::istringstream is(iter->second);
          if (!(is >> temp_cell)) {
            throw gromos::Exception(argv[0],
                    "Cell parameters not numeric");
          }
          cell_data.push_back(temp_cell);
        }
        if (i != 6) {
          throw gromos::Exception(argv[0],
                  "Not enough cell parameters");
        }
      }
      // create the cell objects
      const clipper::Cell_descr cell_descriptor(
              cell_data[0]*10.0, cell_data[1]*10.0, cell_data[2]*10.0,
              cell_data[3], cell_data[4], cell_data[5]);
      const clipper::Cell cell(cell_descriptor);

      // set the box of the final system
      // check whether orthorhombic or monoclinic/triclinic
      Box::boxshape_enum ntb =
              (cell_data[3] == 90.0 && cell_data[4] == 90.0 && cell_data[5] == 90.0) ?
                Box::rectangular : Box::triclinic;

      finalsys.hasBox = true;
      finalsys.box() = Box(ntb, cell_data[0], cell_data[1], cell_data[2],
              cell_data[3], cell_data[4], cell_data[5], 0.0, 0.0, 0.0);
      put_into_box = true;

      // concatenate the spacegroup symbol from the arguments
      std::string spacegroup_symbol;
      {
        Arguments::const_iterator iter = args.lower_bound("spacegroup");
        Arguments::const_iterator to = args.upper_bound("spacegroup");
        spacegroup_symbol = iter->second;
        for (++iter; iter != to; ++iter) {
          spacegroup_symbol += " ";
          spacegroup_symbol += iter->second;
        }
      }
      // get the spacegroup
      clipper::Spacegroup spacegroup;
      try {
        const clipper::Spgr_descr spacegroup_descriptor(spacegroup_symbol);
        spacegroup.init(spacegroup_descriptor);
      } catch (clipper::Message_fatal & msg) {
        throw gromos::Exception(argv[0], msg.text());
      }

      // loop over symmetry operations
      for(int i = 0; i < spacegroup.num_symops(); ++i) {
        // get the rotation/translation operator from the cell and the sym op.
        const clipper::RTop_orth & rt_operator = spacegroup.symop(i).rtop_orth(cell);
        // get the matrix
        const clipper::Mat33<double> & clipper_matrix = rt_operator.rot();
        Matrix rot;
        // convert the matrix to the GROMOS format
        for(unsigned int row = 0; row < 3; ++row) {
          for(unsigned int col = 0; col < 3; ++col) {
            rot(row, col) = clipper_matrix(row, col);
          }
        }
        rotation.push_back(rot);

        // get the translation vector and convert to nm.
        const clipper::Vec3<double> & clipper_vec = rt_operator.trn();
        const Vec trans(clipper_vec[0]/10.0, clipper_vec[1]/10.0, clipper_vec[2]/10.0);
        translation.push_back(trans);
      } // loop over symmetry operations
#else
      throw gromos::Exception(argv[0], "You have to compile GROMOS++ with CCP4"
              " and clipper libraries in order to use the spacegroup feature"
              " of cry.\nUse --with-ccp4 and --with-clipper for configuration.");
#endif

    } // if spacegroup

    
    // tell it already about the solvent that we have
    finalsys.addSolvent(sys.sol(0));

    // read single coordinates
    InG96 ic;
    ic.open(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    ic.close();

    // create one more system to keep the coordinates
    System refsys(sys);

    Triclinic pbc(&finalsys);

    int num = rotation.size();
    for (int i = 0; i < num; ++i) {
      //solute
      for (int m = 0; m < sys.numMolecules(); ++m) {
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a) {
          Vec newPos = rotation[i] * refsys.mol(m).pos(a)
                  + translation[i];
          if (put_into_box) {
            // put the atom back into the unit cell
            newPos = pbc.nearestImage(Vec(0.0, 0.0, 0.0), newPos, finalsys.box());
          }
          sys.mol(m).pos(a) = newPos;
        }
        finalsys.addMolecule(sys.mol(m));
      }
      //solvent
      for (int s = 0; s < sys.sol(0).numPos(); ++s) {
        Vec newPos = rotation[i] * refsys.sol(0).pos(s)
                + translation[i];
        if (put_into_box) {
            // put the atom back into the unit cell
            newPos = pbc.nearestImage(Vec(0.0, 0.0, 0.0), newPos, finalsys.box());
          }
        finalsys.sol(0).addPos(newPos);
      }
    }

    // Print the new set to cout
    OutG96S oc;
    ostringstream title;
    title << "cry applied the following " << num << " symmetry operations: ";
    title.precision(4);
    for(int i = 0; i < num; ++i) {
      title << endl;
      for(unsigned int row = 0; row < 3; ++row) {
        title << "[";
        for(unsigned int col = 0; col < 3; ++col) {
          // looks bad due to numerics
          const double val = rotation[i](row, col) < 1.0e-12 ? 0.0 : rotation[i](row, col);
          title << setw(8) << val;
        }
        title << "] ";
        title << ((row == 1) ? "* ri +" : "      ");
        // dito
        const double val = translation[i][row] < 1.0e-12 ? 0.0 : translation[i][row];
        title << " [" << setw(8) << val << "]" << endl;
      }
    }

    oc.open(cout);
    oc.select("ALL");

    oc.writeTitle(title.str());
    oc << finalsys;
  }  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void read_spec(std::string name,
        vector<Matrix> & rotation,
        vector<Vec> &translation,
        double factor) {
  Ginstream file(name);
  vector<string> buffer;
  file.getblock(buffer);
  file.close();

  if (buffer[0] != "TRANSFORM")
    throw gromos::Exception("cry", "Could not read TRANSFORM block in specification file");
  if (buffer[buffer.size() - 1].find("END") != 0)
    throw gromos::Exception("cry", "Specification file " + file.name() +
          " is corrupted. No END in " + buffer[0] +
          " block. Got\n"
          + buffer[buffer.size() - 1]);
  int num = 0;
  vector<string>::iterator iter = buffer.begin() + 1;

  istringstream is(*iter);
  ++iter;

  if (!(is >> num) || num <= 0)
    throw gromos::Exception("cry", "Need some transformations");
  if (buffer.size() - 3 != unsigned (num * 3))
    throw gromos::Exception("cry", "Line count wrong in " + file.name());

  Matrix rot(3, 3);
  Vec v;

  for (int i = 0; i < num; i++) {
    for (int j = 0; j < 3; j++, ++iter) {
      is.clear();
      is.str(*iter);
      for (int k = 0; k < 3; ++k) {
        if (!(is >> rot(j, k)))
          throw gromos::Exception("cry", "error reading file");
      }
      if (!(is >> v[j]))
        throw gromos::Exception("cry", "error reading file");
    }
    rotation.push_back(rot);
    translation.push_back(factor * v);
  }

}




