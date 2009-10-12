/**
 * @file make_sasa_top.cc
 * add sasa block to molecular topology file
 * required in order to use SASA/VOL implicit solvent model
 */

/**
 * @page programs Program Documentation
 *
 * @anchor make_sasa_top
 * @section make_sasa_top add sasa block to molecular topology file
 * @author @ref kb @ref ja
 * @date 23. 4. 2009
 *
 * Program make_sasa_top adds the atom-specific information required to use the
 * SASA/VOL implicit solvent model to the molecular topology file. It reads in
 * an existing molecular topology file created using @ref make_top, along with a
 * SASA/VOL specification library file, which contains the atom-specific SASA/VOL parameters.
 * The specification library file must be for the same force field as was used to create
 * the molecular topology file.
 * 
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@sasaspec</td><td>&lt;sasa specification library file&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   make_sasa_top
     @topo       ex.top
     @sasaspec   sasaspec45b3.lib
   @endverbatim

 * <hr>
 */

#include <cassert>
#include <locale>
#include <map>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <set>

#include "../src/args/Arguments.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Improper.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gcore/CrossDihedral.h"
#include "../src/gcore/LinearTopology.h"
#include "../src/gio/OutTopology.h"
#include "../src/gio/Ginstream.h"

using namespace gcore;
using namespace gio;
using namespace args;

using namespace std;

struct sasa_parameter {
  double radius;
  double probability;
  double sigma;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "sasa_spec";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo       <molecular topology file>\n";
  usage += "\t@sasaspec   <sasa specification file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());
    LinearTopology topo(sys);

    if (args.count("sasaspec") != 1)
      throw gromos::Exception("sasaspec", "No sasa specification file");

    map<int, sasa_parameter> sasa_spec;
    {
      Ginstream spec_file(args["sasaspec"]);
      vector<string> buffer;
      spec_file.getblock(buffer);
      if (buffer[0] != "SASASPEC")
        throw gromos::Exception("sasaspec",
          "sasa specification file does not contain a SASASPEC block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("sasaspec", "sasa specification file " + spec_file.name() +
          " is corrupted. No END in SASASPEC"
          " block. Got\n"
          + buffer[buffer.size() - 1]);

      vector<string>::const_iterator it = buffer.begin() + 1, to = buffer.end() - 1;

      for (; it != to; ++it) {
        sasa_parameter param;
        istringstream line(*it);
        unsigned int num;
        line >> param.radius >> param.probability >> param.sigma >> num;
        if (line.fail())
          throw gromos::Exception("sasaspec",
            "bad line in SASASPEC block!");

        for (unsigned int i = 0; i < num; ++i) {
          int iac;
          line >> iac;
          if (line.fail()) {
            ostringstream msg;
            msg << "bad line in SASASPEC block: could not read " << num
                << " IACs from line.";
            throw gromos::Exception("sasaspec", msg.str());
          }
          sasa_spec[iac - 1] = param;
        } // for iacs
      } // SASASPEC block
    }

    int totNumAt = 0;
    for (int i = 0; i < sys.numMolecules(); ++i) {
      totNumAt += sys.mol(i).numAtoms();
    }

    // write out original topology
    OutTopology ot(cout);
    ot.setTitle(it.title());
    ot.write(sys, it.forceField());

    // write out sasa block
    cout << "SASAPARAMETERS" << endl;
    std::vector<AtomTopology>::const_iterator atom = topo.atoms().begin(),
        atom_to = topo.atoms().end();
    cout << "# number of atoms\n";
    cout << totNumAt << endl;
    cout << "#atom      radius         prob          sigma" << endl;
    
    for (int i = 1; atom != atom_to; ++atom, ++i) {
      map<int, sasa_parameter>::const_iterator result = sasa_spec.find(atom->iac());
      if (result == sasa_spec.end()) {
        ostringstream out;
        out << "No SASA parameters for atom type: " << atom->iac() + 1;
        throw gromos::Exception("sasaspec", out.str());
      }

      const sasa_parameter & s = result->second;
      cout.precision(8);
      cout << setw(5) << i;
      cout << setw(15) << s.radius << setw(15) << s.probability << setw(15) << s.sigma << endl;
    }
    cout << "END" << endl;
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
