/**
 * @file make_rdc_spec.cc
 * converts RDCs listed per residue into GROMOS input format
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor make_rdc_spec
 * @section make_rdc_spec converts RDCs listed per residue into GROMOS input forma
 * @author @ref gp, ja
 * @date 27. 5. 2009
 *
 * how to use
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@rdc</td><td>&lt;rdc data file&gt; </td></tr>
 * <tr><td> [\@weights</td><td>&lt;rdc-specific weight factors&gt;] </td></tr>
 * <tr><td> [\@nmf</td><td>&lt;number of magnetic field vectors (default: 1)&gt;] </td></tr>
 * <tr><td> \@type</td><td>&lt;type of rdc (as in library file)&gt; </td></tr>
 * <tr><td> \@lib</td><td>&lt;rdc library file&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo         topo.top
     @rdc          NH.rdc
     @type         1
     @lib          rdc.lib
   @endverbatim

 * <hr>
 *
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
#include "../src/gcore/LJException.h"
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

struct rdc_struct {
  int type;
  string name_i;
  string name_j;
  string name_k;
  string name_l;
  double gyr_i;
  double gyr_j;
  double rij;
  double rik;
};

struct rdc_data_struct {
  int residue;
  int num_i;
  int num_j;
  int num_k;
  int num_l;
  double rdc;
  double weight;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "rdc" << "weights" << "nmf" << "type" << "lib";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo        <molecular topology file>\n";
  usage += "\t@rdc         <rdc data file>\n";
  usage += "\t[@weights     <rdc-specific weight factors>]\n";
  usage += "\t[@nmf         <number of magnetic field vectors (default: 1)>]\n";
  usage += "\t@type        <type of rdc>\n";
  usage += "\t@lib         <rdc library file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());
    // check number of molecules
    if (sys.numMolecules() > 1) {
      throw gromos::Exception("make_rdc_spec", "cannot have more than one"
              "molecule in topology: run separately for each molecule");
    }


    // get the type of RDC
    int rdctype;
    if (args.count("type") != 1) {
      throw gromos::Exception("make_rdc_spec", "No or too many rdc type(s) given");
    } else {
      istringstream is(args["type"]);
      is >> rdctype;
    }

    // get information for this rdc type from the library
    if (args.count("lib") != 1)
      throw gromos::Exception("make_rdc_spec", "No rdc library file");
    Ginstream lib_file(args["lib"]);
    vector<string> buffer;
    lib_file.getblock(buffer);
    if (buffer[0] != "RDCSPEC")
      throw gromos::Exception("make_rdc_spec",
            "RDC library file does not contain a RDCSPEC block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("make_rdc_spec", "RDC library file " + lib_file.name() +
            " is corrupted. No END in RDCSPEC"
            " block. Got\n"
            + buffer[buffer.size() - 1]);

    rdc_struct rdc_spec;
    bool found_rdc_type = false;
    vector<string>::const_iterator il = buffer.begin() + 1, lo = buffer.end() - 1;
    for (; il != lo; ++il) {
      istringstream line(*il);
      line >> rdc_spec.type;
      if (rdc_spec.type == rdctype) {
        found_rdc_type = true;
        line >> rdc_spec.name_i >> rdc_spec.name_j >>
                rdc_spec.name_k >> rdc_spec.name_l >> rdc_spec.gyr_i >>
                rdc_spec.gyr_j >> rdc_spec.rij >> rdc_spec.rik;
        if (line.fail())
          throw gromos::Exception("make_rdc_spec", "bad line in RDCSPEC block!");
        break;
      }
    } // RDCSPEC block
    if (!found_rdc_type) {
      throw gromos::Exception("make_rdc_spec", "Unknown rdc type given");
    }


    // read in the rdc data
    vector<rdc_data_struct> rdc_data;
    string rdc_filename;
    if (args.count("rdc") != 1) {
      throw gromos::Exception("make_rdc_spec", "No rdc data file");
    } else {
      rdc_filename = args["rdc"].c_str();
      ifstream rdc_file(args["rdc"].c_str());
      string line;
      std::string::size_type iterator;
      istringstream is;
      while (true) {
        std::getline(rdc_file, line);
        if (rdc_file.eof()) break;
        if ((iterator = line.find('#')) != std::string::npos) {
          if (iterator == 0) continue;
          line = line.substr(0, iterator);
        }
        rdc_data_struct data;
        is.clear();
        is.str(line);
        if (!(is >> data.residue >> data.rdc))
          throw gromos::Exception("make_rdc_spec", "Bad line or non-existent RDC data file!");
        // set weights to one
        data.weight = 1.0;
        rdc_data.push_back(data);
      }
    }


    // read the rdc-specific weights, if present
    if (args.count("weights") == 1) {
      ifstream weights_file(args["weights"].c_str());
      string line;
      std::string::size_type iterator;
      istringstream is;
      while (true) {
        std::getline(weights_file, line);
        if (weights_file.eof()) break;
        if ((iterator = line.find('#')) != std::string::npos) {
          if (iterator == 0) continue;
          line = line.substr(0, iterator);
        }
        is.clear();
        is.str(line);
        int resnum;
        double weight;
        if (!(is >> resnum >> weight)) {
          throw gromos::Exception("make_rdc_spec", "bad line in RDC data file!");
        } else {
          // find right place in data struct
          bool found = false;
          for (unsigned int i = 0; i < rdc_data.size(); i++) {
            int residue = rdc_data[i].residue;
            if (residue == resnum) {
              rdc_data[i].weight = weight;
              found = true;
            }
          }
          if (!found)
            throw gromos::Exception("make_rdc_spec", "residues with weights "
                  "do not match residues with RDCs");
        }
      }
    }

    // get the number of magnetic field vectors
    int nmf = 1;
    if (args.count("nmf") == 1) {
      istringstream is(args["nmf"]);
      is >> nmf;
    }

    // initialise the output
    ostringstream out;
    out << "TITLE\n"
            << "RDC specifications created from " << rdc_filename << " using make_rdc_spec\n"
            << "END\n";
    out << "CONVERSION\n"
            << "# factors to convert the frequency from [RDC]=(s-1) to (ps-1)\n"
            << "# and to convert the gyromagnetic ratios from [gamma]=10^7*(rad/T s) to (e/u)\n"
            << "0.000000000001\n"
            << "0.10375\n"
            << "END\n"
            << "MAGFIELD\n"
            << "# NGF: number of magnetic field vectors\n"
            << "# x/y/z coordinates and mass of the atom representing each magnetic field vector\n"
            << "# NGF\n";
    if (nmf > 1) {
      out << "    " << nmf << "\n"
              << "#      x      y      z    mass\n";
      for (unsigned int i = 0; i < nmf; i++) {
        out << "     0.0  0.001    1.0     1.0\n";
      }
    } else {
      out << "    1\n"
              << "#      x      y      z    mass\n"
              << "     0.0  0.001    1.0     1.0\n";
    }
    out << "END\n"
            << "RDCRESSPEC\n"
            << "# For each RDC restraint the following should be specified:\n"
            << "# IPRDCR, JPRDCR, KPRDCR, LPRDCR, atom numbers (defining the vector that forms the angle with the magnetic field)\n"
            << "# WRDCR                           weight factor (for weighting some RDCs higher than others)\n"
            << "# PRDCR0                          RDC restraint value (i.e. the experimental data)\n"
            << "# RDCGI, RDCGJ                    gyromagnetic ratios for atoms i and j\n"
            << "# RDCRIJ, RDCRIK                  distance between atoms i and j or i and k (RIJ = RCH for CA:HA)\n"
            << "# TYPE                            code to define type of RDC\n"
            << "# IPRDCR JPRDCR KPRDCR LPRDCR   WRDCR    PRDCR0      RDCGI      RDCGJ     RDCRIJ     RDCRIK    RDCTYPE\n";


    // loop over all atoms in molecule
    for (unsigned int i = 0; i < sys.mol(0).topology().numAtoms(); ++i) {

      // for C:N RDCs, the N is from the next residue
      int thisRes = sys.mol(0).topology().resNum(i) + 1;
      int prevRes = 1;
      if (thisRes >= 1)
        prevRes = thisRes - 1;

      // loop over RDC data
      for (unsigned int j = 0; j < rdc_data.size(); ++j) {

        // N:H or CA:C
        switch(rdctype){
	case 1:
        case 2:
        // if residue matches
        {
          if (thisRes == rdc_data[j].residue) {
            // check atom name
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i+1;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i+1;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
        }
          break;
          // C:N
          case 3:
          // find the C
          {
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i+1;
            }
          } else if (prevRes == rdc_data[j].residue) {
            // find the N (from residue i+1)
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i+1;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
          }
          break;
          // CA:HA - need to get CA, N, C and CB (all from same residue)
          case 4:
          {
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i+1;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i+1;
            } else if (atom_name == rdc_spec.name_k) {
              rdc_data[j].num_k = i+1;
            } else if (atom_name == rdc_spec.name_l) {
              rdc_data[j].num_l = i+1;
            }
          }
          }
          break;
          // side-chain N:H (two possible hydrogens)
          case 5:
          case 6:
          {
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i+1;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i+1;
            } else if (atom_name == rdc_spec.name_k) {
              rdc_data[j].num_k = i+1;
            }
          }
          rdc_data[j].num_l = 0;
          }
          break;
          // side-chain H:H
          case 7:
          case 8:
          {
          if (thisRes == rdc_data[j].residue) {
            string atom_name = sys.mol(0).topology().atom(i).name();
            if (atom_name == rdc_spec.name_i) {
              rdc_data[j].num_i = i+1;
            } else if (atom_name == rdc_spec.name_j) {
              rdc_data[j].num_j = i+1;
            }
          }
          rdc_data[j].num_k = 0;
          rdc_data[j].num_l = 0;
          }
          break;
        } // switch
      } // RDC data
    } // atoms

    // loop over RDC data again to write output
    for (unsigned int i = 0; i < rdc_data.size(); ++i) {

      out << setw(8) << rdc_data[i].num_i << setw(7) << rdc_data[i].num_j <<
              setw(7) << rdc_data[i].num_k << setw(7) << rdc_data[i].num_l <<
              setw(8) << rdc_data[i].weight << setw(10) << rdc_data[i].rdc <<
              setw(11) << rdc_spec.gyr_i << setw(11) << rdc_spec.gyr_j <<
              setw(11) << rdc_spec.rij << setw(11) << rdc_spec.rik <<
              setw(11) << rdc_spec.type << endl;
    }

    out << "END\n";
    cout << out.str();

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
