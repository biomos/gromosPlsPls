/**
 * @file make_rdc_spec.cc
 * add gyromagnetic ratios to rdc specification input file
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor make_rdc_spec
 * @section make_rdc_spec adds gyromagnetic ratios to rdc specification input file
 * @author @ref gp, ja
 * @date 27. 5. 2009
 *
 * how to use
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@rdc</td><td>&lt;rdc specification input file&gt; </td></tr>
 * <tr><td> \@gyrolib</td><td>&lt;gyromagnetic ratios library file&gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
   atominfo
     @topo         topo.top
     @rdcspec      rdc.rdc
     @gyrolib      gyr.lib
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
  int atomnri;
  int atomnrj;
  int atomnrk;
  double wrdc;
  double rdcval;
  int IACi;
  int IACj;
  int IACk;
  string namei;
  string namej;
  string namek;
  double gyri;
  double gyrj;
  double rij;
  double rik;
  int type;
};

struct gyr_struct {
  int nr;
  string name;
  double gyr;
};

int main(int argc, char **argv) {

  Argument_List knowns;
  knowns << "topo" << "rdcspec" << "gyrolib";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo      <molecular topology file>\n";
  usage += "\t@rdcspec     <rdc specification input file>\n";
  usage += "\t@gyrolib     <gyromagnetic ratios library file>\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    gcore::System sys(it.system());
    LinearTopology topo(sys);

    if (args.count("rdcspec") != 1)
      throw gromos::Exception("rdc_spec", "No rdc specification input file");

    if (args.count("gyrolib") != 1)
      throw gromos::Exception("rdc_spec", "No gyromagnetic ratios library file");

    {
      ////////////////////////////////////////
      // Read in gyromagnetic ratio library //
      ////////////////////////////////////////

      Ginstream gyrospec_file(args["gyrolib"]);
      vector<string> gbuffer;
      gyrospec_file.getblock(gbuffer);
      if (gbuffer[0] != "GYROSPEC")
        throw gromos::Exception("gyro_spec",
              "Gyromagnetic values library file does not contain a GYROSPEC block!");
      if (gbuffer[gbuffer.size() - 1].find("END") != 0)
        throw gromos::Exception("gyro_spec", "Gyromagnetic values library file " + gyrospec_file.name() +
              " is corrupted. No END in GYROSPEC"
              " block. Got\n"
              + gbuffer[gbuffer.size() - 1]);

      vector<string>::const_iterator itg = gbuffer.begin() + 1, tog = gbuffer.end() - 1;
      map<int, double> gyromap;

      for (; itg != tog; ++itg) {
        gyr_struct gyros;
        istringstream lineg(*itg);
        lineg >> gyros.nr >> gyros.name >> gyros.gyr;
        if (lineg.fail())
          throw gromos::Exception("gyro_spec", "bad line in GYROSPEC block!");

        gyromap[gyros.nr] = gyros.gyr;
      } // GYRSPEC block


      //////////////////////////////////////////////////////////////////
      // look up IAC, look up gyromagnetic ratio, write rdc spec file //
      //////////////////////////////////////////////////////////////////
      int extendedformat = 0;

      Ginstream rdcspec_file(args["rdcspec"]);
      ostringstream out;
      vector<string> rbuffer;

      out << "TITLE\n"
              << rdcspec_file.title()
              << "END\n\n";

      rdcspec_file.getblock(rbuffer);
      if (rbuffer[0] == "CONVERSION") {
        extendedformat = 1;
        out << "CONVERSION\n"
                << "# factors\n"
                << "# to convert the frequency from [RDC]=(Hz) to (ps-1)\n"
                << "# and to convert gyromagnetic ratios from [gamma]=10^7*(rad/T s) to (e/u)\n";

        vector<string>::const_iterator itc = rbuffer.begin() + 1, toc = rbuffer.end() - 1;
        for (; itc != toc; ++itc) {
          stringstream line(*itc);
          out << line.str() << endl;
        }

        out << "END\n\n";

        rdcspec_file.getblock(rbuffer);
        out << "MAGFIELD\n"
                << "# The x/y/z coordinates and the mass of the atom representing the magnetic field vector have to be specified:\n"
                << "#x   y   z      mass\n";
        vector<string>::const_iterator itm = rbuffer.begin() + 1, tom = rbuffer.end() - 1;
        for (; itm != tom; ++itm) {
          stringstream line(*itm);
          out << line.str() << endl;
        }
        out << "END\n\n";

        rdcspec_file.getblock(rbuffer);
      }

      if (rbuffer[0] != "RDCVALRESSPEC")
        throw gromos::Exception("rdc_spec",
              "RDC values specification file does not contain a RDCVALRESSPEC block!");
      if (rbuffer[rbuffer.size() - 1].find("END") != 0)
        throw gromos::Exception("rdc_spec", "RDC values specification file " + rdcspec_file.name() +
              " is corrupted. No END in RDCVALRESSPEC"
              " block. Got\n"
              + rbuffer[rbuffer.size() - 1]);

      if (!extendedformat) {
        out << "CONVERSION\n"
                << "# factors\n"
                << "# to convert the frequency from [RDC]=(s-1) to (ps-1)\n"
                << "# and to convert gyromagnetic ratios from [gamma]=10^7*(rad/T s) to (e/u)\n"
                << "0.000000000001\n"
                << "0.10375\n"
                << "END\n\n"
                << "MAGFIELD\n"
                << "# The x/y/z coordinates and the mass of the atom representing the magnetic field vector have to be specified:\n"
                << "#x   y   z      mass\n"
                << "0.0 0.001 1.0     1\n"
                << "END\n\n";
      }

      out << "RDCVALRESSPEC\n"
              << "# For each RDC restraint the following should be specified:\n"
              << "# IPRDCR, JPRDCR, KPRDCR atom numbers (defining the vector that forms the angle with the magnetic field)\n"
              << "# WRDCR                  weight factor (for weighting some RDCs higher than others)\n"
              << "# PRDCR0                 RDC restraint value (i.e. the experimental data)\n"
              << "#\n"
              << "# IPRDCR JPRDCR KPRDCR   WRDCR    PRDCR0      RDCGI      RDCGJ     RDCRIJ     RDCRIK    RDCTYPE\n\n";

      vector<string>::const_iterator it = rbuffer.begin() + 1, to = rbuffer.end() - 1;
      for (; it != to; ++it) {
        rdc_struct param;
        stringstream line(*it);
        line >> param.atomnri >> param.atomnrj >> param.atomnrk >> param.wrdc
                >> param.rdcval >> param.rij >> param.type;
        if (line.fail())
          throw gromos::Exception("rdc_spec", "bad line in RDCVALRESSPEC block!");
        line.clear();

        //get IAC for param.atomnri -> param.IACi
        param.IACi = topo.atoms().at(param.atomnri - 1).iac() + 1;
        param.namei = topo.atoms().at(param.atomnri - 1).name();
        param.IACj = topo.atoms().at(param.atomnrj - 1).iac() + 1;
        param.namej = topo.atoms().at(param.atomnrj - 1).name();

        // if rij is zero, assign them here (otherwise the values read from file are kept)
        if (param.rij == 0.0) {

          // first for "normal" RDCs
          if (param.type == 0) {

            // rik is redundant
            param.rik = 0.0;

            // N:H or H:N
            if ((param.namei == "H" && param.namej == "N") ||
                    (param.namej == "H" && param.namei == "N")) {
              //rij = 1.04e-10;
              param.rij = 0.104;
              // C:H(n) or H(n):C
            } else if ((param.namei == "C" && param.namej == "H") ||
                    (param.namej == "C" && param.namei == "H")) {
              //rij = 2.04e-10;
              param.rij = 0.204;
              // N:C or C:N
            } else if ((param.namei == "N" && param.namej == "C") ||
                    (param.namej == "N" && param.namei == "C")) {
              //rij = 1.33e-10;
              param.rij = 0.133;
              // CA:C or C:CA
            } else if ((param.namei == "CA" && param.namej == "C") ||
                    (param.namej == "CA" && param.namei == "C")) {
              //rij = 1.53e-10;
              param.rij = 0.153;
            } else {
              out << "Internuclear distance not known for atoms " << param.atomnri << " and " <<
                      param.atomnrj;
              throw gromos::Exception("make_rdc_spec", out.str());
            }
          }// check that if k != 0/type >0 then atoms j and k are of the same type
          else if (param.type > 0) {
            if (param.atomnrk <= 0) {
              throw gromos::Exception("make_rdc_spec",
                      "For RDCTYPE > 0 atom k must be >= 0");
            } else {
              param.IACk = topo.atoms().at(param.atomnrk - 1).iac() + 1;
              param.namek = topo.atoms().at(param.atomnrk - 1).name();
              if (param.IACj != param.IACk) {
                throw gromos::Exception("make_rdc_spec",
                        "For RDCTYPE > 0 atoms j and k must be of the same type");
              } else if (!topo.atoms().at(param.atomnrj - 1).isH() ||
                      !topo.atoms().at(param.atomnrk - 1).isH()) {
                throw gromos::Exception("make_rdc_spec",
                        "For RDCTYPE > 0 atoms j and k must be hydrogens");
              } else {
                // assign internuclear distances (if not given)
                if (param.rij != 0.0) {
                  param.rij = 0.104;
                }
                // and set rik = rij
                param.rik = param.rij;
              }
            }
            // for HH
          } else if (param.type < 0) {
            if (!topo.atoms().at(param.atomnri - 1).isH() ||
                    !topo.atoms().at(param.atomnrj - 1).isH()) {
              throw gromos::Exception("make_rdc_spec",
                      "For RDCTYPE < 0 atoms j and k must be hydrogens");
            } else {
              param.rij = 0.0;
              param.rik = 0.0;
            }
          } else {
            throw gromos::Exception("make_rdc_spec",
                    "RDCTYPE not known");
          }
        } else {
          // set rik
          if (param.type > 0) {
            param.rik = param.rij;
          } else {
            param.rik = 0.0;
          }
        }

        //get gyr for param.IACi -> param.gyri
        map<int, double>::const_iterator posgyri = gyromap.find(param.IACi);
        map<int, double>::const_iterator posgyrj = gyromap.find(param.IACj);
        if (posgyri == gyromap.end()) {
          ostringstream out;
          out << "No gyromagnetic ratio in library for atom type: " << posgyri->first;
          throw gromos::Exception("make_rdc_spec", out.str());
        }
        if (posgyrj == gyromap.end()) {
          ostringstream out;
          out << "No gyromagnetic ratio in library for atom type: " << posgyrj->first;
          throw gromos::Exception("make_rdc_spec", out.str());
        }
        param.gyri = posgyri->second;
        param.gyrj = posgyrj->second;

        //write output
        out << setw(8) << param.atomnri << setw(7) << param.atomnrj << setw(7) << param.atomnrk
                << setw(8) << param.wrdc << setw(10) << param.rdcval << setw(11) << param.gyri <<
                setw(11) << param.gyrj << setw(11) << param.rij << setw(11) << param.rik <<
                setw(11) << param.type << endl;

      } // RDCVALSSPECPEC block
      out << "END\n";
      cout << out.str();

    } //routine

  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
