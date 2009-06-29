/**
 * @file xray_map.cc
 * Transformations and filtering of crystallographic maps
 */

/**
 * @page contrib Contrib Program Documentation
 *
 * @anchor xray_map
 * @section xray_map Transformations and filtering of crystallographic maps
 * @author @ref ns
 * @date 29-6-2009
 *
 * Program xray_map is used to transform and/or filter crystallographic maps.
 * It reads given CCP4 map files (\@map) and atomic coordinates (\@pos) and
 * writes the final result to a CCP4 map file (\@out).
 *
 * If requested by \@expression an @ref utils::ExpressionParser expression is 
 * evaluated to calculate every grid points value from the maps, which are
 * available in the expression via the <tt>rho1</tt>, <tt>rho2</tt> etc. symbols. 
 * By default the expression <tt>rho1</tt> is evaluated which corresponds of the
 * value of the first given map. A difference map, for example, can be calculated
 * by giving <tt>rho1 - rho2</tt>.
 *
 * The final map can be filtered by a simple cutoff criterion. All grid points
 * closer than a given distance (\@cutoff) to given atom centres (\@centre) are
 * included in the final map. All other grid points are set to zero.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;coordinate file for filtering and expressions&gt; </td></tr>
 * <tr><td> \@map</td><td>&lt;one of many map files&gt; </td></tr>
 * <tr><td>[\@expression</td><td>&lt;@ref utils::ExpressionParser expression applied to the maps&gt;]</td></tr>
 * <tr><td>[\@centre</td><td>&lt;@ref utils::AtomSpecifier "the atoms to select"&gt;]</td></tr>
 * <tr><td>[\@cutoff</td><td>&lt;grid cell cutoff (in nm)&gt;]</td></tr>
 * <tr><td> \@out</td><td>&lt;output file name&gt;</td></tr>
 * </table>
 *
 * Example:
 * @verbatim
  xray_map
    @topo          ex.top
    @pos           ex.cnf
    @map           ex.ccp4 ex2.ccp4
    @expression    rho1 - rho2
    @centre        1:CA,C,O,N
    @cutoff        0.3
    @out           diff.ccp4
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <string>
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
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/InG96.h"
#include "../src/utils/ExpressionParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/bound/Triclinic.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace std;
using namespace utils;
using namespace bound;

// Additional Clipper Headers
#include "../config.h"

#ifdef HAVE_CLIPPER
#include <clipper/clipper.h>
#include <clipper/clipper-ccp4.h>


string map2var(unsigned int i) {
  ostringstream os;
  os << "rho" << i+1;
  return os.str();
}

int main(int argc, char *argv[]) {
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pos            <reference frame>\n";
  usage += "\t@map            <crystallographic map(s)>\n";
  usage += "\t[@expression    <transformation expressions>]\n";
  usage += "\t[@centre        <AtomSpecifier: atoms to select>]\n";
  usage += "\t[@cutoff        <cutoff for selection, default: 0.5>]\n";
  usage += "\t@out            <output file name>";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "pos" << "map" << "expression" << "centre" << "cutoff" << "out";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {
    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology intopo(args["topo"]);
    System sys(intopo.system());

    // read reference frame
    InG96 ic(args["pos"]);
    ic.select("ALL");
    ic >> sys;

    vector<clipper::Xmap<clipper::ftype64> > maps;
    vector<string> map_names;
    {
      Arguments::const_iterator it = args.lower_bound("map"),
              to = args.upper_bound("map");
      for(; it != to; ++it) {
        try {
          clipper::CCP4MAPfile file;
          file.open_read(it->second);
          clipper::Xmap<clipper::ftype64> map;
          file.import_xmap(map);
          file.close_read();
          maps.push_back(map);
          map_names.push_back(map2var(maps.size()-1));
        } catch (const clipper::Message_fatal & msg) {
          throw gromos::Exception(argv[0], msg.text());
        }
      }
      if (maps.empty())
        throw gromos::Exception(argv[0], "no maps given.");

      // check the maps
      for (unsigned int i = 1; i < maps.size(); ++i) {
        // check spacegroup
        if (maps[0].spacegroup().symbol_hm() != maps[i].spacegroup().symbol_hm()) {
          ostringstream msg;
          msg << "Spacegroups of first map ("
                  << maps[0].spacegroup().symbol_hm()
                  << ") and map " << i + 1 << " ("
                  << maps[i].spacegroup().symbol_hm() << ") do not match.";
          throw gromos::Exception(argv[0], msg.str());
        }
        // check cell
        if (maps[0].cell().a() != maps[i].cell().a() ||
            maps[0].cell().b() != maps[i].cell().b() ||
            maps[0].cell().c() != maps[i].cell().c() ||
            maps[0].cell().alpha() != maps[i].cell().alpha() ||
            maps[0].cell().beta() != maps[i].cell().beta() ||
            maps[0].cell().gamma() != maps[i].cell().gamma()) {
          ostringstream msg;
          msg << "Cell dimensions of first map and map " << i + 1 << " do not match.";
          throw gromos::Exception(argv[0], msg.str());
        }
        // check grid
        if (maps[0].grid_sampling().nu() != maps[i].grid_sampling().nu() ||
            maps[0].grid_sampling().nv() != maps[i].grid_sampling().nv() ||
            maps[0].grid_sampling().nw() != maps[i].grid_sampling().nw()) {
          ostringstream msg;
          msg << "Grid dimensions of first map and map " << i + 1 << " do not match.";
          throw gromos::Exception(argv[0], msg.str());
        }
      } // check maps
    } // map argument

    // get the expression
    string expression_string("");
    for (Arguments::const_iterator it = args.lower_bound("expression"),
            to = args.upper_bound("expression"); it != to; ++it) {
      expression_string += it->second + " ";
    }
    if (expression_string.empty())
      expression_string = map_names[0];

    ExpressionParser<double> parser(&sys, NULL);
    std::vector<ExpressionParser<double>::expr_struct> parsed_expression;
    // define known variables (maps)
    std::map<std::string, double> expression_variables;
    for(unsigned int i = 0; i < maps.size(); i++) {
      expression_variables[map_names[i]] = 0.0;
    }

    // parse the expression
    parser.parse_expression(expression_string, expression_variables, parsed_expression);

    // apply expression
    // get resulting map
    clipper::Grid_sampling resgrid(maps[0].grid_sampling().nu(),
            maps[0].grid_sampling().nv(),
            maps[0].grid_sampling().nw());
    clipper::Xmap<clipper::ftype64> resmap(maps[0].spacegroup(), maps[0].cell(), resgrid);
    // loop over maps
    clipper::Xmap<clipper::ftype32>::Map_reference_index ix = resmap.first();
    vector<clipper::Xmap<clipper::ftype64>::Map_reference_index> map_ix;
    for(vector<clipper::Xmap<clipper::ftype64> >::const_iterator it = maps.begin(),
            to = maps.end(); it != to; ++it) {
      map_ix.push_back(it->first());
    }

    // loop over maps
    for (; !ix.last(); ix.next()) {
      // setup known variables for expression
      for(unsigned int i = 0; i < maps.size(); ++i) {
        expression_variables[map_names[i]] = maps[i][map_ix[i]];
        map_ix[i].next();
      }
      resmap[ix] = parser.calculate(parsed_expression, expression_variables);
    }

    // resmap finished, filter?
    AtomSpecifier atoms(sys);
    for (Arguments::const_iterator it = args.lower_bound("centre"),
            to = args.upper_bound("centre"); it != to; ++it) {
      atoms.addSpecifier(it->second);
    }

    if (args.count("centre") > 0 && atoms.empty())
      cerr << "Warning: no atoms selected, disabling filtering." << endl;

    clipper::Grid_sampling outgrid(maps[0].grid_sampling().nu(),
            maps[0].grid_sampling().nv(),
            maps[0].grid_sampling().nw());
    clipper::Xmap<clipper::ftype64> outmap(maps[0].spacegroup(), maps[0].cell(), outgrid);

    if (atoms.empty()) {
      // do not filter
      clipper::Xmap<clipper::ftype32>::Map_reference_index ix = resmap.first();
      clipper::Xmap<clipper::ftype32>::Map_reference_index out_ix = outmap.first();
      for (; !ix.last(); ix.next(), out_ix.next()) {
        outmap[out_ix] = resmap[ix];
      }
    } else {
      // do filtering
      outmap = 0.0;
      // create a box
      sys.box() = Box(Box::triclinic, maps[0].cell().a()/10.0, maps[0].cell().b()/10.0, maps[0].cell().c()/10.0,
              maps[0].cell().alpha_deg(), maps[0].cell().beta_deg(), maps[0].cell().gamma_deg(), 0.0, 0.0, 0.0);
      sys.hasBox = true;
      // create periodic boudnary conditions
      Triclinic pbc(&sys);
      const Vec & centre = (sys.box().K() * 0.5) + (sys.box().L() * 0.5) + (sys.box().M() * 0.5);

      // get cutoff
      double cutoff = 0.5;
      if (args.count("cutoff") > 0) {
        if(!(istringstream(args["cutoff"]) >> cutoff))
          throw gromos::Exception(argv[0], "Cutoff is not numeric!");
      }
      const double cutoff2 = cutoff * cutoff * 100.0; // to angstrom

      clipper::Grid_range gd(outmap.cell(), outmap.grid_sampling(), cutoff * 10.0);

      // loop over filter atoms
      for(int i = 0; i < atoms.size(); ++i) {
        // put atom into cell
        Vec pos = atoms.pos(i) - pbc.nearestImage(atoms.pos(i), centre, sys.box()) + centre;
        clipper::Coord_orth atom_pos(pos[0] * 10.0, pos[1] * 10.0, pos[2] * 10.0);
        clipper::Coord_frac uvw = atom_pos.coord_frac(outmap.cell());

        clipper::Coord_grid g0 = uvw.coord_grid(outmap.grid_sampling()) + gd.min();
        clipper::Coord_grid g1 = uvw.coord_grid(outmap.grid_sampling()) + gd.max();

        clipper::Xmap<clipper::ftype64>::Map_reference_coord i0, iu, iv, iw;
        clipper::Xmap<clipper::ftype64>::Map_reference_coord i0_res, iu_res, iv_res, iw_res;

        i0 = clipper::Xmap<clipper::ftype64>::Map_reference_coord(outmap, g0);
        i0_res = clipper::Xmap<clipper::ftype64>::Map_reference_coord(resmap, g0);
        // loop over grid and convolve with the atomic density gradient
        for (iu = i0, iu_res = i0_res; iu.coord().u() <= g1.u(); iu.next_u(), iu_res.next_u()) {
          for (iv = iu, iv_res = iu_res; iv.coord().v() <= g1.v(); iv.next_v(), iv_res.next_v()) {
            for (iw = iv, iw_res = iv_res; iw.coord().w() <= g1.w(); iw.next_w(), iw_res.next_w()) {
              const double d2 = (iw.coord_orth() - atom_pos).lengthsq();
              if (d2 <= cutoff2) {
                outmap[iw] = resmap[iw_res];
              } // in cutoff
            }
          }
        } // loop over grid
      } // for filter atoms
      
    }

    try {
      clipper::CCP4MAPfile file;
      file.open_write(args["out"]);
      file.export_xmap(outmap);
      file.close_write();
    } catch (const clipper::Message_fatal & msg) {
      throw gromos::Exception(argv[0], msg.text());
    }


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

