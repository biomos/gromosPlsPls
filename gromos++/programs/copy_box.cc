/*
 * This file is part of GROMOS.
 * 
 * Copyright (c) 2011, 2012, 2016, 2018, 2021, 2023 Biomos b.v.
 * See <https://www.gromos.net> for details.
 * 
 * GROMOS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/**
 * @file copy_box.cc
 * Repeats a simulation box in a given direction
 */

/**
 * @page programs Program Documentation
 *
 * @anchor copy_box
 * @section copy_box Repeats a simulation box in a given direction
 * @author @ref co
 * @date 8-6-07
 *
 * Program copy box can be used to duplicate the size of a system in the x, y 
 * or z direction (or k,l,m for triclinic boxes). This is especially convenient
 * if one wants to multiply the size of a system under periodic boundary
 * conditions in which the central box has a rectangular or triclinic shape.
 * The number of new copies can be set with extra :N option, e.g. x:2 will create
 * two more copies in the x direction. x:2,y:2 will create all combinations
 * totalling in a grid of 3x3 copies.
 * If one wants perform more elaborate transformations, the program @ref cry might
 * be of use (see section V-2.17). Note that program @ref com_top (see section
 * V-2.2) can be useful to additionally duplicate the solute block in the 
 * topology.
 * Note that the \@pbc flag is optional. Only if this flag is given, gathering
 * of the molecules will be performed before copying the box!
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pos</td><td>&lt;input coordinate file&gt; </td></tr>
 * <tr><td> \@dir</td><td>&lt;coordinate to duplicate: x/y/z/k/l/m  with optional :N&gt; </td></tr>
 * <tr><td> [\@pbc</td><td>&lt;boundary type&gt; </td></tr>;] 
 * </table>
 *
 *
 * Example:
 * @verbatim
  copy_box
    @topo ex.top
    @pos  exref.coo
    @dir  x:2,y:2
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cstdlib>
#include <string>
#include <sstream>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "../src/args/Arguments.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Vec.h"
#include "../src/bound/Boundary.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gromos/Exception.h"

using namespace std;
using namespace gcore;
using namespace gio;
using namespace fit;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace utils;

typedef map<char,int> DirMap;

DirMap parseDir(const string &dirstr)
{
  DirMap dirs;
  bool has_xyz = false;
  bool has_klm = false;

  string s = dirstr;
  replace(s.begin(), s.end(), ',', ' ');
  istringstream iss(s);

  string token;
  while (iss >> token) {
    int mult = 1;
    size_t colon = token.find(':');
    if (colon != string::npos) {
      mult = atoi(token.substr(colon + 1).c_str());
      if (mult < 1)
        throw gromos::Exception("copy_box",
          "Multiplicity in @dir must be >= 1");
      token = token.substr(0, colon);
    }

    for (char c : token) {
      if (string("xyz").find(c) != string::npos) has_xyz = true;
      else if (string("klm").find(c) != string::npos) has_klm = true;
      else
        throw gromos::Exception("copy_box",
          "Invalid direction specifier in @dir");

      dirs[c] += mult;
    }
  }

  if (has_xyz && has_klm)
    throw gromos::Exception("copy_box",
      "Cannot mix x/y/z with k/l/m in @dir");

  return dirs;
}

Vec shiftVector(const System &sys, char d)
{
  if (d == 'x') return Vec(sys.box().K()[0], 0, 0);
  if (d == 'y') return Vec(0, sys.box().L()[1], 0);
  if (d == 'z') return Vec(0, 0, sys.box().M()[2]);
  if (d == 'k') return sys.box().K();
  if (d == 'l') return sys.box().L();
  if (d == 'm') return sys.box().M();

  throw gromos::Exception("copy_box", "Invalid direction");
}

int main(int argc, char **argv)
{
  Argument_List knowns;
  knowns << "topo" << "pos" << "dir" << "pbc";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <molecular topology file>\n";
  usage += "\t@pos  <input coordinate file>\n";
  usage += "\t@dir  <coordinates to duplicate: x,y,z,k,l,m with optional :N>\n";
  usage += "\t[@pbc <boundary type> [<gathermethod>]]\n";

  try {
    Arguments args(argc, argv, knowns, usage);

    args.check("dir",1);
    string dirstr = args["dir"];
    DirMap dirs = parseDir(dirstr);

    
    // read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    System sys(it.system());
    System refSys(it.system());
    
    // read single atom coordinates...
    InG96 ic;
    ic.open(args["pos"]);
    ic.select("ALL");
    ic >> sys;
    ostringstream title;
    title << ic.title();
    ic.close();
    // Check if @pbc is given
    if (args.count("pbc") > 0) {
      Boundary *pbc = BoundaryParser::boundary(sys, args);
      Boundary::MemPtr gathmethod =
        args::GatherParser::parse(sys, refSys, args);
      (*pbc.*gathmethod)();

      for (auto &d : dirs) {
        if ((d.first=='x'||d.first=='y'||d.first=='z') &&
            pbc->type()!='r') {
          cerr << "WARNING: copying non-rectangular box along x/y/z\n";
        }
        if ((d.first=='k'||d.first=='l'||d.first=='m') &&
            pbc->type()!='c') {
          cerr << "WARNING: copying non-triclinic box along k/l/m\n";
        }
      }
    }

    vector<System> images;
    images.push_back(sys);

    for (auto &d : dirs) {
      char axis = d.first;
      int n = d.second;
      Vec baseShift = shiftVector(sys, axis);

      vector<System> newImages;
      for (int i = 1; i <= n; i++) {
        Vec shift = baseShift * i;
        for (const auto &ref : images) {
          System tmp(ref);
          PositionUtils::translate(&tmp, shift);
          newImages.push_back(tmp);
        }
      }
      images.insert(images.end(), newImages.begin(), newImages.end());
    }

    System finalSys(images[0]);

    for (size_t im = 1; im < images.size(); im++) {
      for (int i = 0; i < images[im].numMolecules(); i++)
        finalSys.addMolecule(images[im].mol(i));

      for (int i = 0; i < images[im].sol(0).numPos(); i++)
        finalSys.sol(0).addPos(images[im].sol(0).pos(i));
    }

    for (auto &d : dirs) {
      Vec shift = shiftVector(sys, d.first) * d.second;

      if (d.first=='x') finalSys.box().K()[0] += shift[0];
      else if (d.first=='y') finalSys.box().L()[1] += shift[1];
      else if (d.first=='z') finalSys.box().M()[2] += shift[2];
      else if (d.first=='k') finalSys.box().K() += shift;
      else if (d.first=='l') finalSys.box().L() += shift;
      else if (d.first=='m') finalSys.box().M() += shift;
    }

    OutG96S oc;
    title << "\nCopy_box: " << args["pos"]
          << " duplicated with dir=" << dirstr
          << " (" << images.size() << " copies)";

    oc.open(cout);
    oc.select("ALL");
    oc.writeTitle(title.str());
    oc << finalSys;
  }
  catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    return 1;
  }
  return 0;
}




