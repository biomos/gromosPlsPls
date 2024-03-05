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
 * @file pdb2g96.cc
 * Converts coordinate files from pdb to GROMOS format
 */

/**
 * @page programs Program Documentation
 *
 * @anchor pdb2g96
 * @section pdb2g96 Converts coordinate files from pdb to GROMOS format
 * @author @ref vk @ref mp
 * @date 7-6-07, 28-02-2017
 *
 * Converts a pdb-file (Protein Data Bank) into GROMOS coordinates. The unit of
 * the coordinates is converted from Angstrom to nm. The order of the atoms in
 * the pdbfile does not necessarily correspond to the order of the atoms in the
 * molecular topology file, but the residues should come in the proper order.
 * The program identifies atoms and residues based on their names, alternatives
 * to the atom and residue names in the topology can be specified in a library
 * file (see Volume IV). The only requirement on residue numbers in the
 * pdb-file is that the residue number should change when going from one
 * residue to the next. Mismatches between the topology and the pdb-file are
 * treated as follows:
 * <ol>
 * <li> If the expected residue according to the topology is not found, a
 *      warning is written out and the next residue in the pdb-file is read in
 *      until a match with the topology is found.
 * <li> Atoms that are expected according to the topology, but that are not
 *      found in the pdb-file are either generated (if they are hydrogens and
 \@gch is given)
 *       or written out in the coordinate file with
 *      coordinates (0.0, 0.0, 0.0). In that case a warning is written to cerr.
 * <li> Atoms that are present in the pdb-file, but not expected according to
 *      the topology are ignored, a warning is written to cerr.
 * </ol>
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pdb</td><td>&lt;pdb coordinates&gt; </td></tr>
 * <tr><td> \@out</td><td>&lt;resulting GROMOS coordinates&gt; (optional,
 defaults to stdout) </td></tr>
 * <tr><td> \@lib</td><td>&lt;library for atom and residue names&gt; </td></tr>
 * <tr><td>[\@gch</td><td>&lt;(re)generate hydrogen coordinates&gt;] </td></tr>
 * <tr><td>[\@tol</td><td>&lt;tolerance (default 0.1 %)&gt;] </td></tr>
 * <tr><td>[\@pbc</td><td>&lt;boundary type&gt; &lt;gather method&gt;]
 </td></tr>
 * <tr><td>[\@outbf</td><td>&lt;write B factors and occupancies to an additional
 file&gt;]</td></tr>
 * <tr><td>[\@factor</td><td>&lt;factor to convert lentgh unit to
 Angstrom&gt;]</td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  pdb2g96
    @topo  ex.top
    @pdb   exref.pdb
    @pbc   v
    @tol   0.1
    @gch
    @out   ex.g96
    @lib   pdb2g96.lib
 @endverbatim
 *
 * <hr>
 */

/* pdb2g96.cc  This program reads in a topology and a pdb file.
 *             it will then try to generate a gromos-coordinate file
 *             with the atoms in the correct order
 *
 *
 * Only parses the ATOM and HETATM records from the coordinate section.
 * See
 * http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
 * for a specification of the PDB file format.
 * For ATOMs and HETATMs we find:
 */

/* #define ATOM substr(0, 4) */
/* #define HETATM substr(0, 6) */
// the following two are not strictly standard
/* #define ATOMNAME substr(12, 5) */
/* #define RESNAME substr(17, 4) */
//
/* #define RESNUM substr(22, 4) */
/* #define COORDX substr(30, 8) */
/* #define COORDY substr(38, 8) */
/* #define COORDZ substr(46, 8) */
/* #define OCCUPANCY substr(54, 6) */
/* #define BFACTOR substr(60, 6) */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/System.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/OutG96S.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/Gch.h"

using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace args;
using namespace std;
using namespace utils;
using namespace bound;

auto ATOMNAME = [](const string &str) { return str.substr(12, 5); };
auto RESNAME = [](const string &str) { return str.substr(17, 4); };
auto RESNUM = [](const string &str) { return str.substr(22, 4); };

struct InPDBLine {
  string atom() { return line.substr(0, 4); };
  string hetatm() { return line.substr(0, 6); };
  double coordx() { return stod(line.substr(30, 8)); };
  double coordy() { return stod(line.substr(38, 8)); };
  double coordz() { return stod(line.substr(46, 8)); };
  double occupancy() { return stod(line.substr(54, 6)); };
  double bfactor() { return stod(line.substr(60, 6)); };
  string line;
};

/* ugly but functional c++ hack to strip whitespace
 * from a string */
string stripWhite(string s) {

  istringstream bla(s.c_str());
  string fasel;
  bla >> fasel;
  return fasel;
}
/*
 * checks if two names are the same or should be
 * considered the same
 */
bool checkName(multimap<string, string> lib, string nameA, string nameB) {
  nameA = stripWhite(nameA);
  if (nameA == nameB)
    return true;
  for (multimap<string, string>::const_iterator iter = lib.lower_bound(nameA),
                                                to = lib.upper_bound(nameA);
       iter != to; ++iter)
    if (iter->second == nameB)
      return true;
  return false;
}

bool isnan(Vec v) {
  return std::isnan(v[0]) || std::isnan(v[1]) || std::isnan(v[2]);
}

vector<string> nextPdbResidue(list<vector<string>> &pdbResidues) {

  vector<string> pdbResidue;

  if (pdbResidues.begin() != pdbResidues.end())
    pdbResidue = *(pdbResidues.begin());

  return pdbResidue;
}

void checkResidueName(vector<string> pdbResidue, string resName,
                      const multimap<string, string> &libRes) {

  if (!pdbResidue.size()) {
    ostringstream os;
    os << "Error: Empty Residue.\n"
       << "No coordinates in pdb file.";

    throw gromos::Exception("pdb2g96", os.str());
  }

  if (!checkName(libRes, RESNAME(pdbResidue[0]), resName)) {
    ostringstream os;
    os << "Error: Residue names do not match.\n"
       << "\tIn topology: " << resName
       << ", in pdb file: " << RESNAME(pdbResidue[0]);

    throw gromos::Exception("pdb2g96", os.str());
  }
}

void warnNotFoundAtom(int atomNum, string atomName, int resNum,
                      string resName) {

  cerr << "Warning: Could not find atom " << atomNum + 1 << " (" << atomName
       << ")"
       << ","

       << " in residue " << resNum + 1 << " (" << resName << ")"
       << ".\n"

       << "\tSet coordinates to (0.0 0.0 0.0)\n";
}

void warnIgnoredAtoms(vector<string> pdbAtoms) {

  for (unsigned int lineNum = 0; lineNum < pdbAtoms.size(); lineNum++)

    cerr << "Warning: Ignored atom " << ATOMNAME(stripWhite(pdbAtoms[lineNum]))

         << " in residue " << stripWhite(RESNUM(pdbAtoms[lineNum])) << " ("
         << stripWhite(RESNAME(pdbAtoms[lineNum])) << ").\n";
}

void readLibrary(Ginstream &lib, multimap<string, string> &libRes,
                 map<string, multimap<string, string>> &libAtom) {

  typedef multimap<string, string>::value_type mapType;

  std::vector<std::string> buffer;
  std::vector<std::vector<std::string>> content;
  while (!lib.stream().eof()) {
    lib.getblock(buffer);
    if (!lib.stream().eof()) {
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("pdb2g96", "Library file " + lib.name() +
                                               " is corrupted. No END in " +
                                               buffer[0] + " block. Got\n" +
                                               buffer[buffer.size() - 1]);

      content.push_back(buffer);
    }
  }
  // now loop over the content
  std::vector<std::vector<std::string>>::const_iterator iter = content.begin();
  string resNameA, resNameB, atomNameA, atomNameB;

  for (; iter != content.end(); ++iter) {
    if ((*iter)[0] == "RESIDUES" || (*iter)[0] == "RESIDUENAMELIB") {
      for (unsigned int i = 1; i < (*iter).size() - 1; i++) {
        std::istringstream linestream((*iter)[i]);
        linestream >> resNameA >> resNameB;
        libRes.insert(mapType(resNameA, resNameB));
      }

    } else if ((*iter)[0] == "ATOMS" || (*iter)[0] == "ATOMNAMELIB") {
      for (unsigned int i = 1; i < (*iter).size() - 1; i++) {
        std::istringstream linestream((*iter)[i]);
        linestream >> resNameA >> atomNameA >> atomNameB;
        libAtom[resNameA].insert(mapType(atomNameA, atomNameB));
      }
    } else
      throw gromos::Exception("pdb2g96", "Don't know how to handle " +
                                             (*iter)[0] + "-block in library");
  }
}

vector<string> split(const string &s, char delim) {
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

/* loop over all molecules */
void molecules(System &sys, list<vector<string>> &pdbResidues,
               const multimap<string, string> &libRes,
               map<string, multimap<string, string>> &libAtom,
               const double &fromang, bool do_bfactors, ofstream &bf_file,
               const Arguments &args) {
  for (int molNum = 0; molNum < sys.numMolecules(); molNum++) {

    // reserve memory for the coordinates
    sys.mol(molNum).initPos();
    // determine which are hydrogens based on the mass
    sys.mol(molNum).topology().setHmass(1.008);

    // loop over all residues
    int firstAtomNum = 0, lastAtomNum = 0;
    for (int resNum = 0; resNum != sys.mol(molNum).topology().numRes();
         ++resNum) {

      vector<string> pdbResidue = nextPdbResidue(pdbResidues);
      // if the residues in the pdb and the topology are
      // not identical, skip this loop.
      try {
        checkResidueName(pdbResidue, sys.mol(molNum).topology().resName(resNum),
                         libRes);
        pdbResidues.pop_front();
      } catch (gromos::Exception &e) {
        cerr << e.what() << endl;
        cerr << " Could not read residue number " << resNum + 1;
        cerr << " of molecule " << molNum + 1;
        cerr << " from pdb file." << endl;
        cerr << "Skipped" << endl;
      }

      /*
       * determine the first and the last atom number of
       * this residue in the topology
       */
      for (firstAtomNum = lastAtomNum;
           sys.mol(molNum).topology().resNum(firstAtomNum) != resNum;
           firstAtomNum++)
        ;

      for (lastAtomNum = firstAtomNum;
           lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
           sys.mol(molNum).topology().resNum(lastAtomNum) == resNum;
           lastAtomNum++)
        ;

      /*
       * for every atom in the topology residue,
       * look for an atom in the pdb residue with the same name,
       * and import its coordinates. If we can't find one,
       * set the coords to 0,0,0 and issue a warning.
       */
      for (int atomNum = firstAtomNum; atomNum < lastAtomNum; atomNum++) {

        bool foundAtom = false;

        for (unsigned int pdbAtomNum = 0; pdbAtomNum < pdbResidue.size();
             pdbAtomNum++) {

          InPDBLine inPdbLine{pdbResidue[pdbAtomNum]};

          if (checkName(libAtom[sys.mol(molNum).topology().resName(resNum)],
                        ATOMNAME(inPdbLine.line),
                        sys.mol(molNum).topology().atom(atomNum).name()) &&
              !foundAtom) {

            foundAtom = true;
            sys.mol(molNum).pos(atomNum) =
                Vec(fromang * inPdbLine.coordx(), fromang * inPdbLine.coordy(),
                    fromang * inPdbLine.coordz());

            if (do_bfactors) {
              int res = sys.mol(molNum).topology().resNum(atomNum);
              bf_file << "# " << setw(5) << molNum + 1 << setw(5) << res + 1
                      << setw(5) << sys.mol(molNum).topology().resName(res)
                      << setw(5) << atomNum + 1 << setw(5)
                      << sys.mol(molNum).topology().atom(atomNum).name()
                      << endl;
              bf_file << setw(15) << fromang * fromang * inPdbLine.bfactor()
                      << setw(15) << inPdbLine.occupancy() << endl;
            }

            pdbResidue.erase(pdbResidue.begin() + pdbAtomNum);
          }
        }
        if (!foundAtom) {
          sys.mol(molNum).pos(atomNum) = Vec(0.0, 0.0, 0.0);
          if (do_bfactors) {
            int res = sys.mol(molNum).topology().resNum(atomNum);
            bf_file << "# " << setw(5) << molNum + 1 << setw(5) << res + 1
                    << setw(5) << sys.mol(molNum).topology().resName(res)
                    << setw(5) << atomNum + 1 << setw(5)
                    << sys.mol(molNum).topology().atom(atomNum).name()
                    << ": not found!" << endl
                    << setw(15) << 0.01 << setw(15) << 0.0 << endl;
          }
          // if we are adding hydrogen positions later, warn only if it is not
          // a hydrogen
          if (args.count("gch") < 0 ||
              !sys.mol(molNum).topology().atom(atomNum).isH()) {
            warnNotFoundAtom(
                atomNum, sys.mol(molNum).topology().atom(atomNum).name(),
                resNum, sys.mol(molNum).topology().resName(resNum));
          }
        }
      }
      // print a warning for the pdb atoms that were ignored
      warnIgnoredAtoms(pdbResidue);
    }
  }
}

/* This should be it */
/* everything that is left over, should be solvent */
void solvent(System &sys, list<vector<string>> &pdbResidues,
             const multimap<string, string> &libRes,
             map<string, multimap<string, string>> &libAtom,
             const double &fromang, bool do_bfactors, ofstream &bf_file,
             const Arguments &args) {
  int resNum = 0;
  sys.sol(0).topology().setHmass(1.008);

  for (int resNum = 0; resNum != pdbResidues.size(); ++resNum) {
    vector<string> pdbResidue = nextPdbResidue(pdbResidues);

    try {
      pdbResidues.pop_front();
      checkResidueName(pdbResidue, "SOLV", libRes);
    } catch (gromos::Exception &e) {
      cerr << e.what() << endl;
      cerr << " Could not read residue number " << resNum + 1;
      cerr << " of the solvent from pdb file." << endl;
      cerr << "Skipped" << endl;
      continue;
    }
    /*
     * for every atom in the topology residue,
     * look for an atom in the pdb residue with the same name,
     * and import its coordinates. If we can't find one,
     * set the coords to 0,0,0 and issue a warning.
     */
    for (int atomNum = 0; atomNum < sys.sol(0).topology().numAtoms();
         atomNum++) {

      bool foundAtom = false;

      for (unsigned int pdbAtomNum = 0; pdbAtomNum < pdbResidue.size();
           pdbAtomNum++) {

        InPDBLine inPdbLine{pdbResidue[pdbAtomNum]};

        if (checkName(libAtom["SOLV"], ATOMNAME(inPdbLine.line),
                      sys.sol(0).topology().atom(atomNum).name()) &&
            !foundAtom) {

          foundAtom = true;
          sys.sol(0).addPos(Vec(fromang * inPdbLine.coordx(),
                                fromang * inPdbLine.coordy(),
                                fromang * inPdbLine.coordz()));
          pdbResidue.erase(pdbResidue.begin() + pdbAtomNum);
          if (do_bfactors) {
            bf_file << "# " << setw(5) << "SOLV" << atomNum + 1 << endl;
            bf_file << setw(15) << fromang * fromang * inPdbLine.bfactor()
                    << setw(15) << inPdbLine.occupancy() << endl;
          }
        }
      }
      if (!foundAtom) {
        sys.sol(0).addPos(Vec(0.0, 0.0, 0.0));
        if (do_bfactors) {
          bf_file << "# " << setw(5) << "SOLV" << atomNum + 1 << ": not found!"
                  << endl;
          bf_file << setw(15) << 0.01 << setw(15) << 0.0 << endl;
        }
        if (args.count("gch") < 0 ||
            !sys.sol(0).topology().atom(atomNum).isH()) {
          warnNotFoundAtom(atomNum, sys.sol(0).topology().atom(atomNum).name(),
                           resNum, "SOLV");
        }
      }
    }

    // print a warning for the pdb atoms that were ignored
    warnIgnoredAtoms(pdbResidue);
  }
}

/* Open and read pdb file
 * Reads the ATOM and HETATM lines from a pdb file into
 * a list<vector<string>>, one vector<string> per residue.
 */

list<vector<string>> readPdbAtoms(const Arguments &args) {
  ifstream pdbFile(args["pdb"]);
  if (!pdbFile.good()) {
    throw gromos::Exception("Ginstream",
                            "Could not open file '" + args["pdb"] + "'");
  }
  if (!pdbFile.is_open()) {
    throw gromos::Exception("Ginstream",
                            "could not open file '" + args["pdb"] + "'");
  }

  string resNum = "    ";
  InPDBLine inPdbLine;
  vector<string> pdbResidue;
  list<vector<string>> pdbResidues;

  while (!pdbFile.eof()) {
    getline(pdbFile, inPdbLine.line);
    if (inPdbLine.atom() == "ATOM" || inPdbLine.hetatm() == "HETATM") {

      // check if we're in a new residue
      if (RESNUM(inPdbLine.line) != resNum) {

        resNum = RESNUM(inPdbLine.line);

        // if we're not in the first residue
        if (pdbResidue.size())
          pdbResidues.push_back(pdbResidue);

        pdbResidue.clear();
      }
      pdbResidue.push_back(inPdbLine.line);
    }
  }
  pdbFile.close();

  // push the last residue
  pdbResidues.push_back(pdbResidue);

  return pdbResidues;
}

void wrap(Arguments &args, System &sys) {

  list<vector<string>> pdbResidues{readPdbAtoms(args)};

  // read the library file
  std::multimap<std::string, std::string> libRes;
  std::map<std::string, std::multimap<std::string, std::string>> libAtom;
  if (args.count("lib") > 0) {
    Ginstream lib(args["lib"]);
    cerr << "# using library file " << args["lib"] << endl;
    readLibrary(lib, libRes, libAtom);
  }

  bool do_bfactors = false;
  ofstream bf_file;
  if (args.count("outbf") > 0) {
    bf_file.open(args["outbf"].c_str());
    if (!bf_file.is_open()) {
      throw gromos::Exception("pdb2g96",
                              "Cannot open @outbf file for writing.");
    }
    do_bfactors = true;
    bf_file << "TITLE\nB-Factors and occupancies\n\nEND\nBFACTOROCCUPANCY\n";
  }
  std::map<std::pair<int, int>, BFactorOccupancyData> bfactors;

  // get the factor
  double factor = args.getValue<double>("factor", false, 10.0);
  double fromang = 1.0 / factor;

  molecules(sys, pdbResidues, libRes, libAtom, fromang, do_bfactors, bf_file,
            args);
  solvent(sys, pdbResidues, libRes, libAtom, fromang, do_bfactors, bf_file,
          args);

  if (do_bfactors) {
    bf_file << "END\n";
  }
}

void out_pdb2g96_file(Arguments &args, System &sys) {
  InTopology it(args["topo"]);
  // that's really it for pdb2g96
  System outSys(sys);
  ostringstream os;
  os << "pdb2g96: Reordered atoms from " << args["pdb"];

  if (args.count("gch") >= 0) {
    //***************
    // gch: generate H coordinates
    //***************
    GromosForceField gff(it.forceField());

    // read in the accuracy
    double eps = args.getValue<double>("tol", false, 0.1) / 100.0;

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys, sys, args);

    // gather the system!
    (*pbc.*gathmethod)();

    // a bit of an ugly hack for solvent molecules. The whole program has
    // been written for solutes, but sometimes we would have crystallographic
    // waters for which we want to do the same thing.
    if (sys.sol(0).numPos()) {
      MoleculeTopology mt;
      for (int a = 0; a < sys.sol(0).topology().numAtoms(); a++) {
        mt.addAtom(sys.sol(0).topology().atom(a));
        mt.setResNum(a, 0);
      }
      mt.setResName(0, "SOLV");

      ConstraintIterator ci(sys.sol(0).topology());
      for (; ci; ++ci) {
        gff.addBondType(BondType(gff.numBondTypes(), 1, ci().dist()));
        Bond b(ci()[0], ci()[1], false);
        b.setType(gff.numBondTypes() - 1);
        mt.addBond(b);
      }

      // add every solvent molecule as a solute
      int numSolvent = sys.sol(0).numPos() / sys.sol(0).topology().numAtoms();
      for (int i = 0; i < numSolvent; i++) {
        Molecule m(mt);
        m.initPos();
        for (int a = 0; a < mt.numAtoms(); a++) {
          m.pos(a) = sys.sol(0).pos(i * sys.sol(0).topology().numAtoms() + a);
        }
        sys.addMolecule(m);
      }
      // and remove the original solvent
      sys.sol(0).setNumPos(0);
    }

    // initialize two counters
    int replaced = 0, kept = 0;

    // loop over all atoms
    for (int m = 0; m < sys.numMolecules(); m++) {

      // flag the atoms with mass 1.008 as hydrogens
      sys.mol(m).topology().setHmass(1.008);
      for (int a = 0; a < sys.mol(m).numAtoms(); a++) {

        if (!sys.mol(m).topology().atom(a).isH()) {

          // divide into hydrogens and non-hydrogens
          vector<int> h;
          vector<int> nh;
          get_h_nh_neighbours(sys, gff, m, a, h, nh);

          // only continue if we have hydrogens
          int numH = h.size();
          int numNH = nh.size();
          int geom = get_geometry(numH, numNH);
          if (numH && !geom) {
            ostringstream os;
            os << "Unexpected geometry for hydrogen bound to atom: " << m + 1
               << ":" << a + 1 << endl;
            throw(gromos::Exception("pdb2g96", os.str()));
          }
          // we have to have a geometry (this means that there are hydrogens)
          // and a should not be a hydrogen itself. (in the case of H2O we
          // have e.g. H-H bonds. These are only treated via the oxygen.
          if (geom) {
            int r = generate_hcoordinates(sys, gff, m, a, h, nh, geom, eps);
            replaced += r;
            kept += (numH - r);
          }
        }
      }
    }

    bool found_nan = false;
    // fix the solvent hack
    int solventIndex = 0;
    for (int m = 0; m < sys.numMolecules(); ++m) {
      if (m < outSys.numMolecules()) {
        // solute
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a) {
          if (isnan(sys.mol(m).pos(a))) {
            std::cerr << "Warning: " << sys.mol(m).topology().resNum(a) + 1
                      << " "
                      << sys.mol(m).topology().resName(
                             sys.mol(m).topology().resNum(a))
                      << " " << sys.mol(m).topology().atom(a).name() << " "
                      << v2s(sys.mol(m).pos(a)) << std::endl;
            found_nan = true;
          }
          outSys.mol(m).pos(a) = sys.mol(m).pos(a);
        }
      } else {
        // solvent
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++solventIndex) {
          if (isnan(sys.mol(m).pos(a))) {
            std::cerr << "Warning: " << sys.mol(m).topology().resNum(a) + 1
                      << " "
                      << sys.mol(m).topology().resName(
                             sys.mol(m).topology().resNum(a))
                      << " " << sys.mol(m).topology().atom(a).name() << " "
                      << v2s(sys.mol(m).pos(a)) << std::endl;
            found_nan = true;
          }
          outSys.sol(0).pos(solventIndex) = sys.mol(m).pos(a);
        }
      }
    }

    if (found_nan) {
      std::cerr << "WARNING: Some positions could not be generated (nan).\n  "
                   "This can e.g. be due to missing or overlapping "
                   "heteroatom positions.\n";
    }

    os << "\nFound " << replaced + kept << " hydrogen atoms " << endl;
    os << kept << " were within " << eps * 100 << "% of minimum energy bond "
       << "length" << endl;
    os << replaced << " were assigned new coordinates based on geometry";
  }

  // now define an output stream and write the coordinates
  OutG96S oc;
  ofstream fout;
  try {
    args.check("out", 1);
    fout.open(args["out"]);
    oc.open(fout);
  } catch (const gromos::Exception &e) {
    oc.open(cout);
  }
  oc.select("ALL");
  oc.writeTitle(os.str());
  oc << outSys;
  oc.close();
}

int main(int argc, char *argv[]) {
  Argument_List knowns;
  knowns << "topo"
         << "out"
         << "lib"
         << "outbf"
         << "factor"
         << "pdb"
         << "pbc"
         << "tol"
         << "gch";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <molecular topology file>\n";
  usage += "\t@pdb   <input coordinate file: pdb or gromos format, determined "
           "by extension>\n";
  usage += "\t[@out    <resulting GROMOS coordinates> (optional, defaults to "
           "stdout)]\n";
  usage += "\t@lib     <library for atom and residue names>\n";
  usage += "\t[@gch   <(re)generate hydrogen coordinates>]\n";
  usage += "\t[@pbc   <boundary type> <gather method>]\n";
  usage += "\t[@tol   <tolerance for gch (default 0.1 %)>]\n";
  usage +=
      "\t[@outbf  <write B factors and occupancies to an additional file>]\n";
  usage += "\t[@factor <factor to convert length unit to Angstrom, 10.0>]\n";

  Arguments args(argc, argv, knowns, usage);
  try {

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    wrap(args, sys);
    out_pdb2g96_file(args, sys);
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
