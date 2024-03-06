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
string stripWhite(string str) {

  istringstream bla(str);
  string fasel;
  bla >> fasel;
  return fasel;
}
/*
 * checks if two names are the same or should be
 * considered the same
 */
bool checkName(const multimap<string, string> &lib, string nameA,
               string nameB) {
  nameA = stripWhite(nameA);
  if (nameA == nameB)
    return true;
  for (auto iter = lib.lower_bound(nameA), to = lib.upper_bound(nameA);
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

using LibRes = multimap<string, string>;
using LibAtom = map<string, multimap<string, string>>;
class Library {
public:
  void read_library(Ginstream &lib) {

    using mapType = multimap<string, string>::value_type;

    vector<string> buffer;
    vector<vector<string>> content;
    while (!lib.stream().eof()) {
      lib.getblock(buffer);
      if (!lib.stream().eof()) {
        if (buffer[buffer.size() - 1].find("END") != 0) {
          throw gromos::Exception("pdb2g96", "Library file " + lib.name() +
                                                 " is corrupted. No END in " +
                                                 buffer[0] + " block. Got\n" +
                                                 buffer[buffer.size() - 1]);
        }
        content.push_back(buffer);
      }
    }
    // now loop over the content
    string resNameA, resNameB, atomNameA, atomNameB;
    for (const auto &line : content) {
      if (line[0] == "RESIDUES" || line[0] == "RESIDUENAMELIB") {
        for (auto iter = ++line.begin(); iter < line.end() - 1; ++iter) {
          istringstream linestream(*iter);
          linestream >> resNameA >> resNameB;
          libRes.insert(mapType(resNameA, resNameB));
        }

      } else if (line[0] == "ATOMS" || line[0] == "ATOMNAMELIB") {
        for (auto iter = ++line.begin(); iter < line.end() - 1; ++iter) {
          istringstream linestream(*iter);
          linestream >> resNameA >> atomNameA >> atomNameB;
          libAtom[resNameA].insert(mapType(atomNameA, atomNameB));
        }
      } else
        throw gromos::Exception("pdb2g96", "Don't know how to handle " +
                                               line[0] + "-block in library");
    }
  }
  const LibRes &get_residue() const { return libRes; }
  LibAtom get_atom() const { return libAtom; }

private:
  LibRes libRes;
  LibAtom libAtom;
};

vector<string> split(const string &s, char delim) {
  vector<string> elems;
  stringstream ss(s);
  string item;
  while (getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

struct Add_position {
  virtual ~Add_position() = default;
  virtual void add(System &sys, int molNum, int atomNum, const Vec &coord) = 0;
};

struct Add_solute_position : Add_position {
  void add(System &sys, int molNum, int atomNum, const Vec &coord) override {
    sys.mol(molNum).pos(atomNum) = coord;
  }
};

struct Add_solvent_position : Add_position {
  void add(System &sys, int molNum, int atomNum, const Vec &coord) override {
    sys.sol(0).addPos(coord);
  }
};

struct BFactor {
  virtual ~BFactor() = default;
  virtual void print(ofstream &bf_file, bool foundAtom, int molNum, int resNum,
                     const string &resName, int atomNum, const string &atomName,
                     const double &bfactor, const double &occupancy) = 0;
};

struct BFactor_solute : BFactor {
  void print(ofstream &bf_file, bool foundAtom, int molNum, int resNum,
             const string &resName, int atomNum, const string &atomName,
             const double &bfactor, const double &occupancy) override {
    bf_file << "# " << setw(5) << molNum + 1 << setw(5) << resNum + 1 << setw(5)
            << resName << setw(5) << atomNum + 1 << setw(5) << atomName;
    if (!foundAtom)
      bf_file << ": not found!";
    bf_file << '\n';
    bf_file << setw(15) << bfactor << setw(15) << occupancy << '\n';
  }
};

struct BFactor_solvent : BFactor {
  void print(ofstream &bf_file, bool foundAtom, int molNum, int resNum,
             const string &resName, int atomNum, const string &atomName,
             const double &bfactor, const double &occupancy) override {
    bf_file << "# " << setw(5) << resName << atomNum + 1;
    if (!foundAtom)
      bf_file << ": not found!";
    bf_file << '\n';
    bf_file << setw(15) << bfactor << setw(15) << occupancy << '\n';
  }
};

struct PdbResidue {
  void checkResidueName(const string &resName, const LibRes &libRes) const {

    if (!pdbAtom.size()) {
      ostringstream os;
      os << "Error: Empty Residue.\n"
         << "No coordinates in pdb file.";

      throw gromos::Exception("pdb2g96", os.str());
    }

    if (!checkName(libRes, resname_of_atom(0), resName)) {
      ostringstream os;
      os << "Error: Residue names do not match.\n"
         << "\tIn topology: " << resName
         << ", in pdb file: " << resname_of_atom(0);

      throw gromos::Exception("pdb2g96", os.str());
    }
  }
  void warnNotFoundAtom(int atomNum, const string &atomName, int resNum) const {

    cerr << "Warning: Could not find atom " << atomNum + 1 << " (" << atomName
         << ")"
         << ","

         << " in residue " << resNum + 1 << " (" << resName << ")"
         << ".\n"

         << "\tSet coordinates to (0.0 0.0 0.0)\n";
  }
  void Match_Atom(ofstream &bf_file, LibAtom libAtom, System &sys,
                  const Arguments &args, int molNum, int resNum, int atomNum,
                  double fromang, bool do_bfactors, bool b_solv) {
    string atomName;
    if (!b_solv) {
      resName = sys.mol(molNum).topology().resName(resNum);
      atomName = sys.mol(molNum).topology().atom(atomNum).name();
    } else {
      resName = "SOLV";
      atomName = sys.sol(0).topology().atom(atomNum).name();
    }
    bool foundAtom = false;

    for (auto it = pdbAtom.begin(); it < pdbAtom.end(); ++it) {

      InPDBLine inPdbLine{*it};

      if (checkName(libAtom[resName], ATOMNAME(inPdbLine.line), atomName) &&
          !foundAtom) {

        foundAtom = true;
        Add_position *ptr_position;
        if (!b_solv)
          ptr_position = new Add_solute_position;
        else
          ptr_position = new Add_solvent_position;

        ptr_position->add(sys, molNum, atomNum,
                          Vec(fromang * inPdbLine.coordx(),
                              fromang * inPdbLine.coordy(),
                              fromang * inPdbLine.coordz()));
        delete ptr_position;
        if (do_bfactors) {
          double bfactor = fromang * fromang * inPdbLine.bfactor();
          double occupancy = inPdbLine.occupancy();
          BFactor *ptr_bfactor;
          if (!b_solv)
            ptr_bfactor = new BFactor_solute;
          else
            ptr_bfactor = new BFactor_solvent;

          ptr_bfactor->print(bf_file, foundAtom, molNum, resNum, resName,
                             atomNum, atomName, bfactor, occupancy);
          delete ptr_bfactor;
        }
        pdbAtom.erase(it);
      }
    }
    if (!foundAtom) {
      Add_position *ptr_position;
      if (!b_solv)
        ptr_position = new Add_solute_position;
      else
        ptr_position = new Add_solvent_position;

      ptr_position->add(sys, molNum, atomNum, Vec(0.0, 0.0, 0.0));
      delete ptr_position;
      if (do_bfactors) {
        double bfactor = 0.01;
        double occupancy = 0.0;
        BFactor *ptr_bfactor;
        if (!b_solv)
          ptr_bfactor = new BFactor_solute;
        else
          ptr_bfactor = new BFactor_solvent;

        ptr_bfactor->print(bf_file, foundAtom, molNum, resNum, resName, atomNum,
                           atomName, bfactor, occupancy);
        delete ptr_bfactor;
      }
      bool is_H;
      if (!b_solv)
        is_H = sys.mol(molNum).topology().atom(atomNum).isH();
      else
        is_H = sys.sol(0).topology().atom(atomNum).isH();

      if (args.count("gch") < 0 || !is_H)
        warnNotFoundAtom(atomNum, atomName, resNum);
    }
  }
  void warnIgnoredAtoms() const {

    for (unsigned int lineNum = 0; lineNum < pdbAtom.size(); lineNum++)

      cerr << "Warning: Ignored atom " << ATOMNAME(stripWhite(pdbAtom[lineNum]))

           << " in residue " << stripWhite(RESNUM(pdbAtom[lineNum])) << " ("
           << stripWhite(resname_of_atom(lineNum)) << ").\n";
  }
  string resname_of_atom(int index) const {
    return pdbAtom[index].substr(17, 4);
  }

  vector<string> pdbAtom;
  string resName;
};

class PdbResidues {
public:
  /* Read atoms from pdb file into lines */
  void readPdbAtoms(const Arguments &args) {
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
  }
  vector<string> nextPdbResidue() {

    vector<string> pdbResidue;

    if (pdbResidues.begin() != pdbResidues.end())
      pdbResidue = *pdbResidues.begin();

    return pdbResidue;
  }
  void parse_molecules(System &sys, const Library &lib, const double &fromang,
                       bool do_bfactors, ofstream &bf_file,
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

        PdbResidue pdbResidue{nextPdbResidue()};
        string resName = sys.mol(molNum).topology().resName(resNum);
        // if the residues in the pdb and the topology are
        // not identical, skip this loop.
        try {
          pdbResidue.checkResidueName(resName, lib.get_residue());
          pdbResidues.pop_front();
        } catch (gromos::Exception &e) {
          cerr << e.what() << endl;
          cerr << " Could not read residue number " << resNum + 1;
          cerr << " of molecule " << molNum + 1;
          cerr << " from pdb file." << endl;
          cerr << "Skipped" << endl;
          continue; /* Missing continue */
        }

        /*
         * determine the first and the last atom number of
         * this residue in the topology
         */
        firstAtomNum = lastAtomNum;
        while (sys.mol(molNum).topology().resNum(firstAtomNum) != resNum)
          firstAtomNum++;

        lastAtomNum = firstAtomNum;
        while (lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
               sys.mol(molNum).topology().resNum(lastAtomNum) == resNum)
          lastAtomNum++;

        /*
         * for every atom in the topology residue,
         * look for an atom in the pdb residue with the same name,
         * and import its coordinates. If we can't find one,
         * set the coords to 0,0,0 and issue a warning.
         */
        for (int atomNum = firstAtomNum; atomNum < lastAtomNum; atomNum++) {
          pdbResidue.Match_Atom(bf_file, lib.get_atom(), sys, args, molNum,
                                resNum, atomNum, fromang, do_bfactors, false);
        }
        // print a warning for the pdb atoms that were ignored
        pdbResidue.warnIgnoredAtoms();
      }
    }
  }
  void parse_solvent(System &sys, const Library &lib, const double &fromang,
                     bool do_bfactors, ofstream &bf_file,
                     const Arguments &args) {

    sys.sol(0).topology().setHmass(1.008);

    int molNum = 0;
    string resName = "SOLV";
    for (int resNum = 0; resNum != pdbResidues.size(); ++resNum) {
      PdbResidue pdbResidue{nextPdbResidue()};

      try {
        pdbResidue.checkResidueName(resName, lib.get_residue());
        pdbResidues.pop_front();
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
        pdbResidue.Match_Atom(bf_file, lib.get_atom(), sys, args, molNum,
                              resNum, atomNum, fromang, do_bfactors, true);
      }
      // print a warning for the pdb atoms that were ignored
      pdbResidue.warnIgnoredAtoms();
    }
  }

private:
  list<vector<string>> pdbResidues;
  Library library;
};

/* loop over all molecules */
void molecules(System &sys, list<vector<string>> &pdbResidues,
               const Library &lib, const double &fromang, bool do_bfactors,
               ofstream &bf_file, const Arguments &args) {

  for (int molNum = 0; molNum < sys.numMolecules(); molNum++) {

    // reserve memory for the coordinates
    sys.mol(molNum).initPos();
    // determine which are hydrogens based on the mass
    sys.mol(molNum).topology().setHmass(1.008);

    // loop over all residues
    int firstAtomNum = 0, lastAtomNum = 0;
    for (int resNum = 0; resNum != sys.mol(molNum).topology().numRes();
         ++resNum) {

      PdbResidue pdbResidue{nextPdbResidue(pdbResidues)};
      string resName = sys.mol(molNum).topology().resName(resNum);
      // if the residues in the pdb and the topology are
      // not identical, skip this loop.
      try {
        pdbResidue.checkResidueName(resName, lib.get_residue());
        pdbResidues.pop_front();
      } catch (gromos::Exception &e) {
        cerr << e.what() << endl;
        cerr << " Could not read residue number " << resNum + 1;
        cerr << " of molecule " << molNum + 1;
        cerr << " from pdb file." << endl;
        cerr << "Skipped" << endl;
        continue; /* Missing continue */
      }

      /*
       * determine the first and the last atom number of
       * this residue in the topology
       */
      firstAtomNum = lastAtomNum;
      while (sys.mol(molNum).topology().resNum(firstAtomNum) != resNum)
        firstAtomNum++;

      lastAtomNum = firstAtomNum;
      while (lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
             sys.mol(molNum).topology().resNum(lastAtomNum) == resNum)
        lastAtomNum++;

      /*
       * for every atom in the topology residue,
       * look for an atom in the pdb residue with the same name,
       * and import its coordinates. If we can't find one,
       * set the coords to 0,0,0 and issue a warning.
       */
      for (int atomNum = firstAtomNum; atomNum < lastAtomNum; atomNum++) {
        pdbResidue.Match_Atom(bf_file, lib.get_atom(), sys, args, molNum,
                              resNum, atomNum, fromang, do_bfactors, false);
      }
      // print a warning for the pdb atoms that were ignored
      pdbResidue.warnIgnoredAtoms();
    }
  }
}

/* This should be it */
/* everything that is left over, should be solvent */
void solvent(System &sys, list<vector<string>> &pdbResidues, const Library &lib,
             const double &fromang, bool do_bfactors, ofstream &bf_file,
             const Arguments &args) {

  sys.sol(0).topology().setHmass(1.008);

  int molNum = 0;
  string resName = "SOLV";
  for (int resNum = 0; resNum != pdbResidues.size(); ++resNum) {
    PdbResidue pdbResidue{nextPdbResidue(pdbResidues)};

    try {
      pdbResidue.checkResidueName(resName, lib.get_residue());
      pdbResidues.pop_front();
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
      pdbResidue.Match_Atom(bf_file, lib.get_atom(), sys, args, molNum, resNum,
                            atomNum, fromang, do_bfactors, true);
    }
    // print a warning for the pdb atoms that were ignored
    pdbResidue.warnIgnoredAtoms();
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

  /* list<vector<string>> pdbResidues{readPdbAtoms(args)}; */
  PdbResidues pdbResidues;
  pdbResidues.readPdbAtoms(args);

  /* read the library file */
  Library library;
  if (args.count("lib") > 0) {
    Ginstream lib(args["lib"]);
    cerr << "# using library file " << args["lib"] << endl;
    library.read_library(lib);
  }

  /* Initiate writing of bfactors */
  bool do_bfactors = false;
  ofstream bf_file;
  if (args.count("outbf") > 0) {
    bf_file.open(args["outbf"]);
    if (!bf_file.is_open()) {
      throw gromos::Exception("pdb2g96",
                              "Cannot open @outbf file for writing.");
    }
    do_bfactors = true;
    bf_file << "TITLE\nB-Factors and occupancies\n\nEND\nBFACTOROCCUPANCY\n";
  }

  /* get the factor */
  double factor = args.getValue<double>("factor", false, 10.0);
  double fromang = 1.0 / factor;

  pdbResidues.parse_molecules(sys, library, fromang, do_bfactors, bf_file,
                              args);
  pdbResidues.parse_solvent(sys, library, fromang, do_bfactors, bf_file, args);
  /* molecules(sys, pdbResidues, library, fromang, do_bfactors, bf_file, args);
   */
  /* solvent(sys, pdbResidues, library, fromang, do_bfactors, bf_file, args); */

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
            cerr << "Warning: " << sys.mol(m).topology().resNum(a) + 1 << " "
                 << sys.mol(m).topology().resName(
                        sys.mol(m).topology().resNum(a))
                 << " " << sys.mol(m).topology().atom(a).name() << " "
                 << v2s(sys.mol(m).pos(a)) << endl;
            found_nan = true;
          }
          outSys.mol(m).pos(a) = sys.mol(m).pos(a);
        }
      } else {
        // solvent
        for (int a = 0; a < sys.mol(m).numAtoms(); ++a, ++solventIndex) {
          if (isnan(sys.mol(m).pos(a))) {
            cerr << "Warning: " << sys.mol(m).topology().resNum(a) + 1 << " "
                 << sys.mol(m).topology().resName(
                        sys.mol(m).topology().resNum(a))
                 << " " << sys.mol(m).topology().atom(a).name() << " "
                 << v2s(sys.mol(m).pos(a)) << endl;
            found_nan = true;
          }
          outSys.sol(0).pos(solventIndex) = sys.mol(m).pos(a);
        }
      }
    }

    if (found_nan) {
      cerr << "WARNING: Some positions could not be generated (nan).\n  "
              "This can e.g. be due to missing or overlapping "
              "heteroatom positions.\n";
    }

    os << "\nFound " << replaced + kept << " hydrogen atoms " << '\n';
    os << kept << " were within " << eps * 100 << "% of minimum energy bond "
       << "length" << '\n';
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
