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

#include <algorithm>
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

class InPDBLine {
public:
  void read_line(ifstream &pdbFile) { getline(pdbFile, line); }
  string atom() const { return line.substr(0, 4); }
  string hetatm() const { return line.substr(0, 6); }
  string atomname() const { return line.substr(12, 5); }
  string resname() const { return line.substr(17, 4); }
  int resnum() const { return stoi(line.substr(22, 4)); }
  double coordx() const { return stod(line.substr(30, 8)); }
  double coordy() const { return stod(line.substr(38, 8)); }
  double coordz() const { return stod(line.substr(46, 8)); }
  Vec coord() const { return {coordx(), coordy(), coordz()}; }
  double occupancy() const { return stod(line.substr(54, 6)); }
  double bfactor() const { return stod(line.substr(60, 6)); }
  const string &get_line() const { return line; }

private:
  string line;
};

struct InCIFLine {
  string atomname() const { return atomName; }
  string resname() const { return resName; }
  int resnum() const { return resNum; }
  Vec coord() const { return {coordx, coordy, coordz}; }
  double occupancy() const { return Occupancy; }
  double bfactor() const { return Bfactor; }

  string atomType;
  string atomName;
  string resName;
  double coordx;
  double coordy;
  double coordz;
  double Occupancy;
  double Bfactor;
  int resNum;
  /* string symbol; //TODO: Member data that could be found in mmCIF files but
   * are not used by gromos */
  /* string chain; */
  /* int atomNum; */
};

struct Atom {
  explicit Atom(const InPDBLine &inPDB)
      : atomName{inPDB.atomname()}, coord{inPDB.coord()},
        occupancy{inPDB.occupancy()}, bfactor{inPDB.bfactor()} {}
  explicit Atom(const InCIFLine &inCIF)
      : atomName{inCIF.atomname()}, coord{inCIF.coord()},
        occupancy{inCIF.occupancy()}, bfactor{inCIF.bfactor()} {}
  ~Atom() = default;
  void convert_units() {
    coord *= fromang;
    bfactor *= fromang * fromang;
  }

  string atomName;
  Vec coord;
  double occupancy;
  double bfactor;
  static double fromang;
};

double Atom::fromang = 1.0 / 10.0;

/* New stripWhite function in line with modern C++  */
string stripWhite(const string &str) {
  std::string result;
  result.reserve(str.size()); // Reserve space for efficiency

  const std::locale loc;
  // Copy only non-whitespace characters
  for (char c : str) {
    if (!std::isspace(c, loc))
      result.push_back(c);
  }
  return result;
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

struct Residue {
  void add_atom(Atom atom) {
    atom.convert_units();
    this->atom.push_back(atom);
  }
  void checkResidueName(const LibRes &libRes, const string &resName) const {

    if (!atom.size()) {
      ostringstream os;
      os << "Error: Empty Residue.\n"
         << "No coordinates in pdb file.";

      throw gromos::Exception("pdb2g96", os.str());
    }

    if (!checkName(libRes, this->resName, resName)) {
      ostringstream os;
      os << "Error: Residue names do not match.\n"
         << "\tIn topology: " << resName
         << ", in coordinate file: " << this->resName;
      throw gromos::Exception("pdb2g96", os.str());
    }
  }
  void Match_Atom(LibAtom libAtom, System &sys, const Arguments &args,
                  int molNum, int resNum, int atomNum, bool b_solv) {
    string atomName;
    if (!b_solv) {
      resName = sys.mol(molNum).topology().resName(resNum);
      atomName = sys.mol(molNum).topology().atom(atomNum).name();
    } else {
      resName = "SOLV";
      atomName = sys.sol(0).topology().atom(atomNum).name();
    }

    bool foundAtom = false;
    for (auto it = atom.begin(); it < atom.end(); ++it) {

      if (checkName(libAtom[resName], it->atomName, atomName) && !foundAtom) {

        foundAtom = true;
        Add_position *ptr_position;
        if (!b_solv)
          ptr_position = new Add_solute_position;
        else
          ptr_position = new Add_solvent_position;

        ptr_position->add(sys, molNum, atomNum, it->coord);
        delete ptr_position;
        atom.erase(it);
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

      bool is_H;
      if (!b_solv)
        is_H = sys.mol(molNum).topology().atom(atomNum).isH();
      else
        is_H = sys.sol(0).topology().atom(atomNum).isH();

      if (args.count("gch") < 0 || !is_H)
        warnNotFoundAtom(atomNum, atomName, resNum);
    }
  }
  void warnNotFoundAtom(int atomNum, const string &atomName, int resNum) const {

    cerr << "Warning: Could not find atom " << atomNum + 1 << " (" << atomName
         << "), in residue " << resNum + 1 << " (" << resName << ").\n"
         << "\tSet coordinates to (0.0 0.0 0.0)\n";
  }
  void warnIgnoredAtoms() const {

    for (unsigned int lineNum = 0; lineNum < atom.size(); lineNum++)

      cerr << "Warning: Ignored atom " << stripWhite(atom[lineNum].atomName)
           << " in residue " << resNum << " (" << resName << ").\n";
  }

  string resName;
  int resNum;
  vector<Atom> atom;
};

void set_positions(const Arguments &args, System &sys, LibAtom libAtom,
                   Residue &residue, int molNum, int resNum, int atomNum,
                   bool b_solv) {
  string atomName;
  if (!b_solv) {
    residue.resName = sys.mol(molNum).topology().resName(resNum);
    atomName = sys.mol(molNum).topology().atom(atomNum).name();
  } else {
    residue.resName = "SOLV";
    atomName = sys.sol(0).topology().atom(atomNum).name();
  }
  bool foundAtom = false;

  for (auto atom = residue.atom.begin(); atom < residue.atom.end(); ++atom) {

    if (checkName(libAtom[residue.resName], atom->atomName, atomName) &&
        !foundAtom) {

      foundAtom = true;
      Add_position *ptr_position;
      if (!b_solv)
        ptr_position = new Add_solute_position;
      else
        ptr_position = new Add_solvent_position;

      ptr_position->add(sys, molNum, atomNum, atom->coord);
      delete ptr_position;
      residue.atom.erase(atom);
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
    bool is_H;
    if (!b_solv)
      is_H = sys.mol(molNum).topology().atom(atomNum).isH();
    else
      is_H = sys.sol(0).topology().atom(atomNum).isH();

    if (args.count("gch") < 0 || !is_H)
      residue.warnNotFoundAtom(atomNum, atomName, resNum);
  }
}

class BFactor {
public:
  explicit BFactor(const string &filename) : bf_file{filename} {}
  ~BFactor() = default;
  void print(System &sys, Library library, list<Residue> residues) {
    if (!bf_file.is_open()) {
      throw gromos::Exception("pdb2g96",
                              "Cannot open @outbf file for writing.");
    }
    bf_file << "TITLE\nB-Factors and occupancies\n\nEND\nBFACTOROCCUPANCY\n";
    parse_molecules(sys, library, residues, false);
    parse_molecules(sys, library, residues, true);
    bf_file << "END\n";
    bf_file.close();
  }

private:
  void parse_molecules(System &sys, Library library, list<Residue> &residues,
                       bool b_solv) {

    int molMax;
    if (!b_solv)
      molMax = sys.numMolecules();
    else
      molMax = 1;

    for (int molNum = 0; molNum < molMax; molNum++) {
      // loop over all residues
      int firstAtomNum;
      int resMax;
      if (!b_solv) {
        firstAtomNum = 0;
        resMax = sys.mol(molNum).topology().numRes();
      } else {
        resMax = residues.size();
      }

      for (int resNum = 0; resNum != resMax; ++resNum) {

        string resName;
        if (!b_solv)
          resName = sys.mol(molNum).topology().resName(resNum);
        else
          resName = "SOLV";

        auto it = residues.begin();

        try {
          it->checkResidueName(library.get_residue(), resName);
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
        int lastAtomNum;
        if (!b_solv) {
          while (sys.mol(molNum).topology().resNum(firstAtomNum) != resNum)
            firstAtomNum++;

          lastAtomNum = firstAtomNum;
          while (lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
                 sys.mol(molNum).topology().resNum(lastAtomNum) == resNum)
            lastAtomNum++;
        } else {
          firstAtomNum = 0;
          lastAtomNum = sys.sol(0).topology().numAtoms();
        }

        /*
         * for every atom in the topology residue,
         * look for an atom in the pdb residue with the same name,
         * and import its coordinates. If we can't find one,
         * set the coords to 0,0,0 and issue a warning.
         */
        for (int atomNum = firstAtomNum; atomNum < lastAtomNum; atomNum++) {
          string atomName;
          if (!b_solv)
            atomName = sys.mol(molNum).topology().atom(atomNum).name();

          print_bfactor(sys, library.get_atom(), *it, molNum, resNum, atomNum,
                        b_solv);
        }
        residues.pop_front();
      }
    }
  }
  void print_bfactor(System &sys, LibAtom libAtom, Residue &residue, int molNum,
                     int resNum, int atomNum, bool b_solv) {
    string atomName;
    if (!b_solv) {
      residue.resName = sys.mol(molNum).topology().resName(resNum);
      atomName = sys.mol(molNum).topology().atom(atomNum).name();
    } else {
      residue.resName = "SOLV";
      atomName = sys.sol(0).topology().atom(atomNum).name();
    }
    bool foundAtom = false;

    for (auto atom = residue.atom.begin(); atom < residue.atom.end(); ++atom) {

      if (checkName(libAtom[residue.resName], atom->atomName, atomName) &&
          !foundAtom) {

        foundAtom = true;
        if (!b_solv) {
          bf_file << "# " << setw(5) << molNum + 1 << setw(5) << resNum + 1
                  << setw(5) << residue.resName << setw(5) << atomNum + 1
                  << setw(5) << atomName;
          bf_file << '\n';
          bf_file << setw(15) << atom->bfactor << setw(15) << atom->occupancy
                  << '\n';
        } else {
          bf_file << "# " << setw(5) << residue.resName << atomNum + 1;
          bf_file << '\n';
          bf_file << setw(15) << atom->bfactor << setw(15) << atom->occupancy
                  << '\n';
        }
      }
    }
    if (!foundAtom) {
      if (!b_solv) {
        bf_file << "# " << setw(5) << molNum + 1 << setw(5) << resNum + 1
                << setw(5) << residue.resName << setw(5) << atomNum + 1
                << setw(5) << atomName;
        bf_file << ": not found!";
        bf_file << '\n';
        bf_file << setw(15) << 0.01 << setw(15) << 0.0 << '\n';
      } else {
        bf_file << "# " << setw(5) << residue.resName << atomNum + 1;
        bf_file << ": not found!";
        bf_file << '\n';
        bf_file << setw(15) << 0.01 << setw(15) << 0.0 << '\n';
      }
    }
  }

  ofstream bf_file;
};

void parse_Atoms(const Arguments &args, System &sys, Library library,
                 list<Residue> &residues, bool b_solv) {

  int molMax;
  if (!b_solv)
    molMax = sys.numMolecules();
  else
    molMax = 1;

  for (int molNum = 0; molNum != molMax; ++molNum) {
    // loop over all residues
    int firstAtomNum;
    int resMax;
    if (!b_solv) {
      firstAtomNum = 0;
      resMax = sys.mol(molNum).topology().numRes();
    } else {
      resMax = residues.size();
    }

    for (int resNum = 0; resNum != resMax; ++resNum) {

      string resName;
      if (!b_solv)
        resName = sys.mol(molNum).topology().resName(resNum);
      else
        resName = "SOLV";

      auto it = residues.begin();
      // if the residues in the pdb and the topology are
      // not identical, skip this loop.
      try {
        it->checkResidueName(library.get_residue(), resName);
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
      int lastAtomNum;
      if (!b_solv) {
        while (sys.mol(molNum).topology().resNum(firstAtomNum) != resNum)
          firstAtomNum++;

        lastAtomNum = firstAtomNum;
        while (lastAtomNum < sys.mol(molNum).topology().numAtoms() &&
               sys.mol(molNum).topology().resNum(lastAtomNum) == resNum)
          lastAtomNum++;
      } else {
        firstAtomNum = 0;
        lastAtomNum = sys.sol(0).topology().numAtoms();
      }

      /*
       * for every atom in the topology residue,
       * look for an atom in the pdb residue with the same name,
       * and import its coordinates. If we can't find one,
       * set the coords to 0,0,0 and issue a warning.
       */
      for (int atomNum = firstAtomNum; atomNum < lastAtomNum; atomNum++) {
        /* it->Match_Atom(library.get_atom(), sys, args, molNum, resNum,
         * atomNum, */
        /*                b_solv); */
        set_positions(args, sys, library.get_atom(), *it, molNum, resNum,
                      atomNum, b_solv);
      }
      /* it->warnIgnoredAtoms(); */
      residues.pop_front();
      // print a warning for the pdb atoms that were ignored
    }
  }
}

void top2pdb(const Arguments &args, System &sys, Library library,
             list<Residue> residues) {

  // reserve memory for the coordinates
  // determine which are hydrogens based on the mass
  for (auto molNum = 0; molNum != sys.numMolecules(); ++molNum) {
    sys.mol(molNum).initPos();
    sys.mol(molNum).topology().setHmass(1.008);
  }
  sys.sol(0).topology().setHmass(1.008);

  parse_Atoms(args, sys, library, residues, false);
  parse_Atoms(args, sys, library, residues, true);
}

class Residues {
public:
  const list<Residue> &get_residues() const { return residues; }
  void add_residue(const Residue &residue) { residues.push_back(residue); }

private:
  list<Residue> residues;
};

Residues readCifAtoms(const string &filename) {
  ifstream cifFile(filename);
  if (!cifFile.good()) {
    throw gromos::Exception("Ginstream",
                            "Could not open file '" + filename + "'");
  }
  if (!cifFile.is_open()) {
    throw gromos::Exception("Ginstream",
                            "could not open file '" + filename + "'");
  }

  Residues residues;
  string line;
  bool foundLoop = false;
  vector<string> keywords;
  while (getline(cifFile, line)) {
    /* Look for a loop block */
    if (!foundLoop && line.find("loop_") != string::npos) {
      foundLoop = true;
      continue;
    }

    /* Look for atom lines */
    if (foundLoop && line.find("_atom_site.") != string::npos) {
      istringstream iss(line);
      string word;
      while (iss >> word) {
        if (word.find("_atom_site.") != string::npos) {
          size_t pos = word.find_last_of(".");
          if (pos != string::npos && pos + 1 < word.length()) {
            // Extract the substring after the last period
            string keyword = word.substr(pos + 1);
            // Store the keyword and its order of appearance
            keywords.push_back(keyword);
          }
        }
      }
      cout << line << endl;
    }
  }

  cifFile.close();
  return residues;
}

Residues readPdbAtoms(const string &filename) {
  ifstream pdbFile(filename);
  if (!pdbFile.good()) {
    throw gromos::Exception("Ginstream",
                            "Could not open file '" + filename + "'");
  }
  if (!pdbFile.is_open()) {
    throw gromos::Exception("Ginstream",
                            "could not open file '" + filename + "'");
  }

  int resNum = 0;
  InPDBLine inPdbLine;
  Residue residue;
  Residues residues;

  while (!pdbFile.eof()) {
    inPdbLine.read_line(pdbFile);
    if (inPdbLine.atom() == "ATOM" || inPdbLine.hetatm() == "HETATM") {

      // check if we're in a new residue
      if (inPdbLine.resnum() != resNum) {

        resNum = inPdbLine.resnum();

        // if we're not in the first residue
        if (!residue.atom.empty())
          residues.add_residue(std::move(residue));

        residue.resNum = inPdbLine.resnum();
        residue.resName = inPdbLine.resname();
        residue.atom.clear();
      }
      residue.add_atom(Atom(inPdbLine));
    }
  }
  // push the last residue
  residues.add_residue(std::move(residue));
  pdbFile.close();
  return residues;
}

string get_Extension(const string &filename) {
  // Find the position of the last occurrence of the period (.)
  size_t dotPos = filename.find_last_of('.');

  // If the period is found and it's not the last character in the string
  if (dotPos != string::npos && dotPos < filename.length() - 1) {
    // Return the substring starting from one position after the period
    return filename.substr(dotPos + 1);
  } else {
    // If no extension is found, return an empty string
    return "";
  }
}

void in_pdb_file(Arguments &args, System &sys) {

  /* read the library file */
  Library library;
  if (args.count("lib") > 0) {
    Ginstream lib(args["lib"]);
    cerr << "# using library file " << args["lib"] << endl;
    library.read_library(lib);
  }

  /* get the factor */
  if (args.count("factor") > 0) {
    double factor = stod(args["factor"]);
    Atom::fromang = 1.0 / factor;
  }

  Residues residues;
  if (args.count("pdb") > 0) {
    /* residues = readPdbAtoms(args["pdb"]); */

    /* Arguments::const_iterator it = args.lower_bound("outformat"), */
    /* to = args.upper_bound("outformat"); */
    string infile = args["pdb"];
    transform(infile.begin(), infile.end(), infile.begin(),
              [](unsigned char c) { return std::tolower(c); });
    string extension = get_Extension(infile);
    if (extension == "pdb") {
      residues = readPdbAtoms(args["pdb"]);
    } else if (extension == "cif") {
      residues = readCifAtoms(args["pdb"]);
    }
  }

  top2pdb(args, sys, library, residues.get_residues());

  /* Initiate writing of bfactors */
  if (args.count("outbf") > 0) {
    BFactor bff(args["outbf"]);
    bff.print(sys, library, residues.get_residues());
  }
}

void out_pdb2g96_file(Arguments &args, System &sys) {
  InTopology it(args["topo"]);
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

  try {
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    in_pdb_file(args, sys);
    // that's really it for pdb2g96
    out_pdb2g96_file(args, sys);
  } catch (const gromos::Exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }

  return 0;
}
