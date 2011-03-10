#include <string>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>
#include <bits/basic_string.h>

#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
#include "InPDB.h"

// the pdb file format definitions (positions of data information)
// for more information see:
// http://www.rcsb.org/pdb/docs/format/pdbguide2.2/guide2.2_frame.html
#define TYPE_POS 0
#define TYPE_LEN 6
#define TYPE substr(TYPE_POS,TYPE_LEN)
#define SERIAL_POS 7
#define SERIAL_LEN 4
#define SERIAL substr(SERIAL_POS,SERIAL_LEN)
#define ATOMNAME_POS 13
#define ATOMNAME_LEN 4
#define ATOMNAME substr(ATOMNAME_POS, ATOMNAME_LEN)
#define RESNAME_POS 17
#define RESNAME_LEN 4
#define RESNAME substr(RESNAME_POS, RESNAME_LEN)
#define ALTLOC_POS 21
#define ALTLOC_LEN 1
#define ALTLOC substr(ALTLOC_POS,ALTLOC_LEN)
#define RESNUM_POS 22
#define RESNUM_LEN 4
#define RESNUM substr(RESNUM_POS, RESNUM_LEN)
#define COORDX_POS 30
#define COORDX_LEN 8
#define COORDX substr(COORDX_POS, COORDX_LEN)
#define COORDY_POS 38
#define COORDY_LEN 8
#define COORDY substr(COORDY_POS, COORDY_LEN)
#define COORDZ_POS 46
#define COORDZ_LEN 8
#define COORDZ substr(COORDZ_POS, COORDZ_LEN)
#define BFACTOR_POS 60
#define BFACTOR_LEN 6
#define BFACTOR substr(BFACTOR_POS,BFACTOR_LEN)
#define OCCUPANCY_POS 54
#define OCCUPANCY_LEN 6
#define OCCUPANCY substr(OCCUPANCY_POS,OCCUPANCY_LEN)
#define ELEMENT_POS 75
#define ELEMENT_LEN 3
#define ELEMENT substr(ELEMENT_POS,ELEMENT_LEN)

using namespace std;

namespace gio {

  /**
   * Class InPDB_ is the implementations class of the class PDB. It holds
   * all the member data so the PDB class needs nothing more than a pointer
   * to this class for its data.
   */
  class InPDB_i {
  public:

    /**
     * The file name of the pdb file which will be read.
     */
    string filemane;
    /**
     * A switch to read or not to read the ATOM atoms of the pdb file.
     */
    bool readATOM;
    /**
     * A switch to read or not to read the HETATOM atoms of the pdb file.
     */
    bool readHETATOM;
    /**
     * A vector to store the residue sequence as it is in the pdb file.
     */
    vector<string> resSeq;
    /**
     * A vector to store the atom types of the pdb atoms.
     */
    vector<string> types;
    /**
     * A vector to store the serial number of the pdb atoms.
     */
    vector<int> serials;
    /**
     * A vector to store the atom names of the pdb atoms.
     */
    vector<string> atoms;
    /**
     * A vector to store the chain number (for multimers) of the pdb atoms.
     */
    vector<string> altLocs;
    /**
     * A vector to store the residue names of the pdb atoms.
     */
    vector<string> resNames;
    /**
     * A vector to store the sequence number the pdb atoms are members of.
     */
    vector<int> seqNos;
    /**
     * A vector to store the x-coordinate of the pdb atoms.
     */
    vector<double> X;
    /**
     * A vector to store the y-coordinate of the pdb atoms.
     */
    vector<double> Y;
    /**
     * A vector to store the z-coordinate of the pdb atoms.
     */
    vector<double> Z;
    /**
     * A vector to store the occupancies of the pdb atoms.
     */
    vector<double> occupancies;
    /**
     * A vector to store the temp factors of the pdb atoms.
     */
    vector<double> tempFactors;
    /**
     * A vector to store the element symbols of the pdb atoms.
     */
    vector<string> elements;

  };

  // Constructor

  InPDB::InPDB(const std::string &filename, bool readATOM, bool readHETATOM) {
    d_this = new InPDB_i;
    d_this->filemane = filename;
    d_this->readATOM = readATOM;
    d_this->readHETATOM = readHETATOM;
  }

  InPDB::~InPDB() {
    if (d_this)delete d_this;
  }

  void InPDB::select(const std::string &thing) {
    if (thing == "ATOM") {
      d_this->readATOM = true;
      d_this->readHETATOM = false;
    } else if (thing == "HETATOM") {
      d_this->readATOM = false;
      d_this->readHETATOM = true;
    } else if (thing == "ALL") {
      d_this->readATOM = true;
      d_this->readHETATOM = true;
    } else {
      stringstream msg;
      msg << "InPDB::select does not know the argument " << thing;
      throw Exception(msg.str());
    }
  }

  void InPDB::read() {

    // open the PDB file
    ifstream fin(d_this->filemane.c_str());

    // check if the PDB file could be opened
    if (!fin.is_open()) {
      stringstream msg;
      msg << "Could not open the PDB file " << d_this->filemane;
      throw InPDB::Exception(msg.str());
    }

    // read and save the contents of the PDB file
    string line;
    stringstream ssline;
    int res = -1;
    while (!fin.eof()) {
      getline(fin, line);
      string type = "";
      if(line.size() > TYPE_POS) {
        ssline << line.TYPE << endl;
        ssline >> type;
      }
      ssline.clear();
      ssline.str("");
      if ((type == "ATOM" && d_this->readATOM) ||
              (type == "HETATOM" && d_this->readHETATOM)) {
        cerr << "1: line.size() = " << line.size() << ", BFACTOR_POS = " << BFACTOR_POS << endl;
        if(line.size() > BFACTOR_POS) { // at least GROMOS PDB file format
          ssline << line.SERIAL << endl;
          ssline << line.ATOMNAME << endl;
          ssline << line.RESNAME << endl;
          ssline << line.ALTLOC << endl;
          ssline << line.RESNUM << endl;
          ssline << line.COORDX << endl;
          ssline << line.COORDY << endl;
          ssline << line.COORDZ << endl;
          ssline << line.OCCUPANCY << endl;
          ssline << line.BFACTOR << endl;
          if(line.size() > ELEMENT_POS) { // officiela PDB file format
            ssline << line.ELEMENT << endl;
          }
          int serial;
          string atom;
          string altLoc;
          string resName;
          int seqNo;
          double x;
          double y;
          double z;
          double occupancy;
          double tempFactor;
          string element;
          // try to convert the read calues to the correct file type
          ssline >> serial >> atom >> resName >> altLoc >> seqNo
                    >> x >> y >> z >> occupancy >> tempFactor;
          // there is an additional column in the official PDF file format
          if(line.size() > ELEMENT_POS) {
            ssline >> element;
          }
          // error message in case of the conversion failed...
          if (ssline.bad() || ssline.fail()) {
            stringstream msg;
            msg << "bad line in PDB file:\n" << line;
            throw gromos::Exception("InPDB.cc", msg.str());
          }
          // memorize the read variables
          d_this->types.push_back(type);
          d_this->serials.push_back(serial);
          d_this->atoms.push_back(atom);
          d_this->altLocs.push_back(altLoc);
          d_this->resNames.push_back(resName);
          d_this->seqNos.push_back(seqNo);
          d_this->X.push_back(x);
          d_this->Y.push_back(y);
          d_this->Z.push_back(z);
          d_this->occupancies.push_back(occupancy);
          d_this->tempFactors.push_back(tempFactor);
          if (line.size() > ELEMENT_POS) {
            d_this->elements.push_back(element);
          }
          // add the residue to the sequence (in case of a new residue)
          if (res != seqNo) {
            res = seqNo;
            d_this->resSeq.push_back(resName);
          }
          ssline.clear();
          ssline.str("");
        } else {
          stringstream msg;
          msg << "bad line in PDB file:\n" << line;
          throw gromos::Exception("InPDB.cc", msg.str());
        }
      }
    } // while
    // close the PDB file
    fin.close();

  }

  vector<string> InPDB::getResSeq() {
    return d_this->resSeq;
  }

  /*void InPDB::changeResSeq(unsigned int i, string newname){
   *  d_this->resSeq[i] = newname;
  }*/


  gmath::Vec InPDB::getAtomPos(unsigned int i) {
    gmath::Vec pos(d_this->X[i], d_this->Y[i], d_this->Z[i]);
    return pos;
  }

  unsigned int InPDB::numAtoms() {
    return d_this->serials.size();
  }

  string InPDB::getResName(unsigned int i) {
    return d_this->resNames[i];
  }

  unsigned int InPDB::getResNumber(unsigned int i) {
    return d_this->seqNos[i];
  }

  string InPDB::getAtomName(unsigned int i) {
    return d_this->atoms[i];
  }

  string InPDB::getChain(unsigned int i) {
    return d_this->altLocs[i];
  }

  void InPDB::renumberRes() {
    int oldPDBnum;
    int newPDBnum;
    int resNum = 1;
    for (unsigned int i = 0; i < numAtoms() - 1; ++i) {
      oldPDBnum = d_this->seqNos[i];
      newPDBnum = d_this->seqNos[i + 1];
      d_this->seqNos[i] = resNum;
      if (oldPDBnum != newPDBnum) {
        resNum++;
      }

      //debug:
      //cout << oldPDBnum << "  " << newPDBnum << "  " << resNum << endl;

    }
    d_this->seqNos[numAtoms() - 1] = resNum;
  }

}
