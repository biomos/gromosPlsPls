#include <string>
#include <list>
#include <vector>
#include <fstream>
#include <sstream>
#include <cassert>
#include <iostream>


#include "../gromos/Exception.h"
#include "../gmath/Vec.h"
#include "InPDB.h"


using namespace std;

namespace gio {

  class InPDB_i {
  public:

    string filemane;

    bool readATOM;
    bool readHETATOM;

    vector<string> resSeq;

    vector<string> types;
    vector<int> serials;
    vector<string> atoms;
    vector<string> altLocs;
    vector<string> resNames;
    vector<int> seqNos;
    vector<double> X;
    vector<double> Y;
    vector<double> Z;
    vector<double> occupancies;
    vector<double> tempFactors;
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
      ssline << line;
      string type;
      ssline >> type;
      if ((type == "ATOM" && d_this->readATOM) ||
              (type == "HETATOM" && d_this->readHETATOM)) {

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

        ssline >> serial >> atom >> resName >> altLoc >> seqNo
                >> x >> y >> z >> occupancy >> tempFactor >> element;

        // error message if failure while reading
        if (ssline.bad() || ssline.fail()) {
          stringstream msg;
          msg << "bad line in " << d_this->filemane << " while reading (" << type << ")";
          throw InPDB::Exception(msg.str());
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
        d_this->elements.push_back(element);

        // add the residue to the sequence (in case of a new residue)
        if (res != seqNo) {
          res = seqNo;
          d_this->resSeq.push_back(resName);
        }

      }

      // clean the stringstream
      ssline.str("");
      ssline.clear();
    }

    // close the PDB file
    fin.close();

  }

  vector<string> InPDB::getResSeq() {
    return d_this->resSeq;
  }
  void InPDB::changeResSeq(unsigned int i, string newname){
    d_this->resSeq[i] = newname;
  }
  
  gmath::Vec InPDB::getAtomPos(unsigned int i){
    gmath::Vec pos(d_this->X[i],d_this->Y[i],d_this->Z[i]);
    return pos;
  }

  unsigned int InPDB::numAtoms() {
    return d_this->serials.size();
  }

  string InPDB::getResName(unsigned int i){
    return d_this->resNames[i];
  }

  string InPDB::getAtomName(unsigned int i){
    return d_this->atoms[i];
  }

}
