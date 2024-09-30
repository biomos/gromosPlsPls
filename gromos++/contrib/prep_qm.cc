/**
 * @file prep_qm.cc
 * Converts coordinate files from pdb into QM zone specification (QMZONE block in qmmm file)
 */

/**
 * @page programs Program Documentation
 *
 * @anchor prep_qm
 * @section prep_qm Converts coordinate files from cnf into QM zone specification (QMZONE block in qmmm file)
 * @author @ref pb
 * @date 14-12-2018
*/

#include <string>
#include <list>
#include <vector>
#include <map>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutG96S.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/BondType.h"
#include "../src/gcore/Constraint.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/AngleType.h"
#include "../src/gcore/Dihedral.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/gmath/Vec.h"
#include "../src/gio/InBFactorOccupancy.h"
#include "../src/utils/Gch.h"


using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace args;
using namespace std;
using namespace utils;
using namespace bound;

std::string filter_digit(std::string instring){
  std::string s;

for (int i = 0; i < instring.size(); i++){
    if (isdigit(instring[i])) break;
     s+=instring[i];
}
  return s;
}

std::vector<Vec> readcoord(ifstream &cnffile,std::vector<std::string> &atoms){
  std::vector<Vec> coords;
  string incnfline;
  string dummy;
  string at;
  string at2;
  bool ended= false;
std::string line;
  while(!cnffile.eof()) {
    getline(cnffile, incnfline);
    line=incnfline.c_str();
    if (line.find("POSITION") != std::string::npos ) {
      while (ended == false) {
        Vec tvec;
        getline(cnffile, incnfline);
        line=incnfline.c_str();
        std::istringstream is(line);
        is >> dummy >> dummy >> at >> dummy >> tvec[0] >> tvec[1] >> tvec[2];
        at2=filter_digit(at);
        coords.push_back(tvec);
        atoms.push_back(at2);
        if (line.find("END") != std::string::npos ) {
          ended=true;
          break;
        }
      }
    }
  }
return coords;
}

int main(int argc, char **argv) {
  Argument_List knowns;
  knowns << "cnf1" << "cnf2" << "out";

  string usage = "# " + string(argv[0]);
  usage += "\t@cnf1   <first input coordinate file: gromos format, determined by extension>\n";
  usage += "\t@cnf2   <second (subsystem) input coordinate file: gromos format, determined by extension>\n";
  usage += "\t[@out  <resulting GROMOS qmmmm file> (optional, defaults to stdout)]\n";

    Arguments args(argc, argv, knowns, usage);

    std::string coordfile = args["cnf1"].c_str();
    std::string coordfile2 = args["cnf2"].c_str();
    std::vector<Vec>  full_qm;
    std::vector<Vec>  sub_qm;
    std::vector<unsigned int> IDQM;
    std::vector<std::string> ATQM;
    std::vector<std::string> ATQM2;
    // open and read pdb file
    ifstream cnffile1(args["cnf1"].c_str());
    ifstream cnffile2(args["cnf2"].c_str());
    std::map<std::string, unsigned int> at_to_num;
    at_to_num["H"]=1;
    at_to_num["C"]=6;
    at_to_num["N"]=7;
    at_to_num["O"]=8;
    at_to_num["F"]=9;
    at_to_num["P"]=15;
    at_to_num["S"]=16;
    at_to_num["Cl"]=17;
    at_to_num["Br"]=35;
    at_to_num["I"]=53;

    if (!cnffile1.good()) {
      throw gromos::Exception("Ginstream", "Could not open file '" + args["cnf1"] + "'");
    }
    if (!cnffile1.is_open()) {
      throw gromos::Exception("Ginstream", "could not open file '" + args["cnf1"] + "'");
    }
  //read QM-zone
    full_qm=readcoord(cnffile1,ATQM);
  //if second cnf file is given, we treat it as the QM-subsystem for qmmmm
    if (cnffile2.good()){
      sub_qm=readcoord(cnffile2,ATQM2);

    }
  std::cout<<"QMZONE" << std::endl;
  for (unsigned i =0;i<full_qm.size()-1;i++) {
    for (unsigned j = 0; j < sub_qm.size()-1; j++) {
        Vec diffV=full_qm[i]-sub_qm[j];
        if (diffV.abs() < 0.0001   ){
          if (  at_to_num[ATQM[i]]  != at_to_num[ATQM2[j]]    ) {//check if atoms are same, ie. if we cut through bond
            std::cout << "            " << ATQM2[i] << "    " << i << "  " << at_to_num[ATQM2[j]] << "  link me! "
                      << std::endl;
          }
          else{
            std::cout << "            "  << ATQM[i] << "    " << i << "  "  <<at_to_num[ATQM[i]] << "  0 " <<std::endl;
          }
          IDQM.push_back(i + 1);

        }
    }
  }
  std::cout<<"END" << std::endl;

    return 0;
  }
