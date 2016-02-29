#include <string>
#include <sstream>
#include <cassert>
#include <iostream>
#include <iomanip>
#include <set>
#include <algorithm>

#include "../gcore/System.h"
#include "../gromos/Exception.h"
#include "../args/Arguments.h"
#include "../args/OutformatParser.h"
#include "../gio/OutCoordinates.h"
#include "../gio/OutG96S.h"
#include "../gio/OutG96.h"
#include "../gio/OutPdb.h"
#include "../gio/Outvmdam.h"

using namespace std;
using namespace gcore;
using namespace args;
using namespace gio;

OutCoordinates * args::OutformatParser::parse(Arguments & args,
        string & ext, string argname) {
  OutCoordinates *oc;
  ext = ".cnf";
  if (args.count("outformat") > 0) {
    Arguments::const_iterator it = args.lower_bound("outformat"),
            to = args.upper_bound("outformat");
    string format = args["outformat"];
    transform(format.begin(), format.end(), format.begin(), static_cast<int (*)(int)> (std::tolower));
    if (format == "pdb") {
      ++it;
      if (it == to) {
        oc = new OutPdb();
      } else {
        double factor = 10.0;
        bool renumber=false;
        while (it != to) {
        istringstream is(it->second);
        if (is.str() == "renumber") { 
          renumber=true;
        }
        else if (!(is >> factor))
          throw gromos::Exception("OutformatParser", "@outformat pdb factor has to be numeric.!");
        ++it;
        }
        oc = new OutPdb(factor,renumber);
      }
      ext = ".pdb";
    } else if (format == "cnf") {
      oc = new OutG96S();
      ext = ".cnf";
    } else if (format == "trc") {
      oc = new OutG96();
      ext = ".trc";
    } else if (format == "por") {
      oc = new OutG96S(true);
      ext = ".por";
    } else if (format == "vmdam") {
      ++it;
      if (it == to) {
        oc = new Outvmdam();
      } else {
        istringstream is(it->second);
        double factor = 10.0;
        if (!(is >> factor))
          throw gromos::Exception("OutformatParser", "@outformat vmdam factor has to be numeric.!");
        oc = new Outvmdam(factor);
      }
      ext = ".vmd";
    } else {
      ostringstream msg;
      msg << "Output format '" << format << "' is unkown." << endl
              << "Known formats are:" << endl
              << "    - cnf" << endl
              << "      Configuration format containing the POSITION block." << endl
              << "    - trc" << endl
              << "      Trajectory format containing the POSITIONRED block." << endl
              << "    - por" << endl
              << "      Position restraints specification format." << endl
              << "    - pdb [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      Protein Data Bank (PDB) format." << endl
              << "    - vmdam [<factor to convert length unit to Angstrom, 10.0>]" << endl
              << "      VMD's Amber Coordinates format." << endl;

      throw gromos::Exception("OutformatParser", msg.str());
    }
  } else {
    oc = new OutG96S();
    ext = ".cnf";
  }

  return oc;
}


