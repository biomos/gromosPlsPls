/**
 * @file prep_noe.cc
 * Converts X-plor NOE data to GROMOS format
 */

#include <fstream>


/**
 * @page programs Program Documentation
 *
 * @anchor prep_noe
 * @section prep_noe Converts X-plor NOE data to GROMOS format
 * @author @ref mk @ref co
 * @date 11-8-2006
 *
 * Program prep_noe converts NOE data from an X-plor like format to GROMOS
 * format, determining the proper choice of pseudo- or virtual atoms based on
 * the topology and a library file. The output can be used to apply distance
 * restraints during a simulation using programs promd or md, or to analyse a
 * molecular trajectory using program @ref noe "noe". For a definition of the
 * different types of pseudo- and virtual atoms see volume 2, page XX. In cases
 * where the library-file specifies a stereospecific CH2 atom (type 4), but
 * does not indicate which of the two protons is specified, NOE upper bounds
 * are created for both protons. Program @ref post_noe can process the output
 * of an NOE analysis to determine the best assignment.
 *
 * The experimentally determined upper bounds are generally listed in a three
 * column format, with distances in Angstrom. Program prep_noe has three types of
 * parsing these three columns. 1) take the first value as the upper bound; 2)
 * take the sum of the first and third values as the upper bound (default); or
 * 3) take the difference between the first and second values (commonly the
 * lower bound).
 *
 * The experimentally determined upper bounds can be corrected for pseudo-atom
 * distances (addition of a geometric constant) or multiplicity factors
 * (typically multiplication with @f$N^{1/p}@f$, where N is the multiplicity of
 * indistinguishable protons involved and p is the averaging power). Such
 * corrections can either be applied to the distances or can be taken out of a
 * set of distances. 
 * 
 * The program can also write a filter file, which can be used to re-evaluate
 * a given analysis over a specific trajectory, without recalculating all
 * distances, through program @ref post_noe "post_noe".
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@title</td><td>&lt;NOE title for output&gt; </td></tr>
 * <tr><td> \@noe</td><td>&lt;X-plor like NOE specification file&gt; </td></tr>
 * <tr><td> \@lib</td><td>&lt;NOE specification library&gt; </td></tr>
 * <tr><td> [\@dish</td><td>&lt;carbon-hydrogen distance; default: 0.1 nm&gt;] </td></tr>
 * <tr><td> [\@disc</td><td>&lt;carbon-carbon distance; default: 0.153 nm&gt;] </td></tr>
 * <tr><td> [\@parsetype</td><td>&lt;Upper bound parse type: 1, 2 or 3&gt; ] </td></tr>
 * <tr><td> [\@correction</td><td>&lt;correction file&gt; [&lt;correction type&gt;] ] </td></tr>
 * <tr><td> [\@action</td><td>&lt;add&gt; or &lt;subtract&gt; correction from upper bound; default: add ] </td></tr>
 * <tr><td> [\@filter</td><td>&lt;discard NOE's above a certain distance [nm]; default 10000 nm&gt;] </td></tr>
 * <tr><td> [\@factor</td><td>&lt;conversion factor Ang to nm; default is 10&gt;] </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  prep_noe
    @topo          ex.top
    @title         octa-alanine
    @noe           noe.prep
    @lib           ../data/noelib.45a3
    @dish          0.1
    @disc          0.153
    @parsetype     2
    @correction    ../data/noecor.gromos96
    @action        add
    @filter        0.8
    @factor        10
 @endverbatim
 *
 * <hr>
 */

#include <cassert>
#include <vector>
#include <map>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <sstream>

#include <args/Arguments.h>
#include <gio/Ginstream.h>
#include <gio/InG96.h>
#include <gcore/System.h>
#include <gio/InTopology.h>
#include <gio/StringTokenizer.h>
#include <bound/Boundary.h>
#include <args/BoundaryParser.h>
#include <utils/VirtualAtom.h>
#include <utils/Neighbours.h>
#include <utils/Noe.h>
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace args;
using namespace gio;
using namespace bound;
using namespace std;
using namespace utils;


class Noeprep {
public:
  int residA;
  string atomA;
  int residB;
  string atomB;
  double dis;

  Noeprep(int resA, string A, int resB, string B, double d) {
    residA = resA;
    atomA = A;
    residB = resB;
    atomB = B;
    dis = d;
  }

  ~Noeprep() {
  }
};

int main(int argc, char *argv[]) {

  vector<int> vacogsubtypeA, vacogsubtypeB;

  // Usage string

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo         <molecular topology file>\n";
  usage += "\t@title        <NOE title for output>\n";
  usage += "\t@noe          <NOE specification file>\n";
  usage += "\t@lib          <NOE specification library>\n";
  usage += "\t[@dish        <carbon-hydrogen distance; default: 0.1 nm>]\n";
  usage += "\t[@disc        <carbon-carbon distance; default: 0.153 nm>]\n";
  usage += "\t[@parsetype   <Upper bound parse type: 1, 2 or 3> ]\n";
  usage += "\t        Choices are:\n";
  usage += "\t        1: Upper bound == first number\n";
  usage += "\t        2: Upper bound == first + third number (most common, default)\n";
  usage += "\t        3: Upper bound == first - second number (commonly the lower bound)\n";
  usage += "\t[@correction  <correction file> [<correction type>] ]\n";
  usage += "\t[@action      <add> or <subtract> correction from upper bound; default: add ]\n";
  usage += "\t[@filter      <discard NOE's above a certain distance [nm]; default 10000 nm>]\n";
  usage += "\t[@factor      <conversion factor Ang to nm; default is 10>]\n";

  // known arguments...
  Argument_List knowns;
  knowns << "topo" << "title" << "filter" << "factor" << "noe" << "lib"
          << "parsetype" << "correction" << "dish" << "disc" << "action";

  // prepare cout for formatted output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);
  cout.precision(3);

  try {

    // Getting arguments and checking if everything is known.
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    //read in title
    string tit;
    {
      Arguments::const_iterator iter = args.lower_bound("title");
      for (; iter != args.upper_bound("title"); ++iter) {
        tit += iter->second + " ";
      }
    }

    //read in filter threshold
    double filt = args.getValue<double>("filter", false, 10000.0);

    //read in conversion factor
    double conv = args.getValue<double>("factor", false, 10.0);

    // Read in and create the NOE list
    Ginstream nf(args["noe"]);
    vector<string> buffer;
    nf.getblock(buffer);
    if (buffer[0] != "NOESPEC")
      throw gromos::Exception("main",
            "NOESPEC file does not contain an NOESPEC block!");
    if (buffer[buffer.size() - 1].find("END") != 0)
      throw gromos::Exception("prep_noe", "NOE file " + nf.name() +
            " is corrupted. No END in NOESPEC"
            " block. Got\n"
            + buffer[buffer.size() - 1]);

    // in noe all noes will be stored.
    vector<Noeprep> noevec;

    // a map with noe's that need to be connected for postnoe
    map<int, vector<int> > connections;

    // get parsetype
    int ptype = args.getValue<int>("parsetype", false, 2);

    for (unsigned int j = 1; j < buffer.size() - 1; j++) {
      StringTokenizer tok(buffer[j], " \t");
      vector<string> tokens = tok.tokenize();
      // check if the input file format has at least 8 columns
      if (tokens.size() < 8) {
        throw gromos::Exception("prep_noe", "Too few columns in \""
                + args["noe"] + "\". Did you set the sequential numbers in column 1 "
                "(see manual for further information)?\n");
      }

      // tokens[0] would be the sequential NOE number not used here but nice
      // to have it in the file for comparision with the output of post_noe
      int a = atoi(tokens[1].c_str());
      int b = atoi(tokens[3].c_str());
      double d = atof(tokens[5].c_str());
      double e = atof(tokens[6].c_str());
      double f = atof(tokens[7].c_str());
      if (tokens.size() > 8) {
        unsigned int g = atoi(tokens[8].c_str());
        if (g != j)
          throw gromos::Exception("prep_noe",
                "Numbering in NOESPEC file (7th column) is not correct");
        int h = atoi(tokens[9].c_str());
        vector<int> links(h - 1);
        for (int ii = 0; ii < h - 1; ii++)
          links[ii] = atoi(tokens[10 + ii].c_str());
        connections[g - 1] = links;
      }
      // apply parse type
      switch (ptype) {
        case 1: d = d;
          break;
        case 2: d = d + f;
          break;
        case 3: d = d - e;
          break;
        default:
          throw gromos::Exception("prep_noe", args["parsetype"] +
                  " unknown. Known types are 1, 2 and 3");
      }

      noevec.push_back(Noeprep(a, tokens[2], b, tokens[4], d));
    }
    nf.close();

    // Read in and create the NOE library
    Ginstream nff(args["lib"]);
    buffer.clear();
    nff.getblock(buffer);
    
    vector<Noelib> noelib;
    parse_noelib(buffer, noelib);
    nff.close();

    //check whether to add or subtract correction
    bool add = true;
    bool sub = false;
    if (args.count("action") > 0) {
      if (args["action"] == "add" || args["action"] == "ADD") {
        add = true;
        sub = false;
      } else if (args["action"] == "sub" || args["action"] == "SUB") {
        add = false;
        sub = true;
      } else
        throw gromos::Exception("prep_noe",
              "action type " + args["action"] + " not known!");
    }

    //read in the correction file if it exists
    map<int, map<int, double> > pseudocorrectiondata; // first key: type
    map<int, map<int, double> > multiplicitydata; // second key: subtype
    bool pseudocorrection = false;
    bool multiplicitycorrection = false;
    try {
      args.check("correction");
      Ginstream corf(args["correction"]);
      //get NOECORGROMOS block
      buffer.clear();
      corf.getblock(buffer);

      if (buffer[0] != "NOECORGROMOS")
        throw gromos::Exception("main",
              "NOE correction file does not contain the "
              "NOECORGROMOS block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("prep_noe", "Correction file " + corf.name() +
              " is corrupted. No END in NOECORGROMOS"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 1; j < buffer.size() - 1; j++) {
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
        pseudocorrectiondata[atoi(tokens[0].c_str())][atoi(tokens[1].c_str())]
                = atof(tokens[2].c_str());
      }


      //get MULTIPLICITY block
      buffer.clear();
      corf.getblock(buffer);

      if (buffer[0] != "MULTIPLICITY")
        throw gromos::Exception("main",
              "NOE correction file does not contain the"
              " MULTIPLICITY block!");
      if (buffer[buffer.size() - 1].find("END") != 0)
        throw gromos::Exception("prep_noe", "Correction file " + corf.name() +
              " is corrupted. No END in MULTIPLICITY"
              " block. Got\n"
              + buffer[buffer.size() - 1]);

      for (unsigned int j = 1; j < buffer.size() - 1; j++) {
        StringTokenizer tok(buffer[j]);
        vector<string> tokens = tok.tokenize();
        multiplicitydata[atoi(tokens[0].c_str())][atoi(tokens[1].c_str())] = atof(tokens[2].c_str());
      }
      corf.close();

      //determine the correction type
      Arguments::const_iterator it = args.lower_bound("correction");
      ++it;

      string ctype;

      if (it == args.upper_bound("correction")) {
        pseudocorrection = true;
        multiplicitycorrection = true;
      } else ctype = it->second;

      if (ctype == "pseudo") pseudocorrection = true;
      else if (ctype == "multiplicity") multiplicitycorrection = true;
      else if (ctype != "")
        throw gromos::Exception("prep_noe",
              "Correction type " + ctype + " not known!" +
              "Use 'pseudo' or 'multiplicity' to apply " +
              "only one type of correction");

    }    catch (Arguments::Exception e) {
      cout << "# No correction file used!" << endl;
    }

    //try for disc and dish
    double dish = args.getValue<double>("dish", false, 0.1);
    double disc = args.getValue<double>("disc", false, 0.153);

    //cout the title and that kind of stuff
    cout << "TITLE" << endl;
    cout << "NOE specification file for: " << tit << endl;
    cout << "END" << endl;
    cout << "NOECALCSPEC" << endl;
    cout << "# DISH: carbon-hydrogen distance" << endl;
    cout << "# DISC: carbon-carbon distance" << endl;
    cout << "#" << setw(9) << "DISH" << setw(10) << "DISC" << endl;
    cout.precision(5);
    cout << setw(10) << dish << setw(10) << disc << endl;
    cout << "# IDR1, JDR1, KDR1, LDR1:\n";
    cout << "#       atom sequence numbers of the real atoms defining the geometric position\n";
    cout << "#       of the first atom of a distance restraint pair\n";
    cout << "# IDR2, JDR2, KDR2, LDR2:\n";
    cout << "#       atom sequence numbers of the real atoms defining the geometric position\n";
    cout << "#       of the second atom of a distance restraint pair\n";
    cout << "# ICDR: geometric code defining the positions of the atoms of a distance\n";
    cout << "#       restraint pair\n";
    cout << "# VACS: subtypes of virual atoms of type COG (-1). Possible subtypes are:\n";
    cout << "#         0: no subtype defined\n";
    cout << "#         1: aromatic flipping ring\n";
    cout << "#         2: non-stereospecific NH2 group\n";
    cout << "# R0:   upper bound\n";
    cout << "#" << setw(4) << "IDR1" << setw(5) << "JDR1"
            << setw(5) << "KDR1" << setw(5) << "LDR1"
            << setw(5) << "ICDR" << setw(5) << "VACS"
            << setw(5) << "IDR2" << setw(5) << "JDR2"
            << setw(5) << "KDR2" << setw(5) << "LDR2"
            << setw(5) << "ICDR" << setw(5) << "VACS"
            << setw(10) << "R0" << endl << "#" << endl;

    //open the filter file...
    ofstream filterfile;
    filterfile.open("noe.filter");
    filterfile << "TITLE" << endl;
    filterfile << "NOE filter file for: " << tit << endl;
    filterfile << "END" << endl;
    filterfile << "NOEFILTER" << endl;
    filterfile << "# noe"
            << setw(4) << "mol"
            << setw(11) << "residue"
            << setw(5) << "atom"
            << setw(5) << "atom"
            << setw(4) << "mol"
            << setw(11) << "residue"
            << setw(5) << "atom"
            << setw(5) << "atom"
            << setw(6) << "r0"
            << " filter noe" << endl;

    ofstream disresfile;
    disresfile.open("noe.dsr");
    disresfile << "TITLE" << endl;
    disresfile << "NOE distance restraints file for: " << tit << endl;
    disresfile << "END" << endl;
    disresfile << "DISTANCERESSPEC" << endl;
    disresfile << "#" << setw(9) << "DISH" << setw(10) << "DISC" << endl;
    disresfile.precision(5);
    disresfile << setw(10) << dish << setw(10) << disc << endl;
    disresfile << "#" << setw(4) << "IDR1" << setw(5) << "JDR1"
            << setw(5) << "KDR1" << setw(5) << "LDR1"
            << setw(5) << "ICDR" << setw(5) << "IDR2"
            << setw(5) << "JDR2" << setw(5) << "KDR2" << setw(5) << "LDR2"
            << setw(5) << "ICDR" << setw(10) << "R0" << setw(10) << "W0"
            << setw(10) << "NRAH" << endl << "#" << endl;

    //here goes the crappy code
    string resnameA, resnameB;
    int molA = 0, molB = 0, resnumA = 0, resnumB = 0, mol = 0, atA = 0, atB = 0, atNumA = 0, atNumB = 0, count = 0;
    int atNOE = 0, totalnoecount = 1;
    double filterbound = 0;

    for (int i = 0; i < int (noevec.size()); ++i) {

      vector<int> links;
      if (connections.count(i)) links = connections[i];

      atNOE = i;
      ostringstream atname;
      Noeprep NOE = noevec[i];
      //first find the residue-name corresponding to
      //your input-number in the topology
      mol = 0;
      atNumA = 0;
      while (NOE.residA > (atNumA += sys.mol(mol).topology().numRes())) {
        ++mol;
        if (mol > sys.numMolecules())
          throw gromos::Exception("prep_noe",
                "Residue number too high in input line:\n");
      }
      atA = NOE.residA;
      atA -= atNumA - sys.mol(mol).topology().numRes();
      resnumA = (atA - 1);
      molA = mol;
      resnameA = (sys.mol(mol).topology().resName(atA - 1));

      mol = 0;
      atNumB = 0;
      while (NOE.residB > (atNumB += sys.mol(mol).topology().numRes())) {
        ++mol;
        if (mol > sys.numMolecules())
          throw gromos::Exception("prep_noe", +"Residue number too high in input line:\n");
      }
      atB = NOE.residB;
      atB -= atNumB - sys.mol(mol).topology().numRes();
      resnumB = (atB - 1);
      molB = mol;
      resnameB = (sys.mol(mol).topology().resName(atB - 1));

      //then map the correct gromos-topology atomname
      //to your input-atomname based on the residue-name
      //and the input-atomname
      bool foundA = false;
      vector<VirtualAtom*> vatomA;
      vector<VirtualAtom*> vatomB;

      int p = 0;
      for (int k = 0; k< int (noelib.size()); ++k) {
        Noelib NOELIB = noelib[k];
        if (NOELIB.resname == resnameA && NOELIB.orgatomname == NOE.atomA) {
          //back to topology to get the atom number
          for (int f = 0; f < sys.mol(molA).numAtoms() && foundA == false; ++f) {
            if (sys.mol(molA).topology().atom(f).name() == NOELIB.gratomname &&
                    sys.mol(molA).topology().resNum(f) == resnumA) {
              int addA = 0;
              foundA = true;
              p = k;

              for (int i = 0; i < molA; ++i) addA += sys.mol(i).numAtoms();

              if (NOE.dis / conv > filt) cout << "#";
              vatomA = getvirtual(f + addA, NOELIB.NOETYPE, NOELIB.NOESUBTYPE, sys,
                      dish, disc);
              // remember the subtype if found VA type is COM (-1)
              if (NOELIB.NOETYPE == -1) {
                vacogsubtypeA.push_back(NOELIB.NOESUBTYPE);
              } else {
                vacogsubtypeA.push_back(0);
              }
            }
          }
        }
      }

      if (!foundA) {
        std::stringstream ss;
        string a;
        ss << NOE.residA;
        ss >> a;
        string b = NOE.atomA;
        string c = " ";
        string d = a + c + b;
        throw gromos::Exception("prep_noe ", d +
                " Noe specification not found in library!");
      }
      Noelib NA = noelib[p];


      bool foundB = false;
      for (int z = 0; z< int (noelib.size()); ++z) {
        Noelib NOELIBB = noelib[z];
        if (NOELIBB.resname == resnameB && NOELIBB.orgatomname == NOE.atomB) {
          //back to topology to get the atom number
          for (int g = 0; g < sys.mol(molB).numAtoms() && foundB == false; ++g) {
            if (sys.mol(molB).topology().atom(g).name() == NOELIBB.gratomname &&
                    sys.mol(molB).topology().resNum(g) == resnumB) {
              int addB = 0;
              foundB = true;
              count += 1;
              for (int i = 0; i < molB; ++i) addB += sys.mol(i).numAtoms();
              vatomB = getvirtual(g + addB, NOELIBB.NOETYPE, NOELIBB.NOESUBTYPE, sys,
                      dish, disc);
              // remember the subtype if found VA type is COM (-1)
              if (NOELIBB.NOETYPE == -1) {
                vacogsubtypeB.push_back(NOELIBB.NOESUBTYPE);
              } else {
                vacogsubtypeB.push_back(0);
              }

              atname << setw(3) << molA + 1 << " "
                      << setw(5) << resnumA + 1 << " "
                      << setw(4) << NA.resname << " "
                      << setw(4) << NA.gratomname << " "
                      << setw(4) << NA.orgatomname << " "
                      << setw(3) << molB + 1 << " "
                      << setw(5) << (resnumB + 1) << " "
                      << setw(4) << NOELIBB.resname << " "
                      << setw(4) << NOELIBB.gratomname << " "
                      << setw(4) << NOELIBB.orgatomname << " ";

            }
          }
        }
      }

      if (!foundB) {
        std::stringstream s;
        string aa;
        s << NOE.residB;
        s >> aa;
        string bb = NOE.atomB;
        string cc = " ";
        string dd = aa + cc + bb;
        throw gromos::Exception("prep_noe ", dd +
                " Noe specification not found in library!");
      }

      //spit out disresblock...
      int atomsA[4], atomsB[4];
      int creatednoe = 0;

      // print out "readable" constraints as a coment
      cout << "# " << count << atname.str() << endl;
      disresfile << "# " << count << atname.str() << endl;

      for (int va = 0; va < (int) vatomA.size(); ++va) {
        int offsetA = 1;
        VirtualAtom VA(*vatomA[va]);

        int mol = VA.conf().mol(0);
        for (int l = 0; l < mol; ++l) offsetA += sys.mol(l).numAtoms();

        for (int aa = 0; aa < 4; ++aa) {
          int att;
          if (VA.conf().size() > aa)
            att = VA.conf().atom(aa);
          else
            att = -1;

          atomsA[aa] = att;
        }

        for (int vb = 0; vb < (int) vatomB.size(); ++vb) {
          int offsetB = 1;
          VirtualAtom VB(*vatomB[vb]);
          int mol = VB.conf().mol(0);
          for (int l = 0; l < mol; ++l) offsetB += sys.mol(l).numAtoms();

          for (int bb = 0; bb < 4; ++bb) {
            int att;
            if (VB.conf().size() > bb)
              att = VB.conf().atom(bb);
            else
              att = -1;


            atomsB[bb] = att;
          }

          ostringstream ss, noeline, disresline;
          ss.setf(ios::right, ios::adjustfield);
          ss.setf(ios::fixed, ios::floatfield);
          ss.precision(5);

          for (int kk = 0; kk < 4; ++kk) {
            if (atomsA[kk] == -1) ss << setw(5) << 0;
            else ss << setw(5) << atomsA[kk] + offsetA;
          }
          ss << setw(5) << VA.type();
          disresline << ss.str();
          ss << setw(5) << vacogsubtypeA[i];
          noeline << ss.str();
          ss.str(""); // clear it

          for (int kk = 0; kk < 4; ++kk) {
            if (atomsB[kk] == -1) ss << setw(5) << 0;
            else ss << setw(5) << atomsB[kk] + offsetB;
          }
          ss << setw(5) << VB.type();
          disresline << ss.str();
          ss << setw(5) << vacogsubtypeB[i];
          noeline << ss.str();
          ss.str(""); // clear it

          double bound = NOE.dis / conv;
          //set up corrections
          vector<int> type;
          vector<int> stype;
          type.push_back(VA.type());
          type.push_back(VB.type());
          stype.push_back(vacogsubtypeA[i]);
          stype.push_back(vacogsubtypeB[i]);

          // first do the multiplicity correction and then the 
          // pseudo atom correction

          double mult = 1;
          double cor = 0;

          //check for multiplicity-correction              
          if (multiplicitycorrection) {
            for (int i = 0; i < (int) type.size(); ++i) {
              if (multiplicitydata.find(type[i])
                      != multiplicitydata.end()) {
                if (multiplicitydata.find(type[i])->second.find(stype[i])
                        != multiplicitydata.find(type[i])->second.end()) {
                  mult *= multiplicitydata.find(type[i])->second.find(stype[i])->second;
                }
              }
            }
          } //end if (multiplicitycorrection) 


          //check for gromos-correction
          if (pseudocorrection) {
            for (int i = 0; i < (int) type.size(); ++i) {
              if (pseudocorrectiondata.find(type[i])
                      != pseudocorrectiondata.end()) {
                if (pseudocorrectiondata.find(type[i])->second.find(stype[i])
                        != pseudocorrectiondata.find(type[i])->second.end()) {
                  cor += pseudocorrectiondata.find(type[i])->second.find(stype[i])->second;
                }
              }
            }
          } //end if (pseudocorrection)
          if (add) bound = bound * mult + cor;
          else if (sub) bound = (bound - cor) / mult;

          // in the filter file I also want the corrected bound
          filterbound = bound;

          disresline << ss.str();
          noeline << ss.str();

          cout << noeline.str() << setw(10) << bound << endl;
          disresfile << disresline.str() << setw(10) << bound
                  << setw(10) << 1.0 << setw(10) << 1 << endl; // half-harmonic attractive
          if (va > 0 || vb > 0) cout << "#automatically generated NOE distance to monitor!" << endl;

          ++creatednoe;
        }
      }

      //smack out the filterfile...
      for (int ii = 0; ii < creatednoe; ++ii) {
        filterfile << setw(5) << totalnoecount << " ";
        filterfile << atname.str();
        filterfile << setw(5) << filterbound << " ";
        int offset = totalnoecount - i + creatednoe - ii - 2;
        filterfile << " " << creatednoe + links.size() << " ";
        for (int iii = 0; iii < creatednoe; ++iii) {
          //                if (ii != iii) filterfile << " " << atNOE+1+iii;
          if (ii != iii) filterfile << " " << totalnoecount - ii + iii;
        }
        for (unsigned int iii = 0; iii < links.size(); ++iii) {
          int iiii = links.size() - 1 - iii;
          filterfile << " " << offset + links[iiii];
        }

        filterfile << endl;
        ++totalnoecount;
      }
    } //end for (int i=0; i < noevec.size()) ++i) ...

    cout << "END" << endl;
    filterfile << "END" << endl;
    disresfile << "END" << endl;
    filterfile.close();
  } catch (gromos::Exception e) {
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


