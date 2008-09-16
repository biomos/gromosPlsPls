/**
 * @file stacking.cc
 * Monitors the occurrence of stacking residues
 */

/**
 * @page programs Program Documentation
 *
 * @anchor stacking
 * @section stacking monitors the occurrence of stacking residues
 * @author @ref ns
 * @date 11-12-2007
 *
 * Program stacking monitors the occurrence of stacking residues over a molecular 
 * trajectory file through geometric criteria.
 *
 * Ring systems were considered to stack if the distance between
 * the centers of geometry of the rings is less than a given distance (typically
 * 0.5 nm) and the angle between the planes through the two rings is
 * less than a user specified angle (typically 30 degree).
 *
 * The user can specify two groups of atoms (donors and acceptors) between which
 * the stacking interactions are to be monitored. Ring systems are identified by
 * a list given in the library file. 
 * 
 * The program calculates average occurances and prints out a time series of the
 * observed residue pairs, distances and angles.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;time and dt&gt; </td></tr>
 * <tr><td> \@donor</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> \@acceptor</td><td>&lt;@ref AtomSpecifier&gt; </td></tr>
 * <tr><td> [\@paras</td><td>&lt;distance [nm] and angle [deg]; default: 0.5, 30&gt;] </td></tr>
 * <tr><td> \@library</td><td>&lt;stacking library file&gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  hbond
    @topo             ex.top
    @pbc              r
    @time             0 1
    @donor            1:a
    @acceptor         2:a
    @paras            0.5 30
    @library          ../data/stacking.lib
    @traj             ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <set>
#include <algorithm>
#include <iterator>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/StringTokenizer.h"
#include "../src/gcore/System.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/utils/AtomSpecifier.h"

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;
using namespace std;


struct Residue {
  System &sys;
  int mol, res;
  string name;
  AtomSpecifier atoms;
  AtomSpecifier planeAtoms;
  
  Residue(System& system, const int m, const int r, const string n) : 
    sys(system), mol(m), res(r), name(n) {
    atoms = AtomSpecifier(sys);
    planeAtoms = AtomSpecifier(sys);
  }
  Residue(const Residue& r) :
    sys(r.sys), mol(r.mol), res(r.res), name(r.name), atoms(r.atoms),
    planeAtoms(r.planeAtoms) { }
  
  Residue & operator=(const Residue &r) {
    mol = r.mol; res = r.res;
    name = r.name;
    atoms = r.atoms;
    planeAtoms = r.planeAtoms;
    return *this;
  }
  
  bool inline operator<(const Residue &r) const {
    return (mol == r.mol) ? res < r.res : mol < r.mol;
  }
  
  bool inline operator==(const Residue &r) const {
    return mol == r.mol && res == r.res;
  }
  
  friend ostream& operator<<(ostream &os,const Residue &obj);
  
  Vec getCentre() const;
  Vec getNorm() const;
  
};

Vec Residue::getCentre() const {
  Vec res(0.0);
  for(int i = 0; i < atoms.size(); i++) 
    res += atoms.pos(i);
  return res / atoms.size();
}

Vec Residue::getNorm() const {
  return (planeAtoms.pos(1) - planeAtoms.pos(0)).cross(planeAtoms.pos(2) - planeAtoms.pos(0));
}

ostream& operator<<(ostream &os,const Residue &obj) {
  os << obj.mol << ":" << obj.name << ":" << obj.res << "[";
  os.width(5);
  vector<string> strs = obj.atoms.toString();
  copy(strs.begin(), strs.end(), ostream_iterator<string>(os));
  os << "]";
  return os;
}

class StackingCalculator {
public:
  StackingCalculator(map<string, vector<string> > lib,
                     AtomSpecifier &donor, AtomSpecifier &acceptor,
                     ostream &o);
  
  void doFrame(double time);
  void inline setParameters(double dist, double ang) {
    distance = dist; angle = ang;
  }
  void printStatistics(ostream &os);
  
protected:
  map<string, vector<string> > library;
  ostream &out;
  vector<Residue> donor;
  vector<Residue> acceptor;
  
  map<pair<Residue, Residue>, unsigned int> index;
  map<pair<Residue, Residue>, unsigned int> counter;
  
  double distance, angle;
  int numFrames;
};

StackingCalculator::StackingCalculator(map<string, vector<string> > lib,
                   AtomSpecifier &d, AtomSpecifier &a,
                   ostream &o) : library(lib), out(o), distance(0.5), angle(30), numFrames(0) {
  out.precision(8);
  out << "#     time        ID      dist     angle" << endl;
  
  System &sys = *(d.sys());
  d.sort(); // sort the atoms.
  
  bool addAtom = false;
  int mol = -1, res = -1;
  vector<string> atoms;
  for(int i = 0; i < d.size(); i++) {
    // when a new residue begins
    if (mol != d.mol(i) || res != d.resnum(i)) {
      mol = d.mol(i); res = d.resnum(i);
      map<string, vector<string> >::iterator result = library.find(d.resname(i));
      // check whether in stacking lib
      if (result != library.end()) {
        atoms = result->second; // store the ring atoms
        donor.push_back(Residue(sys, mol, res, d.resname(i)));   
        addAtom = true; // add the atoms
      } else {
        addAtom = false; // don't add them
      }
    }
    
    vector<string>::iterator res = find(atoms.begin(), atoms.end(), d.name(i));
    if (addAtom && res != atoms.end()) { // only add the ring atoms
      donor.rbegin()->atoms.addAtom(d.mol(i), d.atom(i)); // add them to the last residue's atomspec
      if (std::distance(atoms.begin(), res) < 3)
        donor.rbegin()->planeAtoms.addAtom(d.mol(i), d.atom(i));
    }
  }
  
  // for debug
  // copy(donor.begin(), donor.end(), ostream_iterator<Residue>(cerr));

  if (donor.size() == 0) {
    throw gromos::Exception("stacking", "Could not identify any ring systems in "
            "the donor atom specifier given.");
  }  
  
  //just check wheter all residues have at least 3 atoms!
  vector<Residue>::const_iterator iter = donor.begin(), to = donor.end();
  for(; iter != to; iter++) 
    assert(iter->planeAtoms.size() == 3);
  
  sys = *(a.sys());
  a.sort();
  
  addAtom = false;
  mol = res = -1;
  atoms.clear();
  for(int i = 0; i < a.size(); i++) {
    if (mol != a.mol(i) || res != a.resnum(i)) {
      mol = a.mol(i); res = a.resnum(i);
      map<string, vector<string> >::iterator result = library.find(a.resname(i));
      if (result != library.end()) {
        atoms = result->second;
        acceptor.push_back(Residue(sys, mol, res, a.resname(i)));   
        addAtom = true;
      } else {
        addAtom = false;
      }
    }
      
    vector<string>::iterator res = find(atoms.begin(), atoms.end(), a.name(i));
    if (addAtom && res != atoms.end()) { // only add the ring atoms
      acceptor.rbegin()->atoms.addAtom(a.mol(i), a.atom(i)); // add them to the last residue's atomspec
      if (std::distance(atoms.begin(), res) < 3)
        acceptor.rbegin()->planeAtoms.addAtom(a.mol(i), a.atom(i));
    }
  }
  
  if (acceptor.size() == 0) {
    throw gromos::Exception("stacking", "Could not identify any ring systems in "
            "the acceptor atom specifier given.");
  }
  
  //just check wheter all residues have at least 3 atoms!
  iter = acceptor.begin(), to = acceptor.end(); 
  for(; iter != to; iter++) 
    assert(iter->planeAtoms.size() == 3);
  
  // copy(acceptor.begin(), acceptor.end(), ostream_iterator<Residue>(cerr));
  // cerr << "Having " << donor.size() << " donor residues and " << acceptor.size() << " acceptor residues." << endl;
}

void StackingCalculator::doFrame(double time) {
  numFrames++;
  // loop over all donor/acceptor pairs
  vector<Residue>::const_iterator i = donor.begin();
  vector<Residue>::const_iterator toDonor = donor.end();
  // store pairs that have been done to avoid to do them twice.
  set<pair<Residue, Residue> > done;
  
  for(; i != toDonor; i++) {
    vector<Residue>::const_iterator j = acceptor.begin();
    vector<Residue>::const_iterator toAcceptor = acceptor.end();
    for(; j != toAcceptor; j++) {
      if (*i == *j) continue; // residues are equal
      const pair<Residue, Residue> &p = (*i < *j) ? make_pair(*i, *j) : make_pair(*j, *i);
      
      // already done? Can happen if donor and acceptor overlap
      if (done.find(p) == done.end()) {
        done.insert(p);
      } else {
        continue;
      }
      
      // get the plane's centre and norm vectors
      Vec v1 = i->getCentre();
      Vec v2 = j->getCentre();
      Vec n1 = i->getNorm();
      Vec n2 = j->getNorm();
      double dist = (v1 - v2).abs();
      double ang = acos(n1.dot(n2) / (n1.abs() * n1.abs())) * 57.2958; // * 180 / PI
      if (dist <= distance && ang <= angle) {
        if (index.find(p) == index.end()) {
          index[p] = index.size(); // not size + 1 because size is already raised
          counter[p] = 0;
        }
        
        counter[p]++;
        out << setw(10) << time << setw(10) << index[p];
        out.setf(ios::floatfield, ios::fixed);
        out.precision(3);
        out << setw(10) << dist
            << setw(10) << ang << endl;
      }
    }
  }
}

void StackingCalculator::printStatistics(ostream &os) {
  os << "#  ID    res 1          res 2             #       %" << endl;
  
  vector<unsigned int> is; // copy map->second to a vector that can be sorted
  {
    map<pair<Residue, Residue>, unsigned int>::const_iterator it = index.begin();
    map<pair<Residue, Residue>, unsigned int>::const_iterator to = index.end();
    for(; it != to; it++) 
      is.push_back(it->second);
    sort(is.begin(), is.end());
  }
  
  vector<unsigned int>::const_iterator index_it = is.begin();
  vector<unsigned int>::const_iterator index_to = is.end();
  for(; index_it != index_to; index_it++) { // loop over vector
    map<pair<Residue, Residue>, unsigned int>::const_iterator it = index.begin();
    map<pair<Residue, Residue>, unsigned int>::const_iterator to = index.end();
    for(; it != to; it++) { // loop over residue list
      if (it->second == *index_it) { // found the index -> print
        const Residue &res1 = it->first.first;
        const Residue &res2 = it->first.second;
        os << setw(5) << it->second
           << setw(5) << res1.mol+1 << setw(5) << res1.name << setw(5) << res1.res+1
           << setw(5) << res2.mol+1 << setw(5) << res2.name << setw(5) << res2.res+1
           << setw(8) << counter[it->first];
        os.setf(ios::floatfield, ios::fixed);
        os.precision(2);        
        os << setw(8) << (100.0 * counter[it->first] / numFrames) << endl;
      }
    }
  } 
}


int main(int argc, char **argv){
  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "donor" << "acceptor" << "paras" 
         << "library" << "traj";
  
  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo           <molecular topology file>\n";
  usage += "\t@pbc            <boundary type> [<gathermethod>]\n";
  usage += "\t@time           <time and dt>\n";
  usage += "\t@donor          <atomspecifier>\n";
  usage += "\t@acceptor       <atomspecifier>\n";
  usage += "\t[@paras          <distance [nm] and angle; default: 0.5, 135>]\n";
  usage += "\t@library        <stacking library file\n";
  usage += "\t@traj           <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);
    
    InTopology it(args["topo"]);
    System sys(it.system());
    
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // GatherParser
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    double time=0, dt=1; 
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
	time=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("time"))
        dt=atof(iter->second.c_str());
    }
    
    
    map<string, vector<string> > library;
    Ginstream nf(args["library"]);
    vector<string> buffer;
    nf.getblock(buffer);
    if(buffer[0]!="STACKINGLIB")
      throw gromos::Exception("stacking",
			      "stacking library file does not contain a STACKINGLIB block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("stacking", "library file file " + nf.name() +
			      " is corrupted. No END in STACKINGLIB"
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    // read in the lib
    for(size_t i = 1; i < buffer.size() - 1; i++) {
      string residue;
      vector<string> atoms;
      istringstream ss(buffer[i]);
      ss >> residue;
      for(string atom; ss >> atom;)
        atoms.push_back(atom);
      library[residue] = atoms;
    }
    
    /* // to debug the library parsing
    map<string, vector<string> >::const_iterator iter = library.begin(), to = library.end();
    for(;iter != to; iter++) {
      cerr << iter->first << ": ";
      copy(iter->second.begin(), iter->second.end(), ostream_iterator<string>(cerr));
      cerr << endl;
    } */
    
    AtomSpecifier donor = AtomSpecifier(sys);
    {
      Arguments::const_iterator to = args.upper_bound("donor");
      for(Arguments::const_iterator iter = args.lower_bound("donor"); iter!=to;iter++)
        donor.addSpecifier(iter->second);
    }  
    if(donor.size()==0)
      throw gromos::Exception("stacking", "No donor-atoms specified!");
    
    AtomSpecifier acceptor = AtomSpecifier(sys);
    {
      Arguments::const_iterator to = args.upper_bound("acceptor");
      for(Arguments::const_iterator iter = args.lower_bound("acceptor"); iter!=to;iter++)
        acceptor.addSpecifier(iter->second);
    }  
    if(acceptor.size()==0)
      throw gromos::Exception("stacking", "No acceptor-atoms specified!");
    
    InG96 ic;
    ofstream timeseries("stackingts.dat");
    if(!timeseries.is_open()) {
      throw gromos::Exception("stacking", "Can't open time series file for writing!");
    }
    
    // create the handy calculator
    StackingCalculator sc(library, donor, acceptor, timeseries);
    
    if (args.count("paras") > 0) {
      if (args.count("paras") !=2) {
        throw gromos::Exception("stacking", "@para takes exactly two numbers: "
                "the distance and the angle threshold");
      }
      
      Arguments::const_iterator iter = args.lower_bound("paras");
      double distance = atof(iter->second.c_str()), angle = atof((++iter)->second.c_str());
      sc.setParameters(distance, angle);
    }
    
    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      ic.select("ALL");
      
      // loop over single trajectory
      
      while(!ic.eof()) {
	ic >> sys;
        // gather system
        (*pbc.*gathmethod)();
        
        // to a frame
        sc.doFrame(time);
        
        time += dt;
      }
      ic.close();
    }
    // yes. print the statistics -.-
    sc.printStatistics(cout);
    
  } catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

