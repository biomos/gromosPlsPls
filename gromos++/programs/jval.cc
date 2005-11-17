/**
 * @file jval.cc
 * time series of J-values
 */

/**
 * @page programs Program Documentation
 *
 * @anchor jval
 * @section jval time series of J-values
 * @author @ref co @ref mc
 * @date 22. 11. 2004
 *
 * calculate the time series of J-values
 * 
 * arguments:
 * - topo topology
 * - pbc [v,r,t,c] [gathermethod]
 * - time t0 dt
 * - jval <jvalue specification file>
 * - traj trajectory
 * 
 * Example:
 * @verbatim
  jval
    @topo ex.top
    @pbc  r
    @time 0 0.1
    @jval ex.jval
    @traj ex.tr

    @endverbatim
 *
 * @bug Mar 22 2005: nearestImage calls in properties were missing
 * <hr>
 */

#include <cassert>
#include <sstream>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/fit/Reference.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gmath/Stat.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/utils/PropertyContainer.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;


class karplus{
 public:
  int m_mol;
  int m_i;
  int m_j;
  int m_k;
  int m_l;
  double j0;
  double delta;
  double A;
  double B;
  double C;
  karplus(const System &sys, int i, int j, int k, int l)
    {
      //determine molecule
      int m=0, offset=0;
      
      for(; i>=sys.mol(m).numAtoms()+offset; m++)
	offset += sys.mol(m).numAtoms();
      m_mol=m;
      
      m_i=i-offset;
      m_j=j-offset;
      m_k=k-offset;
      m_l=l-offset;
    }
  karplus(const karplus &k):
    m_mol(k.m_mol), m_i(k.m_i), m_j(k.m_j), m_k(k.m_k), m_l(k.m_l), j0(k.j0),
    delta(k.delta), A(k.A), B(k.B), C(k.C) {}
};


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "jval", "time", "traj"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@jval   <jvalue specifications>\n";
  usage += "\t[@time  <t> <dt>] (optional, only printing time series if given)\n";
  
  usage += "\t@traj   <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
    
    
    //   get simulation time
    double time=0, dt=1; 
    bool do_time=false;
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
	time=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("time")){
        dt=atof(iter->second.c_str());
	do_time=true;
      }
    }
    
    
    //  read topology
    args.check("topo",1);
    InTopology it(args["topo"]);
    
    System sys(it.system());
    
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // read in the j-value specifications
    
    Ginstream jf(args["jval"]);
    vector<string> buffer;
    jf.getblock(buffer);
    
    if(buffer[0]!="JVALRESSPEC")
      throw gromos::Exception("main","jval file does not contain an JVALRESSPEC block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw gromos::Exception("jval", "J-value file " + jf.name() +
			      " is corrupted. No END in "+buffer[0]+
			      " block. Got\n"
			      + buffer[buffer.size()-1]);
    // kps all karplusses
    vector<karplus> kps;
    
    for(unsigned int jj=1; jj< buffer.size()-1; jj++){
      istringstream is(buffer[jj]);
      int i, j,k ,l;
      double fdum;
      is >> i >> j >> k >> l >> fdum;
      karplus kp(sys, i-1,j-1,k-1,l-1);
      // cout << i << " " << j << " " << k << " " << l << " " << fdum;
      
      is >> kp.j0 >> kp.delta >> kp.A >> kp.B >> kp.C;
      //cout << kp.m_i << " " << kp.m_j << " " << kp.m_k << " " << kp.m_l << 
      //	" " << kp.j0 << " " << kp.delta << " " << kp.A << " " << kp.B 
      //   << " " << kp.C << " " << endl;
      
      if(is.fail())
	throw gromos::Exception("jval", "Bad line in jval-file\n"+buffer[jj]);
      kps.push_back(kp);
    }
    jf.close();
    
    
    // now we have to create the properties
    PropertyContainer props(sys, pbc);
    for(unsigned int i=0; i< kps.size(); i++){
      ostringstream os;
      os << "t%" << kps[i].m_mol+1 << ":" 
	 << kps[i].m_i+1 << ","<< kps[i].m_j+1 << ","
	 << kps[i].m_k+1 << ","<< kps[i].m_l+1;
      props.addSpecifier(os.str());
    }

    // define input coordinate
    InG96 ic;

    // define enough statistic classes
    vector<gmath::Stat<double> > stat;
    stat.resize(2*kps.size());

    // loop over all trajectories
    for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys;
	(*pbc.*gathmethod)();
      
	// calculate the props
	props.calc();

	if(do_time){
	  cout << "#\n#\n# TIME\t" << time << endl;	
	  cout << "#" << setw(4) << "num" 
	       << setw(12) << "j" << endl;
	}
	
	// from this calculate the j-value and store
	for(unsigned int i=0; i< kps.size(); i++){
	  double cosphi=cos((props[i]->getValue().scalar()+kps[i].delta) * M_PI /180.0);
	  double J = kps[i].A * cosphi * cosphi +
	    kps[i].B * cosphi +
	    kps[i].C;
	  stat[i].addval(J);
	  stat[kps.size()+i].addval(props[i]->getValue().scalar());
	  
	  if(do_time){
	    
	    cout << setw(5) << i 
		 << setw(12) << J
		 << endl;
	  }
	}
	time+=dt;
      }
      ic.close();
      
    }

    cout << "#\n#\n# JVALUE AVERAGES\n#\n"
	   << "# (print with a2ps --landscape --columns=1 --font-size=9)\n#\n";
    
    // print title
    cout << "#" 
	 << setw(4) << "num" 
	 << setw(4) << "mol"
	 << setw(9) << "residue"
	 << setw(15) << "atom names"
	 << setw(21) << "atom numbers"
	 << setw(11) << "A"
	 << setw(5) << "B"
	 << setw(5) << "C"
	 << setw(9) << "delta"
	 << setw(7) << "j0"
	 << setw(10) << "phi ave"
	 << setw(9) << "rmsd"
	 << setw(9) << "J ave"
	 << setw(8) << "rmsd"
	 << setw(7) << "dev"
	 << endl;

    // now print out and calculate overall performance
    double sum=0;
    double abssum=0;
    double ssum=0;
    
    for(unsigned int i=0; i< kps.size(); i++){
      int m= kps[i].m_mol;
      cout << setw(5) << i+1 
	   << setw(4) << m+1
	   << setw(4) << sys.mol(m).topology().resNum(kps[i].m_i)+1
	   << setw(5) << sys.mol(m).topology().resName(sys.mol(m).topology().resNum(kps[i].m_i))
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_i).name()
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_j).name()
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_k).name()
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_l).name()
	   << setw(7) << kps[i].m_i+1
	   << setw(5) << kps[i].m_j+1
	   << setw(5) << kps[i].m_k+1
	   << setw(5) << kps[i].m_l+1;
      
      cout.setf(ios::fixed, ios::floatfield);
      cout.precision(1);

      cout << setw(7) << kps[i].A
	   << setw(5) << kps[i].B
	   << setw(5) << kps[i].C;
      cout.precision(0);
      cout << setw(7) << kps[i].delta;
      cout.precision(2);
      cout << setw(7) << kps[i].j0;
      cout.precision(1);
      cout << setw(10) << stat[kps.size()+i].ave();
      cout.precision(2);
      cout << setw(9) << stat[kps.size()+i].rmsd()
	   << setw(9) << stat[i].ave();
      cout.precision(3);
      cout << setw(8) << stat[i].rmsd();
      cout.precision(1);
      cout << setw(7) << fabs(stat[i].ave() - kps[i].j0)
	   << endl;
      sum+=stat[i].ave()-kps[i].j0;
      abssum+=fabs(stat[i].ave()-kps[i].j0);
      ssum+=(kps[i].j0-stat[i].ave())*(kps[i].j0-stat[i].ave());
    }
    cout << "\n#"
	 << setw(30) << "average deviation " << sum / kps.size() << endl;
    cout << "#"
	 << setw(30) << "average absolute deviation " 
	 << abssum / kps.size() << endl;
    cout << "#"
	 << setw(30) << "root-mean-square deviation "
	 << sqrt(ssum/kps.size()) << endl;
    cout << "#"
	 << setw(30) << "rmsd over deviations " 
	 << sqrt((ssum-sum*sum/kps.size())/kps.size()) << endl;

  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



