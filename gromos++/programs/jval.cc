//time series tser

#include <cassert>

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
// i will use properties
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
      
      for(; i>sys.mol(m).numAtoms(); m++)
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
  usage += "\t[@time  <t> <dt>] (optionol, only printing time series if given\n";
  
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
    

    // read in the j-value specifications
    
    Ginstream jf(args["jval"]);
    vector<string> buffer;
    jf.getblock(buffer);
    
    if(buffer[0]!="JVALRESSPEC")
      throw gromos::Exception("main","jval file does not contain an JVALRESSPEC block!");
    
    // kps all karplusses
    vector<karplus> kps;
    
    for(unsigned int jj=1; jj< buffer.size()-1; jj++){
      istringstream is(buffer[jj]);
      int i, j,k ,l;
      double fdum;
      is >> i >> j >> k >> l >> fdum;
      karplus kp(sys, i-1,j-1,k-1,l-1);
      //cout << i << " " << j << " " << k << " " << l << " " << fdum;
      
      is >> kp.j0 >> kp.delta >> kp.A >> kp.B >> kp.C;
      //cout << " " << kp.j0 << " " << kp.delta << " " << kp.A << " " << kp.B 
      //	   << " " << kp.C << " " << idum << endl;
      
      if(is.fail())
	throw gromos::Exception("jval", "Bad line in jval-file\n"+buffer[jj]);
      kps.push_back(kp);
    }
    jf.close();
    
    
    // now we have to create the properties
    PropertyContainer props(sys);
    for(unsigned int i=0; i< kps.size(); i++){
      ostringstream os;
      os << "t%" << kps[i].m_mol+1 << ":" 
	 << kps[i].m_i+1 << ","<< kps[i].m_j+1 << ","
	 << kps[i].m_k+1 << ","<< kps[i].m_l+1;
      props.addSpecifier(os.str());
    }

    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // define input coordinate
    InG96 ic;

    // define enough statistic classes
    vector<gmath::stat> stat;
    stat.resize(kps.size());

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
	  double cosphi=cos((props[i]->getValue()+kps[i].delta) * M_PI /180.0);
	  double J = kps[i].A * cosphi * cosphi +
	    kps[i].B * cosphi +
	    kps[i].C;
	  stat[i].addval(J);
	  if(do_time){
	    
	    cout << setw(5) << i 
		 << setw(12) << J
		 << endl;
	  }
	}
	time+=dt;
      }
    }
    if(do_time)
      cout << "#\n#\n# summary\n#\n";
    
    // print title
    cout << "#" 
	 << setw(4) << "num" 
	 << setw(5) << "mol"
	 << setw(10) << "residue"
	 << setw(19) << "atom names"
	 << setw(25) << "atom numbers"
	 << setw(11) << "A"
	 << setw(5) << "B"
	 << setw(5) << "C"
	 << setw(7) << "delta"
	 << setw(7) << "j0"
	 << setw(12) << "J ave"
	 << setw(12) << "rmsd"
	 << endl;

    // now print out and calculate overall performance
    double sum=0;
    double abssum=0;
    double ssum=0;
    
    for(unsigned int i=0; i< kps.size(); i++){
      int m= kps[i].m_mol;
      cout << setw(5) << i+1 
	   << setw(5) << m+1
	   << setw(5) << sys.mol(m).topology().resNum(kps[i].m_i)+1
	   << setw(5) << sys.mol(m).topology().resName(sys.mol(m).topology().resNum(kps[i].m_i))
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_i).name()
	   << "-"
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_j).name()
	   << "-"
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_k).name()
	   << "-"
	   << setw(5) << sys.mol(m).topology().atom(kps[i].m_l).name()
	   << setw(7) << kps[i].m_i+1
	   << "-"
	   << setw(5) << kps[i].m_j+1
	   << "-"
	   << setw(5) << kps[i].m_k+1
	   << "-"
	   << setw(5) << kps[i].m_l+1;
      

      cout << setw(7) << kps[i].A
	   << setw(5) << kps[i].B
	   << setw(5) << kps[i].C
	   << setw(7) << kps[i].delta
	   << setw(7) << kps[i].j0
	   << setw(12) << stat[i].ave()
	   << setw(12) << stat[i].rmsd()
	   << endl;
      sum+=kps[i].j0-stat[i].ave();
      abssum+=fabs(kps[i].j0-stat[i].ave());
      ssum+=(kps[i].j0-stat[i].ave())*(kps[i].j0-stat[i].ave());
    }
    cout << "\n#"
	 << setw(30) << "average deviation " << sum / kps.size() << endl;
    cout << "#"
	 << setw(30) << "average absolute deviation " 
	 << abssum / kps.size() << endl;
    cout << "#"
	 << setw(30) << "rmsd over deviations " 
	 << sqrt((ssum-sum*sum/kps.size())/kps.size()) << endl;
    
    
    ic.close();
  }
  catch (const gromos::Exception &e){
    cerr << "EXCEPTION:\t";
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



