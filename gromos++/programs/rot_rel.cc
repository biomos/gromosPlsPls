/**
 * @file rot_rel.cc
 * Calculates autocorrelation functions for rotational relaxation times
 */

/**
 * @page programs Program Documentation
 *
 * @anchor rot_rel
 * @section rot_rel Calculates autocorrelation functions for rotational relaxation times
 * @author @ref co
 * @date 25-06-07
 *
 * The rotational relaxation time of molecules can be estimated from the
 * autocorrelation function of the Legendre polynomials of molecular axes 
 * @f$\vec{r}_i,\vec{r}_j$ and $\vec{k}@f$.
 *
 * @f[ C_1(t) = \left<\vec{r_i}(\tau) \cdot \vec{r_i}(\tau+t)\right>_{\tau} @f]
 * @f[ C_2(t) = \frac{1}{2} ( 3 \left<\vec{r_i}(\tau) \cdot \vec{r_i}(\tau+t)\right>^2 _{\tau} - 1 ) @f]
 *
 * Program rot_rel calculates the first and second order Legendre polynomials
 * and calculates the time correlation functions. The user specifies two of the
 * molecular axes, the third is defined as the cross product of the first two.
 * The program can average the correlation functions over multiple molecules in
 * the system using the flags @average and @molecules. Note that in the flag @molecules,
 * molecules specified should be separated by spaces. Also, note that if @molecules is
 * not specified the program will average over all molecules except solvent.
 * Note that the output of this program can also be produced by a combination
 * of programs @ref tser and @ref tcf.
 * This program is parallelised.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td>[\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@ax1</td><td>&lt;@ref VectorSpecifier "vector" specifying molecular axis 1&gt; </td></tr>
 * <tr><td> \@ax2</td><td>&lt;@ref VectorSpecifier "vector" specifying molecular axis 2&gt; </td></tr>
 * <tr><td> [\@average</td><td>&lt;average over all molecules&gt;] </td></tr>
 * <tr><td> [\@molecules</td><td>&lt;specify molecules for averaging, separated by spaces] </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files&gt; </td></tr>
 * </table>
 *
 *
 * Example:
 * @verbatim
  rot_rel
    @topo  ex.top
    @pbc   r
    @time  0 1
    @ax1   atom(1:1,3)
    @ax2   atom(1:30,34)
    @average
    @molecules 1 2 3 4 5 6 7 8 9 10
    @traj ex.tr
 @endverbatim
 *
 * <hr>
 */


#include <cassert>
#include <vector>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/LJException.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Value.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/groTime.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "ax1" << "ax2" << "traj" << "average" << "molecules";

  string usage = argv[0];
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@ax1     <vector specifying molecular axis 1>\n";
  usage += "\t@ax2     <vector specifying molecular axis 2>\n";
  usage += "\t@average (average over all molecules)\n";
  usage += "\t@molecules <molecule numbers>(specify molecules for averaging, separated by spaces)\n";
  usage += "\t@traj    <trajectory files>\n";

  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    System refSys(it.system());

    // get simulation time
    Time time(args);
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(sys,refSys,args);

    // read in the axes to follow
    // In this implementation, the user can specify two axes through a 
    // VectorSpecifier, which is either the vector connecting two atoms, or the
    // position of one atom.The third axis is defined as the cross-product of
    // these two.

    vector<utils::VectorSpecifier> vs1;
    vector<utils::VectorSpecifier> vs2;
    
    string s1=args["ax1"];
    string s2=args["ax2"];

    vector<int> molecules;
    {
      Arguments::const_iterator to = args.upper_bound("molecules");
      for(Arguments::const_iterator iter = args.lower_bound("molecules"); iter != to; iter++){
        int tmp=atoi(iter->second.c_str());
        molecules.push_back(tmp);
      }
    }

    //for (int m=0; m<molecules.size(); ++m){
    //  cout << " " << molecules[m] << " " << endl;
    //}

    utils::VectorSpecifier vct1(sys,pbc);
    utils::VectorSpecifier vct2(sys,pbc);
    int nummol=1;

    // do we want to expand to all molecules?
    // This is a bit of an ugly hack
    // Bruno: Well, now I made it a bit uglier, but more flexible.
    // Let's say that the user wants to specify the molecules for
    // which the average has to be calculated.
    // Then, there will be two cases. If the flag @average is activated
    // without the activation of the flag @molecules, then all molecules
    // will be considered for the average, but if specific molecules are
    // defined, then the average runs only over those molecules.
    if(args.count("average")>=0){
      if(args.count("molecules") >= 0) {
        nummol = molecules.size();
        for(unsigned int m = 0; m < molecules.size(); ++m) {
          ostringstream t1, t2;
          t1 << s1.substr(0, s1.find("(") + 1) << molecules[m] << s1.substr(s1.find(":"), s1.size());
          t2 << s2.substr(0, s2.find("(") + 1) << molecules[m] << s2.substr(s2.find(":"), s2.size());
          vct1.setSpecifier(t1.str());
          vs1.push_back(vct1);
          vct2.setSpecifier(t2.str());
          vs2.push_back(vct2);
          
        }
      } 
      else {
        nummol = sys.numMolecules();
      
      for(int m=0; m < nummol; ++m){
	ostringstream t1, t2;
	t1 << s1.substr(0,s1.find("(")+1) << m+1 << s1.substr(s1.find(":"), s1.size());	
	t2 << s2.substr(0,s2.find("(")+1) << m+1 << s2.substr(s2.find(":"), s2.size());
	vct1.setSpecifier(t1.str());
	vs1.push_back(vct1);
	vct2.setSpecifier(t2.str());
	vs2.push_back(vct2);
      }
      }
    }
    else{
      vct1.setSpecifier(s1);
      vct2.setSpecifier(s2);
      vs1.push_back(vct1);
      vs2.push_back(vct2);
    }
    

    // print out a header
    cout << "# Calculating the autocorrelation function of the legendre"
	 << " polynomials" << endl;
    cout << "# of the dot-product of three molecular axes" << endl;
    cout << "# ax1 is defined as \"" << vs1[0].toString() << "\"\n";
    cout << "# ax2 is defined as \"" << vs2[0].toString() << "\"\n";
    cout << "# ax3 is the cross product of these two" << endl;
    cout << "#" << endl;
    if(nummol>1)
      cout << "# All vectors are calculated for " << vs1.size() << " molecules\n"
	   << "# autocorrelation functions are averaged over the molecules\n#\n";
    
    cout << "#" << setw(9) << "time"
	 << setw(14) << "<p_1(ax1)>"
	 << setw(14) << "<p_2(ax1)>"
	 << setw(14) << "<p_1(ax2)>"
	 << setw(14) << "<p_2(ax2)>"
	 << setw(14) << "<p_1(ax3)>"
	 << setw(14) << "<p_2(ax3)>"
	 << endl;
    
    // prepare vector to store all data.
    vector<vector<Vec> > data(nummol*3);
    // define input coordinate
    InG96 ic;

    int numFrames=0;
    vector<double> times;

    // loop over all trajectories
    for(Arguments::const_iterator 
	  iter=args.lower_bound("traj"),
	  to=args.upper_bound("traj");
	iter!=to; ++iter){
      
      // open file
      ic.open((iter->second).c_str());
      
      // loop over single trajectory
      while(!ic.eof()){
	ic >> sys >> time;
        times.push_back(time.time());
	(*pbc.*gathmethod)();

	//loop over the molecules
        #ifdef OMP
        #pragma omp parallel for
        #endif
	for(int m=0; m<nummol; m++){
	  //calculate the vectors for this molecule
          Vec v1 = vs1[m]();
	  Vec v2 = vs2[m]();
	  
	  //normalize these vectos
	  v1 = v1.normalize();
	  v2 = v2.normalize();
	  Vec v3 = v1.cross(v2);
	  v3 = v3.normalize();
	  data[3*m  ].push_back(v1);
	  data[3*m+1].push_back(v2);
	  data[3*m+2].push_back(v3);
	}
	numFrames++;
      }
      ic.close();
    }

    // Now we define a counter for the correct normalization
    int count = 0;

    // now calculate all the autocorrelation functions.
    for(int it=0; it< numFrames; it++){
      double frame_sum1[3]={0.0,0.0,0.0};
      double frame_sum2[3]={0.0,0.0,0.0};
      count = 0;
      #ifdef OMP
      #pragma omp parallel for
      #endif
      for(int j=0; j < numFrames-it; j++){
        double sum1[3]={0.0,0.0,0.0};
        double sum2[3]={0.0,0.0,0.0};
        #ifdef OMP
        #pragma omp critical
        #endif
        {
          count++;
        }
        for(int m=0; m < nummol; m++){
	  for(int k=0; k<3; k++){
	    const double inp = data[3*m+k][j].dot(data[3*m+k][j+it]);
	    sum1[k]+=inp;
	    sum2[k]+=0.5*(3.0*inp*inp-1.0);
           
	  }
	}
        #ifdef OMP
        #pragma omp critical
        #endif
        {
          for(int k=0; k<3; k++) {
            frame_sum1[k] += sum1[k];
            frame_sum2[k] += sum2[k];
          }
        }
      }
      // now print out the information
      cout << times[it];
      for(int k=0; k<3; k++){
        const double ave1 = frame_sum1[k] / count / nummol;
        const double ave2 = frame_sum2[k] / count / nummol;
	cout << setw(14) << ave1
	     << setw(14) << ave2;
      }
      cout << endl;
    }    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

