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
 * @f$\vec{r}_i,j,k@f$.
 *
 * @f[ C_1(t) = <\vec{r_i}(\tau) . \vec{r_i}(\tau+t)>_{\tau} @f]
 * @f[ C_2(t) = \frac{1}{2} ( 3 <\vec{r_i}(\tau) . \vec{r_i}(\tau+t)>^2 _{\tau} - 1 ) @f]
 *
 * Program rot_rel calculates the first and second order Legendre polynomials
 * and calculates the time correlation functions. The user specifies two of the
 * molecular axes, the third is defined as the cross product of the first two.
 * The program can average the correlation functions over multiple molecules in
 * the system.
 * Note that the output of this program can also be produced by a combination
 * of programs @ref tser and @ref tcf.
 *
 * <b>arguments:</b>
 * <table border=0 cellpadding=0>
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;boundary type&gt; </td></tr>
 * <tr><td>[\@time</td><td>&lt;@ref utils::Time "time and dt"&gt;]</td></tr>
 * <tr><td> \@ax1</td><td>&lt;VectorSpecifier: specify molecular axis 1&gt; </td></tr>
 * <tr><td> \@ax2</td><td>&lt;VectorSpecifier: specify molecular axis 2&gt; </td></tr>
 * <tr><td> \@average</td><td>(average over all molecules) </td></tr>
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
    @ax1   1:1,3
    @ax2   1:30,34
    @average
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
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/utils/Value.h"
#include "../src/utils/VectorSpecifier.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/Time.h"

using namespace fit;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;

int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "time" << "ax1" << "ax2" << "traj" << "average";

  string usage = argv[0];
  usage += "\n\t@topo    <molecular topology file>\n";
  usage += "\t@pbc     <boundary type>\n";
  usage += "\t[@time    <time and dt>]\n";
  usage += "\t@ax1     <VectorSpecifier: specify molecular axis 1>\n";
  usage += "\t@ax2     <VectorSpecifier: specify molecular axis 2>\n";
  usage += "\t@average (average over all molecules)\n";
  usage += "\t@traj    <trajectory files>\n";
  
 
  try{
    Arguments args(argc, argv, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    // get simulation time
    Time time(args);
  
    // parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    // read in the axes to follow
    // In this implementation, the user can specify two axes through a 
    // VectorSpecifier, which is either the vector connecting two atoms, or the
    // position of one atom.The third axis is defined as the cross-product of
    // these two.

    vector<utils::VectorSpecifier> vs1;
    vector<utils::VectorSpecifier> vs2;
    
    string s1=args["ax1"];
    string s2=args["ax2"];
    
    utils::VectorSpecifier vct1(sys,pbc);
    utils::VectorSpecifier vct2(sys,pbc);
    int nummol=1;
    
    // do we want to expand to all molecules?
    // This is a bit of an ugly hack
    if(args.count("average")>=0){

      nummol=sys.numMolecules();
      
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
    Vec v1, v2, v3;
    
    // define input coordinate
    InG96 ic;

    int numFrames=0;

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
	(*pbc.*gathmethod)();

	//loop over the molecules
	for(int m=0; m<nummol; m++){
	  //calculate the vectors for this molecule
          v1 = vs1[m]();
	  v2 = vs2[m]();
	  
	  //normalize these vectos
	  v1 = v1.normalize();
	  v2 = v2.normalize();
	  v3 = v1.cross(v2);
	  data[3*m  ].push_back(v1);
	  data[3*m+1].push_back(v2);
	  data[3*m+2].push_back(v3);
	}
	numFrames++;
      }
      ic.close();
    }
    double inp[3], ave1[3], ave2[3];
    
    // now calculate all the autocorrelation functions.
    for(int it=1; it< numFrames; it++){
      double sum1[3]={0.0,0.0,0.0};
      double sum2[3]={0.0,0.0,0.0};
      for(int j=0; j+it < numFrames; j++){
	for(int m=0; m < nummol; m++){
	  for(int k=0; k<3; k++){
	    inp[k]=data[3*m+k][j].dot(data[3*m+k][j+it]);
	    sum1[k]+=inp[k];
	    sum2[k]+=0.5*(3*inp[k]*inp[k]-1);
	  }
	}
      }
      // now print out the information
      cout << time;
      for(int k=0; k<3; k++){
        ave1[k]=sum1[k] / numFrames / nummol;
        ave2[k]=sum2[k] / numFrames / nummol;
	cout << setw(14) << ave1[k]
	     << setw(14) << ave2[k];
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

