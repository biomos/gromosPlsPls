/**
 *@file pairener.cc
 * Calculates non-bonded interaction energies (Van-der-Waals energies
 * and electrostatic energies) for all possible pairs of specified atoms
*/

/**
 *@page contrib Contrib Program Documentation
 *@anchor pairener
 * @section pairener calculates (non-bonded) interaction energies for all
 *          possible pairs of specified atoms
 * @author @ref co
 * @date 31-10-08
 *
 * The program pairener evaluates the Van-der-Waals and the electrostatic
 * interaction energies for each possible pair of selected atoms
 * within a given cut-off distance (default cut=1.4). The selection of atoms 
 * is done by the atomspecifier. In case of a pertubation with soft core 
 * atoms you have to specify them explicitly.
 * Pairener writes out the sum of the calculated Van-der-Waals as well as
 * of the calculated electrostatic interaction energies for each 
 * time step and the average of these interaction energies. 
 * If not explicitly specified kappa and al2 are
 * set to zero while eps is set to 62.0 (SPC water).  
 *
 *<b>arguments:</b>
 * <table border=0 cellpadding=0>  
 * <tr><td> \@topo</td><td>&lt;molecular topology file&gt; </td></tr>
 * <tr><td> \@pbc</td><td>&lt;periodic boundary conditions&gt; </td></tr>
 * <tr><td> \@atoms</td><td>&lt;@ref AtomSpecifier "atoms"&gt; </td></tr>
 * <tr><td> \@time</td><td>&lt;@ref utils::Time "time dt" &gt; </td></tr>
 * <tr><td> \@cut</td><td>&lt;cut-off distance &gt; </td></tr>
 * <tr><td> \@eps</td><td>&lt;epsilon for reaction field correction &gt; </td></tr>
 * <tr><td> \@kap</td><td>&lt;kappa for reaction field correction &gt; </td></tr>
 * <tr><td> \@soft</td><td>&lt;soft @ref AtomSpecifier "atoms" &gt; </td></tr>
 * <tr><td> \@al2</td><td>&lt;alpha * lambda ^2 for soft LJ atoms &gt; </td></tr>
 * <tr><td> \@traj</td><td>&lt;trajectory files &gt; </td></tr>
 * </table>
 *
 * Example:
 * @verbatim
     @topo     ex.top
     @pbc      r
     @atoms    1:1,2,3,4
     @time     0 5
     @cut      1.4
     @eps      2.4
     @kap      0
     @soft
     @al2      0
     @traj     ex.tr
   @endverbatim
 
 * <hr>
 */



#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/GromosForceField.h"
#include "../src/gcore/Exclusion.h"
#include "../src/gcore/LJType.h"
#include "../src/gcore/AtomPair.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;

  
int main(int argc, char **argv){

  Argument_List knowns;
  knowns << "topo" << "pbc" << "atoms" << "time" << "cut" << "eps" << "kap" 
         << "soft" << "al2" << "traj";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <atoms>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@cut <cut-off distance>\n";
  usage += "\t@eps <epsilon for reaction field correction>\n";
  usage += "\t@kap <kappa for reaction field correction>\n";
  usage += "\t@soft <soft atoms>\n";
  usage += "\t@al2 <alpha * lambda ^2 for soft LJ atoms>\n";
  usage += "\t@traj  <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, knowns, usage);

  //   get simulation time
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

  //  read topology
  InTopology it(args["topo"]);
  System sys(it.system());
  GromosForceField gff(it.forceField());

  //  set atoms
  AtomSpecifier atoms(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("atoms");
    Arguments::const_iterator to=args.upper_bound("atoms");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      atoms.addSpecifier(spec);
    }
  }
  int asize=atoms.size();

  // get soft atom list
  AtomSpecifier soft(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("soft");
    Arguments::const_iterator to=args.upper_bound("soft");
    for(;iter!=to;iter++){
      string spec=iter->second.c_str();
      soft.addSpecifier(spec);
    }
  }
  

  //   get cut-off distance
  double cut=1.4;
  {
    Arguments::const_iterator iter=args.lower_bound("cut");
    if(iter!=args.upper_bound("cut"))
      cut=atof(iter->second.c_str());
  }
  
  //  get epsilon
  double eps=62.0;
  {
    Arguments::const_iterator iter=args.lower_bound("eps");
    if(iter!=args.upper_bound("eps"))
      eps=atof(iter->second.c_str());
  }
  
  //  get kappa
  double kap=0;
  {
    Arguments::const_iterator iter=args.lower_bound("kap");
    if(iter!=args.upper_bound("kap"))
      kap=atof(iter->second.c_str());
  }
  double crf=((2-2*eps)*(1+kap*cut)-eps*(kap*kap*cut*cut)) / 
	     ((1+2*eps)*(1+kap*cut)+eps*(kap*kap*cut*cut));

  //  get al2
  double al2=0;
  {
    Arguments::const_iterator iter=args.lower_bound("al2");
    if(iter!=args.upper_bound("al2"))
      al2=atof(iter->second.c_str());
  }
  
  // print titles
  cout << "  Time";
  cout << "         vdw          el";
  cout << endl;

    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // open file for output
  ofstream fout("pair_tot.dat");
  
  // define input coordinate
  InG96 ic;

  // declare some variables
  double d1, d3, d6, qq, c12, c6, c126, drf, cut3=cut*cut*cut;
  int num_frames=0;
  Vec d;

  //   declare energies for averaging
  vector<double> vdw_m, el_m, vdw_s, el_s, vdw_tm, el_tm, vdw_ts, el_ts;
  // the last elements of these vectors contains the sum of all specified atoms
  for(int i=0;i<=asize;i++){
    vdw_tm.push_back(0.0); el_tm.push_back(0.0);
    vdw_m.push_back(0.0);  el_m.push_back(0.0);
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
    while(!ic.eof()){
      ic >> sys;
      //pbc->gather();
      vdw_m[asize]=0; el_m[asize]=0;
      
      // loop over the selected atoms
      for(int i=0;i<asize;i++){
	vdw_m[i]=0; el_m[i]=0;
	int mi=atoms.mol(i), ai=atoms.atom(i);
        int sft=0;
	for(int j=0;j<soft.size();j++)
          if(soft.mol(j) == mi && soft.atom(j) == ai) sft=1;
	
	// find out to which charge group this atom belongs
        Vec chgrp1(0.0,0.0,0.0);
        {
	  int begin=ai-1, end=ai;
          if(ai>0)
	    for(begin=ai-1;
                begin>=0&&sys.mol(mi).topology().atom(begin).chargeGroup()!=1;
                begin--);
	  for(end=ai;
              sys.mol(mi).topology().atom(end).chargeGroup()!=1;
              end++);
	  
	  //charge group goes from begin+1 to end
          for(int k=begin+1;k<=end;k++)
	    chgrp1+=sys.mol(mi).pos(k);
	  chgrp1/=(end-begin);
	}
	
	// loop over the remaining atoms
	for(int j=i+1;j<asize;j++){
          int mj=atoms.mol(j), aj=atoms.atom(j), sf=0;
  	  for(int k=0;k<soft.size();k++)
            if(soft.mol(k) == mj && soft.atom(k) == aj) sf=1;

	  // find out to which charge group this atom belongs
          Vec chgrp2(0.0,0.0,0.0);
          {
	    int begin=aj-1, end=aj;
            if(aj>0)
	      for(begin=aj-1;
                 begin>=0&&sys.mol(mj).topology().atom(begin).chargeGroup()!=1;
                 begin--);
	    for(end=aj;
                sys.mol(mj).topology().atom(end).chargeGroup()!=1;
                end++);
	  
	    //charge group goes from begin+1 to end
            for(int k=begin+1;k<=end;k++)
	      chgrp2+=sys.mol(mj).pos(k);
	    chgrp2/=(end-begin);
	  }
	  chgrp2=pbc->nearestImage(chgrp1, chgrp2, sys.box());
	    
          // check distance
	  d=chgrp2-chgrp1;
	  if(d.abs2()<=cut*cut){
            //determine exclusion / 1-4 or itself
            int nt=0, third=0;
  	    if(mj==mi){
              // aj is excluded from ai
	      for(int e=0; 
                  e<sys.mol(mi).topology().atom(ai).exclusion().size();e++){
	        if(aj==sys.mol(mi).topology().atom(ai).exclusion().atom(e))
                  nt=1;
	      }
	      // ai is excluded from aj
              if(!nt)
	        for(int e=0; 
                    e<sys.mol(mj).topology().atom(aj).exclusion().size();e++){
	  	  if(ai==sys.mol(mj).topology().atom(aj).exclusion().atom(e))
                    nt=1;
	        }
              // a is third neighbour of ai
              if(!nt){
                for(int e=0;
		    e<sys.mol(mi).topology().atom(ai).exclusion14().size();e++){
	          if(aj==sys.mol(mi).topology().atom(ai).exclusion14().atom(e))
                    third=1;
  	        }
	        // ai is third neighbour of a
                if(!third)
		  for(int e=0;
		      e<sys.mol(mj).topology().atom(aj).exclusion14().size();
                      e++){
		    if(ai==sys.mol(mj).topology().atom(aj).exclusion14().atom(e))
		      third=1;
		  }
	      }
	    }
	    // if not nt include this pair
            if(!nt) {
              LJType lj(gff.ljType(AtomPair(
                        sys.mol(mj).topology().atom(aj).iac(),
                        sys.mol(mi).topology().atom(ai).iac())));
              qq=sys.mol(mj).topology().atom(aj).charge()*
	         sys.mol(mi).topology().atom(ai).charge();
	      if(!third){ c12=lj.c12(); c6=lj.c6(); }
	      else{ c12=lj.cs12(); c6=lj.cs6(); }
              Vec dd=pbc->nearestImage(sys.mol(mi).pos(ai),
                                       sys.mol(mj).pos(aj),sys.box());
              d1=(sys.mol(mi).pos(ai) - dd).abs();
	      d3=d1*d1*d1;
	      d6=d3*d3;
	      if(sft||sf){
                if(c6!=0&&c12!=0) c126=c12/c6;
	        else c126=0;
                d6+=al2*c126;
	      }
              drf=1/d1 - 0.5*crf*d1*d1/(cut3) - (1-0.5*crf)/cut;

 	      vdw_m[i]+=(c12/d6-c6)/d6;
	      el_m[i]+=qq*drf*gff.fpepsi();
	    }
	  }
	  
	}
      }
      // calculate sums
      for(int i=0;i<asize;i++){
	vdw_tm[i]+=vdw_m[i];    el_tm[i]+=el_m[i];
        vdw_m[asize]+=vdw_m[i]; el_m[asize]+=el_m[i];
      }
      vdw_tm[asize]+=vdw_m[asize]; el_tm[asize]+=el_m[asize];
      
      //write output
      fout.precision(5);
      fout.setf(ios::right, ios::adjustfield);
      fout << setw(6) << time
           << setw(12) << vdw_m[asize] << setw(12) << el_m[asize]
	   << endl;
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << time
	   << setw(12) << vdw_m[asize] << setw(12) << el_m[asize]
	   << endl;
      time+=dt;
      num_frames++;
    }
  }
  fout << endl;

  fout << "  ave.";
  fout << setw(12) << vdw_tm[asize]/num_frames
       << setw(12) << el_tm[asize]/num_frames
       << endl;
  fout.close();
  
  cout << endl;

  cout << "  ave.";
  cout << setw(12) << vdw_tm[asize]/num_frames
       << setw(12) << el_tm[asize]/num_frames
       << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
