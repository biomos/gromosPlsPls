//ener calculates (non-bonded) interaction energies for specific atoms

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
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
#include <sstream>
#include <fstream>

using namespace std;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;

  
int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "atoms", "time", "cut", 
                    "eps", "kap", "soft", "al2", "traj"};
  int nknowns = 10;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@atoms <atomspecifier>\n";
  usage += "\t@time <time> <dt>\n";
  usage += "\t@cut <cut-off distance>\n";
  usage += "\t@eps <epsilon for reaction field correction>\n";
  usage += "\t@kap <kappa for reaction field correction>\n";
  usage += "\t@soft <atom specifier for soft atoms>\n";
  usage += "\t@al2 <alpha * lambda ^2 for soft LJ atoms>\n";
  usage += "\t@traj  <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

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
  cout<<"             solute                  solvent                 total"; 
  cout<< endl;
  cout << "  Time";
  for(int i=0;i<3;i++)
    cout << "         vdw          el";
  cout << endl;

    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);
  // parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

  // define input coordinate
  InG96 ic;

  //open asize files
  ofstream fout[asize];
  for(int i=0;i<asize;i++){
    ostringstream out;
    out << "ener_" << atoms.mol(i)+1 << ":" << atoms.atom(i)+1
        << ".dat";
    fout[i].open(out.str().c_str());
  }
  

  // declare some variables
  double d1, d3, d6, qq, c12, c6, c126, drf, cut3=cut*cut*cut;
  int num_frames=0;
  Vec d, dd;

  //   declare energies for averaging
  vector<double> vdw_m, el_m, vdw_s, el_s, vdw_tm, el_tm, vdw_ts, el_ts;
  // the last elements of these vectors contains the sum of all specified atoms
  for(int i=0;i<=asize;i++){
    vdw_tm.push_back(0.0); el_tm.push_back(0.0);
    vdw_ts.push_back(0.0); el_ts.push_back(0.0);
    vdw_m.push_back(0.0);  el_m.push_back(0.0);
    vdw_s.push_back(0.0);  el_s.push_back(0.0);
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
      // we gather once, to connect the molecules. Whenever distances
      // are calculated, we take the nearest image explicitly. This to
      // be sure and it might even be faster (you do it once for every
      // charge group and once for all atoms within the cut-off. With
      // large enough box-sizes this is less than doing it for every
      // atom.
      (*pbc.*gathmethod)();
      vdw_m[asize]=0; el_m[asize]=0;
      vdw_s[asize]=0; el_s[asize]=0;
      
      // loop over the selected atoms
      for(int i=0;i<asize;i++){
	vdw_m[i]=0; el_m[i]=0;
	vdw_s[i]=0; el_s[i]=0;
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
	//gather around this chargegroup
        //pbc->setReference(0,chgrp1);
	//pbc->gather();
	
	// loop over all the molecules, first solute
	for(int m=0;m<sys.numMolecules();m++){
          // loop over the charge groups
          int ch_b=0, ch_e=0;
          	  
	  while(ch_b<sys.mol(m).numAtoms()){
            Vec chgrp2(0.0,0.0,0.0);
	    for(ch_e=ch_b;
                sys.mol(m).topology().atom(ch_e).chargeGroup()!=1;
                ch_e++)
	      chgrp2+=sys.mol(m).pos(ch_e);
            chgrp2+=sys.mol(m).pos(ch_e);
	    chgrp2/=(ch_e-ch_b+1);
	    // do things with this charge group
            //cout << ch_b << " " << ch_e << endl;
	    chgrp2=pbc->nearestImage(chgrp1, chgrp2, sys.box());
	    
	    // check distance
	    d=chgrp2-chgrp1;
	    if(d.abs2()<=cut*cut){
              //cout << ch_b << " " << ch_e << endl;
	      for(int a=ch_b;a<=ch_e;a++){
                //determine exclusion / 1-4 or itself
		int nt=0, third=0, sf=0;
                if(!sft)
                  for(int k=0;k<soft.size();k++)
                    if(soft.mol(k)==m&&soft.atom(k)==a)
                      sf=1;
  	        if(m==mi){
		  if(a==ai) nt=1;
		  else{
                    // a is excluded from ai
		    for(int e=0; 
                        e<sys.mol(mi).topology().atom(ai).exclusion().size();
                        e++){
		      if(a==sys.mol(mi).topology().atom(ai).exclusion().atom(e))
                        nt=1;
		    }
		    // ai is excluded from a
                    if(!nt)
		      for(int e=0; 
                          e<sys.mol(m).topology().atom(a).exclusion().size();
                          e++){
		        if(ai==sys.mol(m).topology().atom(a).exclusion().atom(e))
                          nt=1;
		      }

                    // a is third neighbour of ai
                    if(!nt){
		      for(int e=0;
		          e<sys.mol(mi).topology().atom(ai).exclusion14().size();
		          e++){
		        if(a==sys.mol(mi).topology().atom(ai).exclusion14().atom(e))
                          third=1;
  		      }
		      // ai is third neighbour of a
                      if(!third)
		        for(int e=0;
		            e<sys.mol(mi).topology().atom(a).exclusion14().size();
		            e++){
		          if(ai==sys.mol(mi).topology().atom(a).exclusion14().atom(e))
		  	    third=1;
		        }
		    }
	 	  }
	        }
	        // if not nt include this pair
                if(!nt) {
                  LJType lj(gff.ljType(AtomPair(
                     sys.mol(m).topology().atom(a).iac(),
                     sys.mol(mi).topology().atom(ai).iac())));
                  qq=sys.mol(m).topology().atom(a).charge()*
		    sys.mol(mi).topology().atom(ai).charge();
		  if(!third){
		    c12=lj.c12();
		    c6=lj.c6();
	  	  }
		  else{
		    c12=lj.cs12();
		    c6=lj.cs6();
		  }
                  dd=pbc->nearestImage(sys.mol(mi).pos(ai), 
                                   sys.mol(m).pos(a), sys.box());
		  
                  d1=(sys.mol(mi).pos(ai)-dd).abs();
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
            ch_b=ch_e+1;
	  }
	}
        // now do the same for solvent. In GROMOS the solvent charge group is
	// centered on the first atom.
        int count=0;
	
        for(int a=0, na=sys.sol(0).topology().numAtoms(), 
            tna=sys.sol(0).numCoords();
            a<tna;a+=na){
          d=pbc->nearestImage(chgrp1,sys.sol(0).pos(a),sys.box())-chgrp1;
	  
	  if(d.abs2()<=cut*cut){
            count++;

            for(int at=0;at<na;at++){
              LJType lj(gff.ljType(AtomPair(
                sys.mol(mi).topology().atom(ai).iac(),
                sys.sol(0).topology().atom(at).iac())));
              qq=sys.mol(mi).topology().atom(ai).charge()*
	  	        sys.sol(0).topology().atom(at).charge();  
	      dd=pbc->nearestImage(sys.mol(mi).pos(ai), 
                          sys.sol(0).pos(a+at), sys.box());
	      
              d1=(sys.mol(mi).pos(ai)-dd).abs();
	      d3=d1*d1*d1;
	      d6=d3*d3;
              if(sft){
                if(lj.c6()!=0&&lj.c12()!=0) c126=lj.c12()/lj.c6();
                else c126=0;
                d6+=al2*c126;
	      }
              drf=1/d1 - 0.5*crf*d1*d1/(cut3) - (1-0.5*crf)/cut;
	      
	      vdw_s[i]+=(lj.c12()/d6-lj.c6())/d6;
	      el_s[i]+=qq*drf*gff.fpepsi();
	    }
	  }
	}
      }
      // calculate sums
      for(int i=0;i<asize;i++){
	vdw_tm[i]+=vdw_m[i];    el_tm[i]+=el_m[i];
        vdw_ts[i]+=vdw_s[i];    el_ts[i]+=el_s[i];
        vdw_m[asize]+=vdw_m[i]; el_m[asize]+=el_m[i];
        vdw_s[asize]+=vdw_s[i]; el_s[asize]+=el_s[i];
      }
      vdw_tm[asize]+=vdw_m[asize]; el_tm[asize]+=el_m[asize];
      vdw_ts[asize]+=vdw_s[asize]; el_ts[asize]+=el_s[asize];
      
      //write output
      for(int i=0;i<asize;i++){
	fout[i].precision(5);
	fout[i].setf(ios::right, ios::adjustfield);
	fout[i] << setw(6) << time
	        << setw(12) << vdw_m[i] << setw(12) << el_m[i]
                << setw(12) << vdw_s[i] << setw(12) << el_s[i]
                << setw(12) << vdw_m[i]+vdw_s[i]
                << setw(12) << el_m[i]+el_s[i]
	        << endl;
      }
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << time
	   << setw(12) << vdw_m[asize] << setw(12) << el_m[asize]
           << setw(12) << vdw_s[asize] << setw(12) << el_s[asize]
           << setw(12) << vdw_m[asize]+vdw_s[asize]
           << setw(12) << el_m[asize]+el_s[asize]
	   << endl;
      time+=dt;
      num_frames++;
    }
  }
  for(int i=0;i<asize;i++){
    fout[i] << endl;
    fout[i] << "  ave."
            << setw(12) << vdw_tm[i]/num_frames 
            << setw(12) << el_tm[i]/num_frames
            << setw(12) << vdw_ts[i]/num_frames
            << setw(12) << el_ts[i]/num_frames
            << setw(12) << (vdw_tm[i]+vdw_ts[i])/num_frames
            << setw(12) << (el_tm[i]+el_ts[i])/num_frames

            << endl;
    fout[i].close();
  }
  

  cout << endl;

  cout << "  ave.";
  cout << setw(12) << vdw_tm[asize]/num_frames
       << setw(12) << el_tm[asize]/num_frames
       << setw(12) << vdw_ts[asize]/num_frames
       << setw(12) << el_ts[asize]/num_frames
       << setw(12) << (vdw_tm[asize]+vdw_ts[asize])/num_frames
       << setw(12) << (el_tm[asize]+el_ts[asize])/num_frames
       << endl;
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}







