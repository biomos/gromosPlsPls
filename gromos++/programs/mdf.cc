//minimum distance function
//chris

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/InG96.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gmath/Vec.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace std;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "time", "centre", "with", "nsm", "traj"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@centre type <type> or\n";
  usage += "\t        atom <atomspecifier> or\n";
  usage += "\t        cog  <atomspecifier> or\n";
  usage += "\t        all\n";
  usage += "\t@with   type <type> or\n";
  usage += "\t        atom <atomspecifier> or\n";
  usage += "\t        all\n";
  usage += "\t@nsm    <number of solvent molecules>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
try{
  Arguments args(argc, argv, nknowns, knowns, usage);

  // read topology
  args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());

  // read in number of solvent molecules
  int nsm=0;
  {
    Arguments::const_iterator iter=args.lower_bound("nsm");
    if(iter!=args.upper_bound("nsm"))
      nsm=atoi(iter->second.c_str());
  }

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

  // set centre atoms
  int sol_c=0;
  
  AtomSpecifier centre(sys);
  
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    int error=1;
    
    if(iter!=to){
      string s=iter->second.c_str();
      iter++;
      if(s=="type"){
        error=0;
        for(;iter!=to;iter++){
	  string name=iter->second.c_str();
	  for(int i=0;i<sys.numMolecules();i++)
            for(int j=0;j<sys.mol(i).topology().numAtoms();j++)
              if(name==sys.mol(i).topology().atom(j).name())
                centre.addAtom(i,j);
          for(int i=0;i<sys.sol(0).topology().numAtoms();i++)
            if(name==sys.sol(0).topology().atom(i).name())
              for(int j=0;j<nsm;j++){
		int off=j*sys.sol(0).topology().numAtoms();
		centre.addAtom(-1,i+off);
                sol_c++;
	      }
	  
	}
	
      }
      if(s=="atom"){
	error=0;
        for(;iter!=to;iter++){
	  string spec=iter->second.c_str();
	  centre.addSpecifier(spec);
	}
	
      }
      if(s=="cog"){
        error=0;
        AtomSpecifier temp(sys);
	centre.addAtom(-2,0);
	
	
	for(;iter!=to;iter++){
	  string spec=iter->second.c_str();
	  centre.addSpecifier(spec);
	}
      }
      if(s=="all"){
        error=0;
        for(int i=0;i<sys.numMolecules();i++)
	  for(int j=0;j<sys.mol(i).numAtoms();j++)
	    centre.addAtom(i,j);
      }
      if(error==1||centre.size()==0)
        throw gromos::Exception("rdf @centre", s + 
       " unknown or no atoms specified. Give 'type', 'atom', 'cog' or 'all'");
    }
  }
  // set atom to consider
  AtomSpecifier with(sys);
  int sol_w=0;
  
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    
    if(iter!=to){
      string s=iter->second.c_str();
      iter++;
      int error=1;
      
      if(s=="type"){
	error=0;
	
        for(;iter!=to;iter++){
	  string name=iter->second.c_str();
	  for(int i=0;i<sys.numMolecules();i++)
            for(int j=0;j<sys.mol(i).topology().numAtoms();j++)
              if(name==sys.mol(i).topology().atom(j).name())
                with.addAtom(i,j);
          for(int i=0;i<sys.sol(0).topology().numAtoms();i++)
            if(name==sys.sol(0).topology().atom(i).name())
              for(int j=0;j<nsm;j++){
		int off=j*sys.sol(0).topology().numAtoms();
		with.addAtom(-1,i+off);
		sol_w++;
	      }
	}
	
      }
      if(s=="atom"){
        error=0;
	
        for(;iter!=to;iter++){
	  string spec=iter->second.c_str();
	  with.addSpecifier(spec);
	}
	
      }
      if(s=="all"){
        error=0;
	
        for(int i=0;i<sys.numMolecules();i++)
	  for(int j=0;j<sys.mol(i).numAtoms();j++)
	    with.addAtom(i,j);
      }
      if(error==1||with.size()==0)
        throw gromos::Exception("rdf @with", s + 
	   " unknown or no atoms specified.\nGive 'type', 'atom' or 'all'." +
           "(is nsm=0 ?)");
    }
  }
  
  // Parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

  // define input coordinate
  InG96 ic;

  // loop over all trajectories
  int count_frame=0;
  
  // open centre.size() files
  ofstream fout[centre.size()];
  
  for(int i=0; i<centre.size(); i++){
    ostringstream os;
    os << "MIN_" << centre.mol(i)+1 << ":" << centre.atom(i)+1 
       << ".dat";
    fout[i].open(os.str().c_str());
    
  }
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    if(sol_c||sol_w) ic.select("ALL");
    
   
    // loop over single trajectory
    while(!ic.eof()){
      ic >> sys;
      
      if (nsm>sys.sol(0).numPos()/sys.sol(0).topology().numAtoms())
        throw gromos::Exception("rdf", 
	  " nsm specified is more than in coordinate file");  
      else
        sys.sol(0).setNumPos(nsm*sys.sol(0).topology().numAtoms());

      // loop over the centre atoms
      int start=0;
      int minat=0;
      
      Vec cog(0.0,0.0,0.0);
      
      if (centre.mol(0)==-2){
	start=centre.size()-1;
	for(int i=1;i<centre.size();i++)
          cog+=sys.mol(centre.mol(i)).pos(centre.atom(i));
	cog/=(centre.size()-1);
      }

      // now really loop over the centre atoms
      for(int i=start;i<centre.size();i++){
        Vec curr;
	if(centre.mol(i)>=0)
          curr=sys.mol(centre.mol(i)).pos(centre.atom(i));
        else
          curr=sys.sol(0).pos(centre.atom(i));
	double min=sys.box()[0];
	
        if(centre.mol(0)==-2) curr=cog;
        pbc->setReference(0,curr);
        pbc->gather();
	//loop over the atoms to consider
        for(int j=0;j<with.size();j++){
          //calculate distance only if this atom is not the current centre
          if(!(with.mol(j)==centre.mol(i)&&with.atom(j)==centre.atom(i))){
            Vec tmp=curr;
   	    if(with.mol(j)>=0)
              tmp-=sys.mol(with.mol(j)).pos(with.atom(j));
	    else
	      tmp-=sys.sol(0).pos(with.atom(j));
            if(tmp.abs()<min) {
              min=tmp.abs();
              minat=j;
	    }
	    
	  }
	}
	//write out min dist
        fout[i] << time << "\t" << min << "\t# " << with.mol(minat)+1 
                << ":" << with.atom(minat)+1 << endl;
	
      }
      count_frame++;
      time+=dt;
      
    }
    ic.close();
    
  }
  
  //close the files
  for(int i=0;i<centre.size();i++)
    fout[i].close();
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




























