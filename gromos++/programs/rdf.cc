//radial distribution function
//chris
#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/Reference.h"
#include "../src/fit/PositionUtils.h"
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
#include "../src/gmath/Distribution.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>

using namespace fit;
using namespace gcore;
using namespace gmath;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "centre", "with", "cut", "grid", "nsm", "traj"};
  int nknowns = 8;

  string usage = argv[0];
  usage += "\n\t@topo   <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@centre type <type> or\n";
  usage += "\t        atom <atomspecifier> or\n";
  usage += "\t        cog  <atomspecifier> or\n";
  usage += "\t        all\n";
  usage += "\t@with   type <type> or\n";
  usage += "\t        atom <atomspecifier> or\n";
  usage += "\t        all\n";
  usage += "\t@nsm    <number of solvent molecules>\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
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

  // read in cut-off distance
  double cut=1.0;
  {
    Arguments::const_iterator iter=args.lower_bound("cut");
    if(iter!=args.upper_bound("cut"))
      cut=atof(iter->second.c_str());
  }
  
  // read in grid number
  int grid=100;
  {
    Arguments::const_iterator iter=args.lower_bound("grid");
    if(iter!=args.upper_bound("grid"))
      grid=atoi(iter->second.c_str());
  }
  
  // Parse boundary conditions
  double vol_corr=1;
  
  Boundary *pbc;
  try{
    char b=args["pbc"].c_str()[0];
    switch(b){
      case 't':
        pbc=new TruncOct(&sys);
        vol_corr=0.5;
        break;
      case 'v':
        pbc=new Vacuum(&sys);
        break;
      case 'r':
        pbc=new RectBox(&sys);
        break;
      default:
        throw gromos::Exception("Boundary", args["pbc"] + 
				" unknown. Known boundaries are t, r and v");
	
    }
  }
  catch(Arguments::Exception &e){
    pbc = new Vacuum(&sys);
  }

  //parse gather method
  Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

  /*
  int start=0;
  
  if(centre.mol(0)==-2) {
    cout << "centre of geometry for atoms:\n";
    start=1;
  }
  
  for(int i=start; i<centre.size();i++)
    cout << centre.mol(i) << ":" << centre.atom(i) << endl;
  
  cout << "considering atoms\n";
  
  for(int i=0; i<with.size();i++)
    cout << with.mol(i) << ":" << with.atom(i) << endl;
  */

  // define input coordinate
  InG96 ic;

  // set up distribution arrays

  double rdf[grid];
  double correct=4*acos(-1.0)*cut/double(grid);
  double vol,dens, r;
  
  for(int i=0;i<grid;i++) rdf[i]=0;
  
  // loop over all trajectories
  int count_frame=0;
  
  for(Arguments::const_iterator 
        iter=args.lower_bound("traj"),
        to=args.upper_bound("traj");
      iter!=to; ++iter){

    // open file
    ic.open((iter->second).c_str());
    if(sol_c||sol_w) ic.select("ALL");
    
   
    // loop over single trajectory
    while(!ic.eof()){
      gmath::Distribution dist(0, cut, grid);
      ic >> sys;
      
      if (nsm>sys.sol(0).numCoords()/sys.sol(0).topology().numAtoms())
        throw gromos::Exception("rdf", 
	  " nsm specified is more than in coordinate file");  
      else
        sys.sol(0).setnumCoords(nsm*sys.sol(0).topology().numAtoms());

      //pbc->gather();
      // loop over the centre atoms
      int start=0;
      
      Vec cog(0.0,0.0,0.0);
      
      if (centre.mol(0)==-2){
	start=centre.size()-1;
	for(int i=1;i<centre.size();i++)
          cog+=sys.mol(centre.mol(i)).pos(centre.atom(i));
	cog/=(centre.size()-1);
      }
      //calculate the volume
      vol=sys.box()[0]*sys.box()[1]*sys.box()[2]*vol_corr;

      // now really loop over the centre atoms
      for(int i=start;i<centre.size();i++){
        Vec curr;
	if(centre.mol(i)>=0)
          curr=sys.mol(centre.mol(i)).pos(centre.atom(i));
        else
          curr=sys.sol(0).pos(centre.atom(i));
	// see if this atom is also in the with-list
        int inwith=0;
	
        for(int j=0;j<with.size();j++)
	  if(with.atom(j)==centre.atom(i)&&with.mol(j)==centre.atom(j))
            inwith=1;
	  
        if(centre.mol(0)==-2) curr=cog;
        pbc->setReference(0,curr);
        (*pbc.*gathmethod)();
	//loop over the atoms to consider
        for(int j=0;j<with.size();j++){
          //calculate distance only if this atom is not the current centre
          if(!(with.mol(j)==centre.mol(i)&&with.atom(j)==centre.atom(i))){
            Vec tmp=curr;
   	    if(with.mol(j)>=0)
              tmp-=sys.mol(with.mol(j)).pos(with.atom(j));
	    else
	      tmp-=sys.sol(0).pos(with.atom(j));
	    dist.add(tmp.abs());
	  }
	}
	// now calculate the g(r) for this atom
        dens=(with.size()-inwith)/vol;
        for(int i=0; i<grid;i++){
          r=dist.value(i);
	  rdf[i]+=double(dist[i])/(dens*correct*r*r);
	}
	
	  
      }
      count_frame++;
      
    }

  }
  ic.close();
  //now correct the distribution for the number of frames and the number 
  //of centre atoms
  cout << "number of frames considered: " << count_frame << endl;
  int divide=1;
  
  if (centre.mol(0)==-2) divide*=count_frame;
  else divide*=count_frame*centre.size();
  for(int i=0;i<grid;i++){
    double r=(double(i)+0.5)*cut/grid;
    cout << r << "\t" << rdf[i]/double(divide) << endl;
  }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

