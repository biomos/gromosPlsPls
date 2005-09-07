/**
 * @file sdop.cc
 * @page programs Program Documentation
 *
 * @anchor sdop
 * @section solvent density and orientation persistence
 * @author @ref dk/mc
 * @date 25. 08. 2005
 *
 * calculate solvent density and orientation persistence
 * 
 * arguments:
 * - topo topologies
 * - pbc [v,r,t,c] [gathermethod]
 * - atomsfit [atomspecifier that determines the atoms to be used for the fitting]
 * - ref [reference coordinates]
 * - centre [atomspecifier]
 * - [@centrecog] take cog for centre atoms
 * - with [atomspecifier]
 * - cut [maximum distance]
 * - grid [number of point]
 * - traj
 *
 * Example:
 * @verbatim
 sdop 
 @topo TOL.top
 @pbc t
 @atomsfit a:1-12
 @ref refTOL.g96
 @centre a:1-12
 @centrecog cog
 @with   s:OW
 @cut   1.0
 @grid  10
 @traj
 tTOL00298K_500ps.trj.gz

 @endverbatim
 *
 * <hr>
 */

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"


#include <string>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "atomsfit", "pbc", "ref", "centre", "with", "centrecog", "cut", "grid"};
  int nknowns = 10;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc<boundary type>\n";
  usage += "\t@atomsfit <atomspecifier that determines the atoms to be used for the fitting>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@centre <atomspecifier>\n";
  usage += "\t[@centrecog] take cog for centre atoms\n";
  usage += "\t@with   <atomspecifier> or\n";
  usage += "\t@cut    <maximum distance>\n";
  usage += "\t@grid   <number of points>\n";
  usage += "\t@traj <trajectory files>\n";


  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());

    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
    //ref coordinates : only solute
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
     //if ref not given use first step of traj
      ic.open(args["traj"]);
    }
  
    ic >> refSys;
    ic.close();

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    
    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    
    Reference ref(&refSys);
 //////REF OK SOLUTE ONLY

    //get the atoms for the fit
    AtomSpecifier fitatoms(refSys);
    {
      Arguments::const_iterator iter = args.lower_bound("atomsfit");
      Arguments::const_iterator to = args.upper_bound("atomsfit");
      for(;iter!=to;iter++){
	string spec=iter->second.c_str();
	fitatoms.addSpecifier(spec);
      }
    }
    
    ref.addAtomSpecifier(fitatoms);
    //NOW REFERENCE HAVE THE SYST FROM TOPO + THE ATOMS SPEC FOR THE FIT

    
    // set centre atoms : for setting the with atom see below because need to create sys
  AtomSpecifier centre(refSys);
  {
    Arguments::const_iterator iter=args.lower_bound("centre");
    Arguments::const_iterator to=args.upper_bound("centre");
    for(;iter!=to;iter++)
      centre.addSpecifier(iter->second.c_str());
  }
    

  bool cog=false;
  if(args.count("centrecog")>=0) cog=true;
  

  // read in cut-off distance
  double cut=1.0;
  if(args.count("cut")>0) cut=atof(args["cut"].c_str());
  
  // read in grid number
  int grid=10;
  if(args.count("grid")>0) grid=atoi(args["grid"].c_str());
  
   //calculate here the cell side size and other stuffs that are needed for the Index
   double cell = cut/double(grid);
   grid=2*grid+1;
   int grid1=grid*grid;
   int gridsize=grid1*grid;
   int grid05=static_cast<int>(gridsize/2);
   int grid105=static_cast<int>(grid1/2);
   int grid005=static_cast<int>(grid/2);
   
   // System for calculation OK because use the same topo for the system the coordinates are not yet assigned!
    System sys(refSys);

  // set atom to consider with sys that now is able to contain the solvent
  AtomSpecifier with(sys);
  {
    Arguments::const_iterator iter=args.lower_bound("with");
    Arguments::const_iterator to=args.upper_bound("with");
    for(;iter!=to;iter++)
      with.addSpecifier(iter->second.c_str());
  }
    
    
    // Parse boundary conditions for sys and create object rf in order to use the RotationalFit methods
    pbc = BoundaryParser::boundary(sys, args);
    RotationalFit rf(&ref);
    
    
    int numFrames = 0;
  
    //vector to store the occurence of finding solvent in the grid and it the bissector vector 
    //also init the apos for averaging and spos for the temp positions
    vector<int> occur(gridsize,0);
    vector<Vec> bvec(gridsize,Vec(0.0,0.0,0.0));
    vector<Vec> apos;
    Vec vcog(0.0,0.0,0.0);
    Vec zero(0.0,0.0,0.0);
    Vec spos(0.0,0.0,0.0);
    for (int i=0;i<sys.mol(0).numPos();i++) {
	    apos.push_back(zero);
    }

    //loop over trajectory
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      ic.select("ALL");
      // loop over all frames
      while(!ic.eof()){
	numFrames++;
	ic >> sys;


      //gather the system applying the pbc and make the fit
	(*pbc.*gathmethod)();

      //the rotational fit does not work for molecule less than 4 atoms
	if (!(sys.mol(0).numPos()<=3)) {
		rf.fit(&sys);
      }

      //push back the solute position
      for (int i=0;i<sys.mol(0).numPos();i++) {
              spos = sys.mol(0).pos(i);
	      apos[i] += spos;
      }
      

      // just to make absolutely sure that we treat the centre and with
      // always in the same way, sort them
      centre.sort();
      with.sort();

      if(!cog){
              int gx, gy, gz, Index;
	      for(int i=0; i<centre.size();i++){  //loop over the centre atoms
	              for(int j=0; j<with.size();j++){   //loop over the atoms to consider
                              double tmp;
                              Vec tmpwith, tmpvdist;
			      tmpwith=pbc->nearestImage(*centre.coord(i),*with.coord(j),sys.box());
			      tmpvdist=tmpwith-*centre.coord(i);
			      tmp=(tmpvdist).abs();
                              if(tmp<=cut){
				      Vec tmpv2, tmpv3, bvtp1, bvtp2, bv;
				      gx=static_cast<int>(rint(tmpvdist[0]/cell));
				      gy=static_cast<int>(rint(tmpvdist[1]/cell));
				      gz=static_cast<int>(rint(tmpvdist[2]/cell));
                                      Index = gz * grid1 + gy * grid + gx + grid05;
			      
                              //NOW COMPUTE THE BISSECTOR
			      //ATTENTION: calculate the bissector of water ONLY
		              tmpv2=pbc->nearestImage(*centre.coord(i),
					      sys.sol(0).pos(with.atom(j)+1),
					      sys.box()); 		 
		              tmpv3=pbc->nearestImage(*centre.coord(i),
					      sys.sol(0).pos(with.atom(j)+2),
					      sys.box()); 
			      bvtp1 = tmpv3 - tmpwith;
			      bvtp2 = tmpv2 - tmpwith;
			      bv = (bvtp2 + bvtp1)/((bvtp2 + bvtp1).abs());

                              //NOW FILLED THE VECTORS
			      bvec[Index] += bv;
			      occur[Index]++;
			      
			      }
                      }
              }
      }

      else{
	// or if we want to do it for the cog
	(*pbc.*gathmethod)();

	for(int i=0; i<centre.size(); i++)
	  vcog+=*centre.coord(i);
	  vcog/=centre.size();

	// loop over the atoms to consider
              int gx, gy, gz, Index;
	      for(int j=0; j<with.size();j++){   //loop over the atoms to consider
                  double tmp;
                  Vec tmpwith, tmpvdist;
		  tmpwith=pbc->nearestImage(vcog,*with.coord(j),sys.box());
		  tmpvdist=tmpwith-vcog;
		  tmp=(tmpvdist).abs();
                  if(tmp<=cut){
		          Vec tmpv2, tmpv3, bvtp1, bvtp2, bv;
			  gx=static_cast<int>(rint(tmpvdist[0]/cell));
			  gy=static_cast<int>(rint(tmpvdist[1]/cell));
			  gz=static_cast<int>(rint(tmpvdist[2]/cell));
                          Index = gz * grid1 + gy * grid + gx + grid05;
		  
                  //NOW COMPUTE THE BISSECTOR
		  //ATTENTION: calculate the bissector of water ONLY
		  tmpv2=pbc->nearestImage(vcog,
		    	      sys.sol(0).pos(with.atom(j)+1),
		    	      sys.box()); 		 
		  tmpv3=pbc->nearestImage(vcog,
		    	      sys.sol(0).pos(with.atom(j)+2),
		    	      sys.box()); 
		  bvtp1 = tmpv3 - tmpwith;
		  bvtp2 = tmpv2 - tmpwith;
		  bv = (bvtp2 + bvtp1)/((bvtp2 + bvtp1).abs());
		  //bv = (bvtp2 + bvtp1)/2.0;

                  //NOW FILLED THE VECTORS
		  bvec[Index] += bv;
		  occur[Index]++;
		  
		  }
              }
      }
	
	//write out the trajectory
          {
              //OutG96 oc;
              //Outvmdam oc;
              //oc.open(cout);
              //oc.select("ALL");
              //oc.writeTitle("Test system");
              //oc << sys;
          }

      } // end loop over the number of frame in a traj file
      ic.close();
    } //end loop over trajectories files
    
    //average solute positions, write them out following the grid
    ofstream file_solute("solute_pos.dat");
    if(!cog){ 
	    for(int i=0; i<centre.size();i++){  //loop over the centre atoms
            for (int j = 0; j < (int) apos.size(); ++j) { 
                     apos[j] = apos[j]/numFrames;
	             Vec tmpv;
		     int gx, gy, gz;
		     tmpv=apos[j]-*centre.coord(i);
	             gx=static_cast<int>(rint(tmpv[0]/cell));
	             gy=static_cast<int>(rint(tmpv[1]/cell));
	             gz=static_cast<int>(rint(tmpv[2]/cell));
             const string a = sys.mol(0).topology().atom(j).name();
             char b[10]; 
             strcpy(b, a.c_str()); 
	     file_solute << setw(8)<< *b << setw(8) << gx << setw(8) 
	             << gy << setw(8) << gz << endl;
            }
            }
    }
    else {
	    for (int j = 0; j < (int) apos.size(); ++j) { 
                     apos[j] = apos[j]/numFrames;
	             Vec tmpv;
		     int gx, gy, gz;
		     tmpv=apos[j]-vcog;
	             gx=static_cast<int>(rint(tmpv[0]/cell));
	             gy=static_cast<int>(rint(tmpv[1]/cell));
	             gz=static_cast<int>(rint(tmpv[2]/cell));
             const string a = sys.mol(0).topology().atom(j).name();
             char b[10]; 
             strcpy(b, a.c_str()); 
	     file_solute << setw(8)<< *b << setw(8) << gx << setw(8) 
	             << gy << setw(8) << gz << endl;
            }
    }
    file_solute.close();

    //average solute positions, write them out
    //for (int i = 0; i < (int) apos.size(); ++i) { 
    //	    apos[i] = apos[i]/numFrames;
    //      cout.precision(3);
    //      cout << setw(8) << 10*apos[i][0] << setw(8) 
    //      << 10*apos[i][1] << setw(8) << 10*apos[i][2] << endl;
    //}
    
    //write out the density and the bissectors
    //now correct the distribution for the number of frames and the number 
    //of centre atoms
    //cout << "# number of frames considered: " << numFrames << endl;
    int divide=numFrames;

    if (!cog) divide*=centre.size();
    int ggx, ggy, ggz;
    double percent[gridsize];
    for(int i=0; i<=gridsize-1;i++){  
    if (occur[i]!=0) {
       bvec[i] /= double(occur[i]);
       percent[i] = double(occur[i])/double(divide);
       ggz = (i-grid105-grid05)/(grid1); 
       if(ggz>=0) ggz = (i+grid105-grid05)/(grid1); 
       ggy = (i-grid005-grid05-ggz*grid1)/(grid); 
       if(ggy>=0) ggy = (i+grid005-grid05-ggz*grid1)/(grid);
       ggx = i-grid05-ggz*grid1-ggy*grid; 
       //IIndex = ggz * grid1 + ggy * grid + ggx + grid05;
       cout.precision(8);
       cout <<setw(10)<< i << "\t"<<setw(8) << percent[i]<< setw(14) 
	    << bvec[i][0]<< setw(14) << bvec[i][1]<< setw(14) << bvec[i][2] 
	    << setw(10) << ggx << setw(5) << ggy<< setw(5) << ggz<< endl;
       }
    }
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


