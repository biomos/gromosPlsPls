// frameout.cc

#include "../src/args/Arguments.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/fit/TranslationalFit.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96S.h"
#include "../src/gio/OutPdb.h"
#include "../src/gio/Outvmdam.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/Solvent.h"
#include "../src/gcore/SolventTopology.h"
#include "../src/gmath/Vec.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cstdlib>
#include <iostream>

using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj",  "gather",  
		    "pbc", "spec", "frames", "outformat", "include", "center"};
  int nknowns = 9;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@gather <gather type, can be GNORM (default), GGR or COGG>\n";
  usage += "\t@center <Vector to center mol(0) in case of COGG>\n";
  usage += "\t@spec <specification for writing out frames. either ALL, EVERY or SPEC>\n";
  usage += "\t@frames <frames to be written out>\n";
  usage += "\t@outformat <output format. either pdb, g96 or vmdam>\n"; 
  usage += "\t@include <either SOLUTE (default), SOLVENT or ALL>\n";
  usage += "\t@traj <trajectory files>\n";
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);
   
    //smack in the framenumbers
    vector<int> fnum;
    for(Arguments::const_iterator it=args.lower_bound("frames");
	  it != args.upper_bound("frames"); ++it){
	int bla=atoi(it->second.c_str());
        fnum.push_back(bla);
    }      

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    // get centering vector    
    Vec center (0.0, 0.0, 0.0);

    {
    Arguments::const_iterator iter=args.lower_bound("center");
   if(iter!=args.upper_bound("center")){
     center[0]=atof(iter->second.c_str());
     ++iter; 
   }
   if(iter!=args.upper_bound("center")){
     center[1]=atof(iter->second.c_str());
     ++iter; 
   }
   if(iter!=args.upper_bound("center")){
     center[2]=atof(iter->second.c_str());
   }
    }

    // Parse boundary conditions
    Boundary *pbc;
    try{
      char b=args["pbc"].c_str()[0];
      switch(b){
      case 't':
	pbc=new TruncOct(&sys);
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

    // parse includes
    string inc = "SOLUTE";
           try{
    string format = args["include"];
    if(format == "ALL")
      inc = "ALL";
    else if(format == "SOLUTE")
      inc = "SOLUTE";
    else if (format == "SOLVENT")
      inc = "SOLVENT";
    else 
      throw gromos::Exception("frameout","include format "+format+" unknown.\n");
    }
	    catch(Arguments::Exception &e){
      //  pbc = new Vacuum(&sys);
     }  
    // parse spec

    string spec = "ALL";
            try{
     string format = args["spec"];
    if(format == "ALL")
      spec = "ALL";
    else if(format == "EVERY")
      spec = "EVERY";
    else if (format == "SPEC")
      spec = "SPEC";
    else 
      throw gromos::Exception("frameout","spec format "+format+" unknown.\n");
    }
      catch(Arguments::Exception &e){
      //  pbc = new Vacuum(&sys);
    }  

    // outformat
     OutCoordinates *oc;
     string ext = ".g96";
    try{
      string format = args["outformat"];
      if(format == "pdb"){
       oc = new OutPdb();
       ext = ".pdb";}
      else if(format == "g96"){
        oc = new OutG96S();
        ext = ".g96";}
      else if (format == "vmdam"){
        oc = new Outvmdam();
        ext = ".vmd";}
      else 
        throw gromos::Exception("frameout","output format "+format+" unknown.\n");
    }
    catch(Arguments::Exception &){
      oc = new OutG96S();
    }

    // gather type
    string ga = "GNORM";
    try{
      string format = args["gather"];
      if(format == "GNORM")
        ga = "GNORM";
      else if(format == "COGG")
        ga = "COGG";
      else if(format == "GGR")
        ga = "GGR";
      else 
        throw gromos::Exception("frameout","gather format "+format+" unknown.\n");
    }
    catch(Arguments::Exception &){
      oc = new OutG96S();
    }
 
           
    // loop over all trajectories
    InG96 ic;
    int numFrames = 0;
    int every = 0;
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);  
      // loop over all frames
      while(!ic.eof()){
	numFrames++;every++;
	ic.select(inc.c_str());
	ic >> sys;
        if (ga == "GNORM"){
		pbc->gather();
	}
        else if(ga == "GGR") {pbc->gathergr();}
        else {pbc->coggather(center);}

	   if (spec == "ALL"){	
                                char outFile[]="FRAME";
	       ostrstream out;
	       if (numFrames < 10){
		 out <<outFile<<"_"<<"0000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 100){
                 out <<outFile<<"_"<<"000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 1000){
                 out <<outFile<<"_"<<"00"<<numFrames<<ext<<ends;
	       }
               else {
                out <<outFile<<"_"<<"0"<<numFrames<<ext<<ends;
	       }
               ofstream os(out.str()); 
	       oc->open (os);       
               oc-> select(inc.c_str());
	       oc->writeTitle(out.str());
	       
	       *oc << sys;
	       os.close();
	   }
		else if (spec == "EVERY" && (every-fnum[0] == 0)){
                 char outFile[]="FRAME";
	       ostrstream out;
	       if (numFrames < 10){
		 out <<outFile<<"_"<<"0000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 100){
                 out <<outFile<<"_"<<"000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 1000){
                 out <<outFile<<"_"<<"00"<<numFrames<<ext<<ends;
	       }
               else {
                out <<outFile<<"_"<<"0"<<numFrames<<ext<<ends;
	       } 
               ofstream os(out.str()); 
	       oc->open (os);       
               oc-> select(inc.c_str());
	       oc->writeTitle(out.str());
	       
	       *oc << sys;
	       os.close();	       
               every=0;
		}

		else if (spec == "SPEC"){
                 for (int i=0; i< int (fnum.size()); ++i){
		  if (numFrames == fnum[i]){
                  char outFile[]="FRAME";
	       ostrstream out;
	       if (numFrames < 10){
		 out <<outFile<<"_"<<"0000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 100){
                 out <<outFile<<"_"<<"000"<<numFrames<<ext<<ends;
	       }
               else if (numFrames < 1000){
                 out <<outFile<<"_"<<"00"<<numFrames<<ext<<ends;
	       }
               else {
                out <<outFile<<"_"<<"0"<<numFrames<<ext<<ends;
	       } 
	        ofstream os(out.str()); 
	       oc->open (os);       
               oc-> select(inc.c_str());
	       oc->writeTitle(out.str());
	       
	       *oc << sys;
	       os.close();
		  }
		 }
		}
    	
      }
   ic.close();
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
