// oparam.cc

#include "../src/args/Arguments.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"


#include <vector>
#include <iomanip>
#include <sstream>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "atoms", "type",  
		    "pbc",  "moln", "refvec"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@moln <number of molecules to include>\n";
  usage += "\t@refvec <x y z>\n";
  usage += "\t@atoms <atom(s) to calculate order parameter from>\n";
  usage += "\t@type <0 1>
             Choices are:
                   0: middle atom specification
                   1: all atom specification\n";
  usage += "\t@traj <trajectory files>\n";
  


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    //get reference vector, normalize it
    Vec refvec (0.0, 0.0, 0.0);

    {
    Arguments::const_iterator iter=args.lower_bound("refvec");
   if(iter!=args.upper_bound("refvec")){
     refvec[0]=atof(iter->second.c_str());
     ++iter; 
   }
   if(iter!=args.upper_bound("refvec")){
     refvec[1]=atof(iter->second.c_str());
     ++iter; 
   }
   if(iter!=args.upper_bound("refvec")){
     refvec[2]=atof(iter->second.c_str());
   }
    }
    //normalize 
    refvec = refvec.normalize();
   
    // read topology
    InTopology it(args["topo"]);
    //make system out of topology
    System sys(it.system());
    // if you want to work with a topology that only contains one molecule:
      // do this and it will add the appropriate number to the system -> actually not true
        //  for (i=0;i<moln-1;++i){
        //  sys.addMolecule(it);
        //	}
     // set molecule number
   int  moln=sys.numMolecules();
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     moln=atoi(iter->second.c_str());
   }
    }
   if (moln <= 0) {throw gromos::Exception("oparam", "molecule number cannot be <= 0!\n");}

   //get specification
   int  type=0;
   {
   Arguments::const_iterator iter=args.lower_bound("type");
   if(iter!=args.upper_bound("type")){
     type=atoi(iter->second.c_str());
   }
    }
   if (type < 0 || type > 1) {throw gromos::Exception("oparam", "type cannot be other than 0 or 1!\n");}

    //add the atoms to calculate the oparams

      vector<int> atoms; 
     for(Arguments::const_iterator iter=args.lower_bound("atoms");
	iter!=args.upper_bound("atoms"); ++iter){
       int arse = atoi(iter->second.c_str())-1;
       atoms.push_back(arse);
     }
   
    if (int (atoms.size()) == 0) {throw gromos::Exception("oparam", "at least one atom needs to be defined!\n");}
 
    // put stuff in other vector
    vector<int> at;    
    if (type == 0) {
      for (int i=0; i < int (atoms.size()); ++i){
        at.push_back(atoms[i]-1);
      if (atoms[i-1] == 0) {throw gromos::Exception("oparam", "cannot calculate the oparam for the 1st atom!\n");}
        at.push_back(atoms[i]);
        at.push_back(atoms[i]+1);
      }
    }
    else if (type == 1) {
     if ((int (atoms.size()))%3 != 0) {      
     throw gromos::Exception("oparam", "defined atoms are not a multiple of 3!\n");
    }   
      at = atoms;
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

    
    // loop over all trajectories
    int numFrames = 0; 
    Vec z (0.0,0.0,0.0);
    Vec y (0.0,0.0,0.0); 
    Vec x (0.0,0.0,0.0); 
    Vec cos (0.0,0.0,0.0);
    Vec scos2 (0.0,0.0,0.0);
    vector<Vec> S; for (int i=0;i< int (at.size()/3);++i) {S.push_back(z);} 
    vector<Vec> avcos2; for (int i=0;i< int (at.size()/3);++i) {avcos2.push_back(z);}   
   
    // define input coordinate
    InG96 ic;
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      
      // loop over all frames
      while(!ic.eof()){
      numFrames++;
	ic >> sys;
	pbc->gathergr();
        
	// calculate the z-vector between atom i-1 and i+1, normalize
        int cc=-2;     
	for (int i=0;i<moln;++i){
          cc=-2;
	  for (int j=0;j< int (at.size()/3); j++){
	    cc+=2;
	  z=sys.mol(i).pos(at[j+2+cc])-sys.mol(i).pos(at[j+cc]);
          z = z.normalize();
	//calculate y, normalize
          y=(sys.mol(i).pos(at[j+1+cc])-sys.mol(i).pos(at[j+cc]));
          y=y-((sys.mol(i).pos(at[j+1+cc])-(sys.mol(i).pos(at[j+cc]))).dot(z))*z;
	  y = y.normalize();
	  //calculate x
          x=z.cross(y);
	    x = x.normalize(); //is this nescessary?
          // determine the angle
          cos[0]=refvec.dot(x);
          cos[1]=refvec.dot(y);
          cos[2]=refvec.dot(z);
          scos2=cos*cos;
          avcos2[j]+=scos2;
	 }
	}

      }
    }
      ic.close();
    
    // average the avcos2 finally
    for (int i=0; i < int (at.size()/3); ++i){
     avcos2[i]/=(numFrames*moln);
    }

    //get the S components
    for (int i=0; i< int (at.size()/3); ++i){
       Vec tmp=avcos2[i];
       tmp[0] = ((3*tmp[0])-1)/2;
       tmp[1] = ((3*tmp[1])-1)/2;
       tmp[2] = ((3*tmp[2])-1)/2;
      S[i] = tmp; 
    }

    //determine the  SCDOP and SCDAP, SCDEP
    vector<double> SCDOP, SCDAP, SCDEB;
    for (int i=0; i< int (at.size()/3); ++i){
      Vec tmp = avcos2[i];
      Vec tmp1 = S[i];
      SCDOP.push_back(-(tmp[0]+(tmp[1]-1)/2));
      SCDAP.push_back(tmp1[2]/2);
      SCDEB.push_back((2*tmp1[0]/3)+tmp1[1]/3);
    }
    
      //print out results
    cout << setw(4) << "Atom" 
         << setw(8) << "SX" << setw(8) << "SY" << setw(8) << "SZ" 
         << setw(8) << "SCDOP" 
         << setw(8) << "SCDAP" 
         << setw(8) << "SCDEB" << endl;
      cout << endl;

      int c=-1; 
      for (int i=0; i < int (at.size()/3); ++i) {
        c+=2;
        Vec tmp = S[i];
	cout.setf(ios::floatfield, ios_base::fixed);
        cout.setf(ios::right, ios::adjustfield);
        cout.precision(4); 
        cout << setw(4) << at[i+c]+1 
             << setw(8) <<   tmp[0] << setw(8) << tmp[1] << setw(8) << tmp[2] 
             << setw(8) << SCDOP[i] << setw(8) << SCDAP[i] << setw(8) << SCDEB[i] << endl;
      }      
  } 
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

