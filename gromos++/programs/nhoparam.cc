// nhoparam.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/Vacuum.h"
#include "../src/bound/RectBox.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Exclusion.h"
#include "../src/utils/Neighbours.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gmath/Vec.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cstdlib>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace utils;
using namespace fit;

vector <double> Matbuil(Vec S, vector<double> sum);
double S2calc(vector<double> S, int frnum);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "atoms", "type", "ref", "class", "oatoms", 
		    "pbc",  "moln", "mol", "time", "winframe"};
  int nknowns = 12;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@moln <molecule number for S2 calculation>\n";
  usage += "\t@mol <molecules considered in fit>\n";
  usage += "\t@oatoms <atom(s) to calculate order parameter from>\n";
  usage += "\t@type <atom type of the calculation<for now only NH>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@atoms <atoms to consider for fit>\n";
  usage += "\t@time <time and dt>\n";
  usage += "\t@winframe <averaging window [# of frames]>\n";
  usage += "\t@traj <trajectory files>\n";
  


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    //make system out of topology
    System refSys(it.system());

   int moln = 0; 
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     moln=atoi(iter->second.c_str())-1;
   }
    }

   int winframe=1;
   {
   Arguments::const_iterator iter=args.lower_bound("winframe");
   if(iter!=args.upper_bound("winframe")){
     winframe=atoi(iter->second.c_str());
   }
    }



   //pop in ref coordinates
       InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();
  
    // get simulation time
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

   

    // parse type
    string t = "NH";
    try{
     string format = args["type"];
     if(format == "NH")
      t = "NH";
     else 
      throw gromos::Exception("nhoparam","include format " +format+ " unknown.\n");
    }
      catch(Arguments::Exception &e){}




    //add the atoms to calculate the oparams

      vector<int> atoms; 
     for(Arguments::const_iterator iter=args.lower_bound("oatoms");
	iter!=args.upper_bound("oatoms"); ++iter){
       int arse = atoi(iter->second.c_str())-1;
       if (arse == 0) {throw gromos::Exception("oparam", "cannot calculate the oparam for the 1st atom!\n");}
       atoms.push_back(arse);
     }
   
     //  if (int (atoms.size()) == 0) {throw gromos::Exception("oparam", "at least one atom needs to be defined!\n");}
    

    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    // gather reference system
    pbc->gathergr();
    delete pbc;
    
    Reference ref(&refSys);

    // Adding references
    int added=0;
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<refSys.numMolecules();++i)
        mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
          it!=args.upper_bound("mol");++it){
        if(atoi(it->second.c_str())>refSys.numMolecules())
          throw Arguments::Exception(usage);
        mols.push_back(atoi(it->second.c_str())-1);
      }
    // add classes
    for(Arguments::const_iterator it=args.lower_bound("class");
        it != args.upper_bound("class"); ++it){
      for(vector<int>::const_iterator mol=mols.begin();
          mol!=mols.end();++mol)
        ref.addClass(*mol,it->second);
      added=1;
    }
    // add single atoms
    for(Arguments::const_iterator it=args.lower_bound("atoms");
        it != args.upper_bound("atoms"); ++it){
      int atom=atoi(it->second.c_str())-1, mol=0;
      while(atom >= refSys.mol(mol).numAtoms()){
        atom-=refSys.mol(mol).numAtoms();
        ++mol;
        if(mol==refSys.numMolecules())
          throw Arguments::Exception(usage);
      }
      ref.addAtom(mol,atom);
      added=1;
    }
    // did we add anything at all?
    if(!added)
      throw Arguments::Exception(usage);

    // System for calculation
    System sys(refSys);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&ref);

    //see whether no atoms are supplied
    if (atoms[0] == -1) {
     atoms.resize(0);
      for (int i=0; i < sys.mol(moln).numAtoms(); ++i){
       if (sys.mol(moln).topology().atom(i).mass()==14.00670){
       Neighbours neigh(sys,moln,i);
       for(int j=0; j< int (neigh.size()); ++j){ 
        if (sys.mol(moln).topology().atom(neigh[j]).mass()==1.00800){
	  atoms.push_back(i);
          break;
	}
       }
       }
      }
    } 

    //set up atom check, smack in the bounded H atoms
    vector<Exclusion> Hbound;

    for (int i=0; i< int (atoms.size()); ++i){
      if (sys.mol(moln).topology().atom(atoms[i]).mass()!=14.00670){
	throw gromos::Exception("nhoparam", "atom not a Nitrogen! Check your topology!\n");}
      Neighbours neigh(sys,moln,atoms[i]);
      vector<int> neiat; Exclusion tmp, tmp1;
      for(int j=0; j< int (neigh.size()); ++j){ 
       if (sys.mol(moln).topology().atom(neigh[j]).mass()==1.00800){
	 tmp.insert(neigh[j]);
       }
      }
      Hbound.push_back(tmp);
      tmp=tmp1;
    }


    
    // loop over all trajectories
    int numFrames = 0, frcount=0; 
    Vec nh (0.0,0.0,0.0);
    Vec nh1, nh2; nh1=nh; nh2=nh;
    vector<double> S2; for (int i=0;i< int (atoms.size());++i) {S2.push_back(0);}
    vector<double> S2_win; for (int i=0;i< int (atoms.size());++i) {S2_win.push_back(0);}
    vector<double> S2_w; for (int i=0;i< int (atoms.size());++i) {S2_w.push_back(0);}
    vector<double> sum, tmp, erase; for (int i=0;i< 9;++i) {sum.push_back(0); tmp.push_back(0); erase.push_back(0);}
    vector<vector<double> > S2_sum;  for (int i=0;i< int (atoms.size());++i) {S2_sum.push_back(sum);}
    vector<vector<double> > S2_sumwin;  for (int i=0;i< int (atoms.size());++i) {S2_sumwin.push_back(sum);} 
    double d_ts=0, d_tsw=0;
    ofstream ts; ts.open("OPts.out");
    ofstream tsw; tsw.open("OPwints.out");
    ts.setf(ios::floatfield, ios_base::fixed);
        ts.setf(ios::right, ios::adjustfield);
        ts.precision(4);
    tsw.setf(ios::floatfield, ios_base::fixed);
        tsw.setf(ios::right, ios::adjustfield);
        tsw.precision(4);
    for (int i=0; i < int (atoms.size()); ++i) {
      ts << setw(5) << atoms[i]+1;
      tsw << setw(5) << atoms[i]+1;
      }
    ts << endl; tsw << endl;
    
    
    // define input coordinate
    for(Arguments::const_iterator iter=args.lower_bound("traj");
	iter!=args.upper_bound("traj"); ++iter){
      ic.open(iter->second);
      // loop over all frames
      while(!ic.eof()){
	numFrames++; frcount++;
	ic >> sys;
	pbc->gathergr();
        rf.fit(&sys);
        ts << setw(6) << time;
	// calculate the z-vector between atom i-1 and i+1, normalize
 	 for (int i=0;i< int (atoms.size()); i++){
	   Exclusion tmpex = Hbound[i];
          
	   switch(tmpex.size()) {
	    case 1: 
             nh=sys.mol(moln).pos(tmpex.atom(0))-sys.mol(moln).pos(atoms[i]);
	     nh=nh.normalize();             
             S2_sum[i] = Matbuil(nh, S2_sum[i]);
             S2_sumwin[i] = Matbuil(nh, S2_sumwin[i]);
             tmp = S2_sum[i]; 
             d_ts= S2calc(tmp, numFrames); 
             ts << setw(8) << d_ts;
             break;
            case 2: 
	     nh=sys.mol(moln).pos(tmpex.atom(0))-sys.mol(moln).pos(atoms[i]);
             nh1=sys.mol(moln).pos(tmpex.atom(1))-sys.mol(moln).pos(atoms[i]);
             nh=nh.normalize();  nh1=nh1.normalize();
             nh=(nh+nh1)/2;
             nh=nh.normalize();
             S2_sum[i] = Matbuil(nh, S2_sum[i]);
             S2_sumwin[i] = Matbuil(nh, S2_sumwin[i]);
             tmp = S2_sum[i]; 
             d_ts= S2calc(tmp, numFrames);          
             ts << setw(8) << d_ts;
             break;
            case 3: 
	     nh=sys.mol(moln).pos(tmpex.atom(0))-sys.mol(moln).pos(atoms[i]);
             nh1=sys.mol(moln).pos(tmpex.atom(1))-sys.mol(moln).pos(atoms[i]);
             nh2=sys.mol(moln).pos(tmpex.atom(2))-sys.mol(moln).pos(atoms[i]);
             nh=nh.normalize();  nh1=nh1.normalize();  nh2=nh2.normalize();
             nh=(nh+nh1+nh2)/3; nh=nh.normalize();
             S2_sum[i] = Matbuil(nh, S2_sum[i]);
             S2_sumwin[i] = Matbuil(nh, S2_sumwin[i]);
             tmp = S2_sum[i]; 
             d_ts= S2calc(tmp, numFrames);   
             ts << setw(8) << d_ts;
             break;
	    case 4:
             throw gromos::Exception("nhoparam", "i think you have an error in your topology, e.g. 4 hydrogens bound to a Nitrogen!\n");
	   }
	 }
         ts << endl;
         time+=dt;
         if (frcount == winframe){
	   tsw << setw(6) << time;
         for (int i=0; i < int (atoms.size()); ++i){ 
          tmp = S2_sumwin[i]; 
          d_tsw= S2calc(tmp, winframe);
          S2_w[i]+=d_tsw;
          tsw << setw(8) << d_tsw; 
          S2_sumwin[i]=erase;          
	 }
	 frcount = 0; tsw << endl;
         } 
      }	    
      ic.close();
    }

//calculate S2, S2_win
    int ntimes=numFrames/winframe;
 for (int i=0; i < int (atoms.size()); ++i){
  tmp=S2_sum[i];
  S2[i]= S2calc(tmp, numFrames);
  S2_w[i]/=ntimes;
}

    
      //print out results
 cout << "S2 over all frames:" << endl;
   cout << setw(4) << "Atom" << setw(8) << "Type" << setw(11) << "#  AS"
         << setw(8) << "S2" << setw(8)
         << endl;
    cout << endl;

      for (int i=0; i < int (atoms.size()); ++i) {
        cout.setf(ios::floatfield, ios_base::fixed);
        cout.setf(ios::right, ios::adjustfield);
        cout.precision(4);
        cout << setw(4) << atoms[i]+1
             << setw(8) << sys.mol(moln).topology().atom(atoms[i]).name()
             << setw(8) << sys.mol(moln).topology().resNum(atoms[i])+1
             << setw(4) << sys.mol(moln).topology().resName(sys.mol(moln).topology().resNum(atoms[i]))
             << setw(10) <<   S2[i]
             << endl;
      }
  cout << endl; 
 ts.close();
 tsw.close();

 cout << "S2 over all windows:" << endl;
   cout << setw(4) << "Atom" << setw(8) << "Type" << setw(11) << "#  AS"
         << setw(8) << "S2" << setw(8)
         << endl;
    cout << endl;

      for (int i=0; i < int (atoms.size()); ++i) {
        cout.setf(ios::floatfield, ios_base::fixed);
        cout.setf(ios::right, ios::adjustfield);
        cout.precision(4);
        cout << setw(4) << atoms[i]+1
             << setw(8) << sys.mol(moln).topology().atom(atoms[i]).name()
             << setw(8) << sys.mol(moln).topology().resNum(atoms[i])+1
             << setw(4) << sys.mol(moln).topology().resName(sys.mol(moln).topology().resNum(atoms[i]))
             << setw(10) << S2_w[i]
             << endl;
      }



  }

	     
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


vector<double> Matbuil(Vec S, vector<double> sum){
 
  
 for (int i=0; i<3;i++){
  for (int j=0; j<3;j++){ 
   sum[3*i+j]+=(S[i]*S[j]);
}
 }

return (sum);
}

double S2calc(vector<double> S, int frnum){
 double Ssq=0;
  
 for (int i=0; i<9;i++){
   S[i]=S[i]/frnum;
   Ssq+= S[i]*S[i];
 }

return (1.5*Ssq-0.5);
}
