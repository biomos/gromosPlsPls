//Intramolecular Hydrogen bond analysis Hbond 
//dedicated to wilfred: "you know, mika: proahb IS my favorite program"
//--mika

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/Ginstream.h"
#include "../src/utils/Neighbours.h"
#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>

using namespace utils;
using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;

ostream &operator<<(ostream &os, const vector<int>& v){
  for (vector<int>::const_iterator i = v.begin(); i !=v.end();++i)
    os << *i << " ";
    return os;
  }
ostream &operator<<(ostream &os, const vector<double>& v){
  for (vector<double>::const_iterator i = v.begin(); i !=v.end();++i)
    os << *i << " ";
    return os;
  }

int main(int argc, char **argv){

  char *knowns[] = {"topo", "pbc", "moln", "Hbparas", "time", "traj"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc    <boundary type>\n";
  usage += "\t@moln   <molecule number in topology 0..X>\n";
  usage += "\t@Hbparas <distance [nm] and angle; default: 0.2, 135>\n";
  usage += "\t@time   <time and dt>\n";
  usage += "\t@traj   <trajectory files>\n";
  
 
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
  

// set molecule number
   int  molecule=0;
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     molecule=atoi(iter->second.c_str())-1;
   }
    }
  // get the paras
  double maxdist=0.2, minangle = 135;
  {
    Arguments::const_iterator iter=args.lower_bound("Hbparas");
    if(iter!=args.upper_bound("Hbparas")){
      maxdist=atof(iter->second.c_str());
      ++iter;
    }
    if(iter!=args.upper_bound("Hbparas"))
        minangle=atof(iter->second.c_str());
  }

  //  read topology
   args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
    // parse boundary conditions
  Boundary *pbc = BoundaryParser::boundary(sys, args);

    // define input coordinate
  InG96 ic;

    //find all H, O, N and S atoms 
  vector<int> Hatoms, Acceptors, Oatoms, Natoms, Satoms, hto;
  for (int j=0;j<sys.mol(molecule).numAtoms();++j){
       if (sys.mol(molecule).topology().atom(j).mass()==1.00800){
	 Hatoms.push_back(j); Neighbours neigh(sys,molecule,j); hto.push_back(neigh[0]);
	  }
      else if (sys.mol(molecule).topology().atom(j).mass()==15.99940){
	Acceptors.push_back(j);
        Oatoms.push_back(j);
	  }
      else if (sys.mol(molecule).topology().atom(j).mass()==14.00670){
	Acceptors.push_back(j);
        Natoms.push_back(j);
	  }
      else if (sys.mol(molecule).topology().atom(j).mass()==32.06000){
	Acceptors.push_back(j);
        Satoms.push_back(j);
	  } 
  }
 
  //build donor, acceptor list, spit it out
   if (Hatoms.size() == 0) {throw gromos::Exception("Hbond", "No Donors found!\n");}
   if (Acceptors.size() == 0) {throw gromos::Exception("Hbond", "No Acceptors found!\n");}


     int num =1;
     for(int i=0;i< int (Hatoms.size());++i){
      for(int j=0;j< int (Acceptors.size());++j){   
	if (hto[i] != Acceptors[j]){
	  num +=1;
	    }
      }
     }             

      cout << endl;
      cout << "Potential Hydrogen Bonds to be monitored:" << num <<endl;
      cout << endl;
      cout << "Starting the run..." << endl;
      vector<double> angtot, disttot, occur;
      for (int i=0;i<num-1;++i){angtot.push_back(0.0);disttot.push_back(0.0);occur.push_back(0);}
      int numFrames = 0, numHbpframe=0;
      vector<double> ti, hn;
      ofstream tsh; tsh.open("Hbnumts.out");
    // loop over all trajectories
      for(Arguments::const_iterator 
      iter=args.lower_bound("traj"),
      to=args.upper_bound("traj");
      iter!=to; ++iter){

     // open file
      ic.open((iter->second).c_str());

   // loop over single trajectory
   
      while(!ic.eof()){
	numFrames++;
      ic >> sys;
      pbc->gather();
      

      double dist = 0, angle=0;
      numHbpframe=0;
      num =0; 
    for(int i=0;i< int (Hatoms.size());++i){
     for(int j=0;j<int (Acceptors.size()) ;++j){   
     if (hto[i] != Acceptors[j]){
       num += 1;
      Vec tmp = (sys.mol(molecule).pos(Hatoms[i])-sys.mol(molecule).pos(Acceptors[j])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){
       Vec tmpA = (sys.mol(molecule).pos(Acceptors[j])-sys.mol(molecule).pos(Hatoms[i])); 
       Vec tmpB = (sys.mol(molecule).pos(hto[i])-sys.mol(molecule).pos(Hatoms[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
        disttot[num-1]+= dist; 
        angtot[num-1]+= angle;
        occur[num-1]+= 1;
        numHbpframe+=1;
        ti.push_back(time); hn.push_back(num);
       }
      }
     }
     }
    }    
        tsh.precision(6);
        tsh << setw(10) << time;
        tsh.precision(5);
        tsh << setw(10) << numHbpframe << endl;
        time += dt;
      }

      }
        ic.close();
        tsh.close();
	//place the averaging in here
      cout << endl;
    cout << "Statistics of the run:" << endl;
    cout << "    #     #HB"   <<   "  Donor" << "  " << "Acceptor" << "       " 
         << "D-        " << "H...        " << "A  " 
         <<" DIST "<< "  ANGLE "<< "   OCCUR "
         << "   % " << endl;

    ofstream per90("90.per");
    ofstream per70("70.per");
    ofstream per50("50.per");
    ofstream per30("30.per");
    ofstream per10("10.per");    
    ofstream per0("0.per");
    num=1;
    int k=0;
    vector<double> cnum,ck;
     for(int i=0;i< int (Hatoms.size());++i){
      for(int j=0;j< int (Acceptors.size());++j){   
	if (hto[i] != Acceptors[j]){
	  k +=1;
         if (occur[k-1]!=0){
        cout.setf(ios::floatfield, ios_base::fixed);
        cout.setf(ios::right, ios::adjustfield);
        cout.precision(3); 
	cout << setw(5) << num << setw(7) << k <<  
         setw(4) << sys.mol(molecule).topology().resNum(Hatoms[i])+1
             << sys.mol(molecule).topology().resName(sys.mol(molecule).topology().resNum(Hatoms[i])) << setw(2) << "-" <<
         setw(4) << sys.mol(molecule).topology().resNum(Acceptors[j])+1
             << sys.mol(molecule).topology().resName(sys.mol(molecule).topology().resNum(Acceptors[j])) 
             << setw(6) << hto[i]+1 << setw(4) << sys.mol(molecule).topology().atom(hto[i]).name() << "-"
             << setw(6) << Hatoms[i]+1  << setw(4) << sys.mol(molecule).topology().atom(Hatoms[i]).name() << "-"
             << setw(6) << Acceptors[j]+1 << setw(4) << sys.mol(molecule).topology().atom(Acceptors[j]).name() 
             << setw (7) << (disttot[k-1]/occur[k-1]);
        cout.precision(2);
        cout << setw(8)  << (angtot[k-1]/occur[k-1]);
        cout.precision(0); 
        cout << setw(7)  << (occur[k-1]);
        cout.setf(ios::floatfield, ios_base::fixed);
        cout.precision(2);             
        double per = (occur[k-1]/numFrames)*100;
        cout << setw(8)  << per 
             << endl;
        
        if (per > 90){
	  per90 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
        else if (per > 70){
          per70 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
        else if (per > 50){
          per50 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
        else if (per > 30){
          per30 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
        else if (per > 10){
          per10 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
        else if (per > 0){
          per0 << sys.mol(molecule).topology().resNum(Hatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(Acceptors[j])+1 << endl;}
         
          cnum.push_back(num);
          ck.push_back(k);
          num +=1;
	 }
	}
      }
     }
   


     for (int i=0;i< int (ti.size());++i){
       for (int j=0;j<int (ck.size());++j){
	 if (hn[i] == ck[j]){
	   hn[i]=cnum[j];
           break;
	 }
       }
     }
 
     ofstream ts; ts.open("Hbts.out");
     for (int i=0;i<int (ti.size());++i){
       ts << ti[i] << " " << hn[i] << endl;
     }
     ts.close(); 
     per90.close();per70.close();per50.close();per30.close();per10.close();per0.close();    
    }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
