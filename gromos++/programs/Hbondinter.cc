//Hydrogen bond analysis Hbondinter 
//dedicated to wilfred: "you know, mika: proahb IS my favorite program"
//--mika

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/fit/PositionUtils.h"
#include "../src/gmath/Vec.h"
#include "../src/gcore/Angle.h"
#include "../src/gcore/Bond.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/utils/Neighbours.h"

#include <vector>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <string>
#include <fstream>


using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;

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
  usage += "\t@moln   <molecule numbers in topology 1..X>\n";
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
  

// set molecule numbers
   int  m1=0, m2=1;
   {
   Arguments::const_iterator iter=args.lower_bound("moln");
   if(iter!=args.upper_bound("moln")){
     m1=atoi(iter->second.c_str())-1;
     ++iter;
   }
   if(iter!=args.upper_bound("moln"))
        m2=atoi(iter->second.c_str())-1;
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
  vector<int> Hm1, Hm2, Am1, Am2, Htom1, Htom2;  
   for (int j=0;j<sys.mol(m1).numAtoms();++j){  
       if (sys.mol(m1).topology().atom(j).mass()==1.00800){
	 Hm1.push_back(j); Neighbours neigh(sys,m1,j); Htom1.push_back(neigh[0]); 
       }
      else if (sys.mol(m1).topology().atom(j).mass()==15.99940){
	Am1.push_back(j);
      }  
        else if (sys.mol(m1).topology().atom(j).mass()==14.00670){
	Am1.push_back(j); 
	  }
      else if (sys.mol(m1).topology().atom(j).mass()==32.06000){
	Am1.push_back(j); 
	  } 
   }
   for (int j=0;j<sys.mol(m2).numAtoms();++j){  
       if (sys.mol(m2).topology().atom(j).mass()==1.00800){
	 Hm2.push_back(j);Neighbours neigh(sys,m2,j); Htom2.push_back(neigh[0]);   
       }
      else if (sys.mol(m2).topology().atom(j).mass()==15.99940){
	Am2.push_back(j);
      }  
        else if (sys.mol(m2).topology().atom(j).mass()==14.00670){
	Am2.push_back(j); 
	  }
      else if (sys.mol(m2).topology().atom(j).mass()==32.06000){
	Am2.push_back(j); 
	  } 
   }
  
 
  //build donor, acceptor list, spit it out
      cout << "Starting the run..." << endl;
      vector<double> angtotm1, disttotm1, occurm1,
                     angtotm2, disttotm2, occurm2;
      for (int i=0;i< int (Hm1.size()* Am2.size());++i)
       {angtotm1.push_back(0.0);disttotm1.push_back(0.0);occurm1.push_back(0);}
      for (int i=0;i< int (Hm2.size()* Am1.size());++i)
       {angtotm2.push_back(0.0);disttotm2.push_back(0.0);occurm2.push_back(0);}
      cout << "Potential Hydrogen Bonds to be monitored: " 
           << (Hm1.size()* Am2.size())+(Hm2.size()* Am1.size()) << endl;



      ofstream tsh; tsh.open("Hbnumts.out");
      int numFrames = 0, numHbpframe=0;
      Vec blaa(0.0,0.0,0.0);
      vector<double> ti; vector<int> ts;
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
      pbc->coggather(blaa);
      
      double dist = 0, angle=0;
      int num =0;numHbpframe=0;
  
    for(int i=0;i< int (Hm1.size());++i){
     for(int j=0;j<int (Am2.size());++j){
       num +=1;
      Vec tmp = (sys.mol(m1).pos(Hm1[i])-sys.mol(m2).pos(Am2[j])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (sys.mol(m2).pos(Am2[j])-sys.mol(m1).pos(Hm1[i])); 
       Vec tmpB = (sys.mol(m1).pos(Htom1[i])-sys.mol(m1).pos(Hm1[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
        disttotm1[num-1]+= dist; 
        angtotm1[num-1]+=angle;
        occurm1[num-1]+= 1;
	ti.push_back(time); ts.push_back(num);
        numHbpframe+=1;  
       }
      }            
     }  
    }

      int nu = num; num =0;
     for(int i=0;i<int (Hm2.size());++i){
      for(int j=0;j<int (Am1.size());++j){
       num +=1;nu+=1;  
     Vec tmp = (sys.mol(m2).pos(Hm2[i])-sys.mol(m1).pos(Am1[j])); 
      dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (sys.mol(m1).pos(Am1[j])-sys.mol(m2).pos(Hm2[i])); 
       Vec tmpB = (sys.mol(m2).pos(Htom2[i])-sys.mol(m2).pos(Hm2[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
       disttotm2[num-1]+= dist; 
       angtotm2[num-1]+= angle;
       occurm2[num-1]+= 1;
       ti.push_back(time); ts.push_back(nu);
       numHbpframe+=1;
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
    cout << setw(5) << "Number" 
         << setw(8) << "Donor" 
         << setw(16) << "Acceptor" << "       " << "D-    " << "H...    " << "A   " 
         <<" Av. Distance [nm]: "<< " Av. Angle: "<< " Occurencies: "
         << " Percentage: " << endl;
    int num=1;int nu=-1;
    vector<int> trans;
    for(int i=0;i<int (Hm1.size());++i){
     for(int j=0;j<int (Am2.size());++j){   
       nu+=1;
      if (occurm1[nu]!=0){
	cout.precision(3);
        cout.setf(ios::right, ios::adjustfield); 
          cout << setw(5) << num  
               << setw(3) << m1+1 
	       << setw(4) << sys.mol(m1).topology().resNum(Hm1[i])+1
               << setw(4) << sys.mol(m1).topology().resName(sys.mol(m1).topology().resNum(Hm1[i])) 
               << "-"
               << setw(3) << m2+1
               << setw(4) << sys.mol(m2).topology().resNum(Am2[j])+1
               << setw(4) << sys.mol(m2).topology().resName(sys.mol(m2).topology().resNum(Am2[j])) 
               << setw(4) << Htom1[i]+1 
               << setw(3) << sys.mol(m1).topology().atom(Htom1[i]).name() 
               << "-"
               << setw(4) << Hm1[i]+1 << sys.mol(m1).topology().atom(Hm1[i]).name() 
               << "..."
               << setw(4) << Am2[j]+1   << sys.mol(m2).topology().atom(Am2[j]).name() 
	       << setw(7) << (disttotm1[nu]/occurm1[nu]);
          cout.precision(5);
          cout << setw(7)  << (angtotm1[nu]/occurm1[nu])   
               << setw(6)  << (occurm1[nu])   
               << setw(6)  << (occurm1[nu]/numFrames)*100  
               << endl;
	  trans.push_back(nu+1);
	  num +=1;
      }
     }
    }
    int bla=nu; nu=-1; 
    for(int i=0;i<int (Hm2.size());++i){
     for(int j=0;j<int (Am1.size());++j){   
       nu+=1; bla+=1;
      if (occurm2[nu]!=0){
	cout.precision(3);
        cout.setf(ios::right, ios::adjustfield); 
          cout << setw(5) << num  
               << setw(3) << m2+1 
	       << setw(4) << sys.mol(m2).topology().resNum(Hm2[i])+1
               << setw(4) << sys.mol(m2).topology().resName(sys.mol(m1).topology().resNum(Hm2[i])) 
               << "-"
               << setw(3) << m1+1
               << setw(4) << sys.mol(m1).topology().resNum(Am1[j])+1
               << setw(4) << sys.mol(m1).topology().resName(sys.mol(m1).topology().resNum(Am1[j])) 
               << setw(4) << Htom2[i]+1 
               << setw(3) << sys.mol(m2).topology().atom(Htom2[i]).name() 
               << "-"
               << setw(4) << Hm2[i]+1 << sys.mol(m2).topology().atom(Hm2[i]).name() 
               << "..."
               << setw(4) << Am1[j]+1   << sys.mol(m1).topology().atom(Am1[j]).name() 
	       << setw(7) << (disttotm2[nu]/occurm2[nu]);
          cout.precision(5);
          cout << setw(7)  << (angtotm2[nu]/occurm2[nu])   
               << setw(6)  << (occurm2[nu])   
               << setw(6)  << (occurm2[nu]/numFrames)*100  
               << endl;
          trans.push_back(bla+1);
	  num +=1;
      }
     }
    }

    //get the timeseries right
    for (int i=0; i< int (ts.size()); ++i){
      for (int j=0; j < int (trans.size()); ++j){
	if (ts[i]==trans[j]){ ts[i]=j+1;}
      }
    }
    
    //smack out the timeseries       
    ofstream tfile; tfile.open("Hbts.out");
    for (int i=0; i < int (ti.size()); ++i){
      tfile << setw(5) << ti[i] << setw(6) << ts[i] << endl;
    }
    tfile.close();
     



    }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
