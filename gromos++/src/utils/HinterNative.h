#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif

#ifndef INCLUDED_FSTREAM
#include <fstream>
#define INCLUDED_FSTREAM
#endif

#ifndef INCLUDED_IOMANIP
#include <iomanip>
#define INCLUDED_IOMANIP
#endif

#ifndef INCLUDED_MATHH
#include <math.h>
#define INCLUDED_MATHH
#endif

#ifndef INCLUDED_IOSTREAM
#include <iostream>
#define INCLUDED_IOSTREAM
#endif


using namespace gcore;
using namespace gio;
using namespace bound;
using namespace args;
using namespace gmath;
using namespace utils;




void natIntermolecular(const Arguments &args, int m1, int m2, double time, double dt, double maxdist, double minangle){
//  read topology
   args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
 System refSys(it.system());

// parse boundary conditions
Boundary *pbc = BoundaryParser::boundary(sys, args);
  InG96 ic;

 // read reference coordinates...
   try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;

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
  

      //build donor, acceptor list based on reference structure
   double dist = 0, angle=0;
   vector<int> Hm1nat, Hm2nat, Am1nat, Am2nat, Htom1nat, Htom2nat;  

    for(int i=0;i< int (Hm1.size());++i){
     for(int j=0;j<int (Am2.size());++j){
      Vec tmp = (refSys.mol(m1).pos(Hm1[i])-refSys.mol(m2).pos(Am2[j])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (refSys.mol(m2).pos(Am2[j])-refSys.mol(m1).pos(Hm1[i])); 
       Vec tmpB = (refSys.mol(m1).pos(Htom1[i])-refSys.mol(m1).pos(Hm1[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
	 Hm1nat.push_back(Hm1[i]); Htom1nat.push_back(Htom1[i]);
         Am2nat.push_back(Am2[j]); 
       }
      }            
     }
    }

     for(int i=0;i<int (Hm2.size());++i){
      for(int j=0;j<int (Am1.size());++j){
     Vec tmp = (refSys.mol(m2).pos(Hm2[i])-refSys.mol(m1).pos(Am1[j])); 
      dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (refSys.mol(m1).pos(Am1[j])-refSys.mol(m2).pos(Hm2[i])); 
       Vec tmpB = (refSys.mol(m2).pos(Htom2[i])-refSys.mol(m2).pos(Hm2[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
       Hm2nat.push_back(Hm2[i]); Htom2nat.push_back(Htom2[i]);
         Am1nat.push_back(Am1[j]); 
      }
     }
      }
     }

  //build donor, acceptor list, spit it out
 cout << "Doing an INTERmolecular native H-bond analysis using molecule" << m1+1
      << " and " << m2+1   
      << " from the topology..." << endl;

      cout << "Starting the run..." << endl;
      vector<double> angtotm1, disttotm1, occurm1,
                     angtotm2, disttotm2, occurm2;
      for (int i=0;i< int (Hm1nat.size()* Am2nat.size());++i)
       {angtotm1.push_back(0.0);disttotm1.push_back(0.0);occurm1.push_back(0);}
      for (int i=0;i< int (Hm2nat.size()* Am1nat.size());++i)
       {angtotm2.push_back(0.0);disttotm2.push_back(0.0);occurm2.push_back(0);}
      cout << "Hydrogen Bonds to be monitored: " 
           << (Hm1nat.size())+(Hm2nat.size()) << endl;



     //prepare for run

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
      

      int num =0;numHbpframe=0;
  
    for(int i=0;i< int (Hm1nat.size());++i){
       num +=1;
      Vec tmp = (sys.mol(m1).pos(Hm1nat[i])-sys.mol(m2).pos(Am2nat[i])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (sys.mol(m2).pos(Am2nat[i])-sys.mol(m1).pos(Hm1nat[i])); 
       Vec tmpB = (sys.mol(m1).pos(Htom1nat[i])-sys.mol(m1).pos(Hm1nat[i]));
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
    

      int nu = num; num =0;
     for(int i=0;i<int (Hm2nat.size());++i){
       num +=1;nu+=1;  
     Vec tmp = (sys.mol(m2).pos(Hm2nat[i])-sys.mol(m1).pos(Am1nat[i])); 
      dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){  
       Vec tmpA = (sys.mol(m1).pos(Am1nat[i])-sys.mol(m2).pos(Hm2nat[i])); 
       Vec tmpB = (sys.mol(m2).pos(Htom2nat[i])-sys.mol(m2).pos(Hm2nat[i]));
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
    for(int i=0;i<int (Hm1nat.size());++i){ 
       nu+=1;
      if (occurm1[nu]!=0){
	cout.precision(3);
        cout.setf(ios::right, ios::adjustfield); 
          cout << setw(5) << num  
               << setw(3) << m1+1 
	       << setw(4) << sys.mol(m1).topology().resNum(Hm1nat[i])+1
               << setw(4) << sys.mol(m1).topology().resName(sys.mol(m1).topology().resNum(Hm1nat[i])) 
               << "-"
               << setw(3) << m2+1
               << setw(4) << sys.mol(m2).topology().resNum(Am2nat[i])+1
               << setw(4) << sys.mol(m2).topology().resName(sys.mol(m2).topology().resNum(Am2nat[i])) 
               << setw(4) << Htom1nat[i]+1 
               << setw(3) << sys.mol(m1).topology().atom(Htom1nat[i]).name() 
               << "-"
               << setw(4) << Hm1nat[i]+1 << sys.mol(m1).topology().atom(Hm1nat[i]).name() 
               << "..."
               << setw(4) << Am2nat[i]+1   << sys.mol(m2).topology().atom(Am2nat[i]).name() 
	       << setw(7) << (disttotm1[nu]/occurm1[nu]);
          cout.precision(5);
          cout << setw(7)  << (angtotm1[nu]/occurm1[nu])   
               << setw(10)  << (occurm1[nu])   
               << setw(10)  << (occurm1[nu]/numFrames)*100  
               << endl;
	  trans.push_back(nu+1);
	  num +=1;
      }
     }
    
    int bla=nu; nu=-1; 
    for(int i=0;i<int (Hm2nat.size());++i){
       nu+=1; bla+=1;
      if (occurm2[nu]!=0){
	cout.precision(3);
        cout.setf(ios::right, ios::adjustfield); 
          cout << setw(5) << num  
               << setw(3) << m2+1 
	       << setw(4) << sys.mol(m2).topology().resNum(Hm2nat[i])+1
               << setw(4) << sys.mol(m2).topology().resName(sys.mol(m2).topology().resNum(Hm2nat[i])) 
               << "-"
               << setw(3) << m1+1
               << setw(4) << sys.mol(m1).topology().resNum(Am1nat[i])+1
               << setw(4) << sys.mol(m1).topology().resName(sys.mol(m1).topology().resNum(Am1nat[i])) 
               << setw(4) << Htom2nat[i]+1 
               << setw(3) << sys.mol(m2).topology().atom(Htom2nat[i]).name() 
               << "-"
               << setw(4) << Hm2nat[i]+1 << sys.mol(m2).topology().atom(Hm2nat[i]).name() 
               << "..."
               << setw(4) << Am1nat[i]+1   << sys.mol(m1).topology().atom(Am1nat[i]).name() 
	       << setw(7) << (disttotm2[nu]/occurm2[nu]);
          cout.precision(5);
          cout << setw(7)  << (angtotm2[nu]/occurm2[nu])   
               << setw(10)  << (occurm2[nu])   
               << setw(10)  << (occurm2[nu]/numFrames)*100  
               << endl;
          trans.push_back(bla+1);
	  num +=1;
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
