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




void SoluSolv(const Arguments &args, int m1, double time, double dt, double maxdist, double minangle){
//  read topology
   args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
Boundary *pbc = BoundaryParser::boundary(sys, args);
Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
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

   for (int j=0;j<sys.sol(0).topology().numAtoms();++j){  
       if (sys.sol(0).topology().atom(j).mass()==1.00800){
	 Hm2.push_back(j);Neighbours neigh(sys,0,j,0); Htom2.push_back(neigh[0]);
       }
      else if (sys.sol(0).topology().atom(j).mass()==15.99940){
	Am2.push_back(j);
      }  
        else if (sys.sol(0).topology().atom(j).mass()==14.00670){
	Am2.push_back(j); 
	  }
      else if (sys.sol(0).topology().atom(j).mass()==32.06000){
	Am2.push_back(j); 
	  } 
   }
  
 
  //build donor, acceptor list, spit it out
      cout << "Doing a SOLUTE <-> SOLVENT H-bond analysis using molecule" << m1+1 
           << " from the topology..." << endl;

      cout << "Starting the run..." << endl;
      vector<double> angtotm1, disttotm1, occurm1,
                     angtotm2, disttotm2, occurm2;
      for (int i=0;i< int (Hm1.size()* Am2.size());++i)
       {angtotm1.push_back(0.0);disttotm1.push_back(0.0);occurm1.push_back(0);}
      for (int i=0;i< int (Hm2.size()* Am1.size());++i)
       {angtotm2.push_back(0.0);disttotm2.push_back(0.0);occurm2.push_back(0);}
      cout << "Potential Hydrogen Bonds to be monitored per solvent molecule: " 
           << (Hm1.size()* Am2.size())+(Hm2.size()* Am1.size()) << endl;



      ofstream tsh; tsh.open("Hbnumts.out");
      int numFrames = 0, numHbpframe=0;
  
      vector<double> ti; vector<int> ts;
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
	numFrames++;
        ic.select("ALL");

      ic >> sys;
       
      (*pbc.*gathmethod)();

      double dist = 0, angle=0;
      int num = 0;numHbpframe=0;
      // do protein -> water
    for (int a=0,na=sys.sol(0).topology().numAtoms(), tna=sys.sol(0).numCoords();a<tna;a+=na){
      num=0;
     for(int i=0;i< int (Hm1.size());++i){
      for(int j=0;j<int (Am2.size());++j){
       num +=1; 
       Vec tmp = (sys.mol(m1).pos(Hm1[i])-sys.sol(0).pos(a+Am2[j])); 
       dist = sqrt(tmp.dot(tmp));
       if (dist <= maxdist){  
        Vec tmpA = (sys.sol(0).pos(a+Am2[j])-sys.mol(m1).pos(Hm1[i])); 
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
      
    
      int nu = num; int numm = 0;
      // do water -> protein 
     for(int i=0;i<int (Hm2.size());++i){
      for(int j=0;j<int (Am1.size());++j){
       numm +=1;nu+=1;  
       Vec tmp = (sys.sol(0).pos(a+Hm2[i])-sys.mol(m1).pos(Am1[j])); 
       dist = sqrt(tmp.dot(tmp));
       if (dist <= maxdist){  
        Vec tmpA = (sys.mol(m1).pos(Am1[j])-sys.sol(0).pos(a+Hm2[i])); 
        Vec tmpB = (sys.sol(0).pos(a+Htom2[i])-sys.sol(0).pos(a+Hm2[i]));
        angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
        if (angle >= minangle){
	 disttotm2[numm-1]+= dist; 
	 angtotm2[numm-1]+= angle;
	 occurm2[numm-1]+= 1;
         ti.push_back(time); ts.push_back(nu);
         numHbpframe+=1;
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
      
      
        ic.close();
      }
        tsh.close();
      
	//place the averaging in here
      cout << endl;
    cout << "Statistics of the run:" << endl;
    cout << setw(5)  << "Number"
         << setw(10) << "Donor" << "-"
         << setw(11) << "Acceptor"
         << setw(8)  << "D" << " - "
         << setw(8)  << "H" << " ... "
         << setw(8)  << "A"
         << setw(16) << "Av. Dist. [nm]"
         << setw(10) << "Av. Angle"
         << setw(12) << "Occurencies"
         << setw(12) << "Percentage" << endl;
         
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
               << setw(4) << Htom1[i]+1 
               << setw(4) << sys.mol(m1).topology().atom(Htom1[i]).name() 
               << " - "
               << setw(4) << Hm1[i]+1 
               << setw(4) << sys.mol(m1).topology().atom(Hm1[i]).name() 
               << " ... "
               << setw(4) << Am2[j]+1   
               << setw(4) << sys.sol(0).topology().atom(Am2[j]).name() 
	       << setw(16) << (disttotm1[nu]/occurm1[nu]);
          cout.precision(5);
          cout << setw(10)  << (angtotm1[nu]/occurm1[nu])   
               << setw(12)  << (occurm1[nu])   
               << setw(12)  << (occurm1[nu]/numFrames)*100  
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
               << setw(3) << m1+1
               << setw(4) << sys.mol(m1).topology().resNum(Am1[j])+1
               << setw(4) << sys.mol(m1).topology().resName(sys.mol(m1).topology().resNum(Am1[j])) 
               << setw(4) << Htom2[i]+1 
               << setw(4) << sys.sol(0).topology().atom(Htom2[i]).name() 
               << " - "
               << setw(4) << Hm2[i]+1 
               << setw(4) << sys.sol(0).topology().atom(Hm2[i]).name() 
               << " ... "
               << setw(4) << Am1[j]+1   
               << setw(4) << sys.mol(m1).topology().atom(Am1[j]).name() 
	       << setw(16) << (disttotm2[nu]/occurm2[nu]);
          cout.precision(5);
          cout << setw(10)  << (angtotm2[nu]/occurm2[nu])   
               << setw(12)  << (occurm2[nu])   
               << setw(12)  << (occurm2[nu]/numFrames)*100  
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
