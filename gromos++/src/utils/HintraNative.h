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




void natIntramolecular(const Arguments &args, int molecule, double time, double dt, double maxdist, double minangle){
//  read topology
   args.check("topo",1);
  InTopology it(args["topo"]);
  System sys(it.system());
   System refSys(it.system());

   //parse boundary conditions
Boundary *pbc = BoundaryParser::boundary(sys, args);
  InG96 ic;

  //read reference coordinates
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

  //build donor, acceptor list based on reference system
  double dist = 0, angle=0;
  vector<int> refHatoms, refAcceptors;

    for(int i=0;i< int (Hatoms.size());++i){
     for(int j=0;j<int (Acceptors.size()) ;++j){   
     if (hto[i] != Acceptors[j]){
      Vec tmp = (refSys.mol(molecule).pos(Hatoms[i])-refSys.mol(molecule).pos(Acceptors[j])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){
       Vec tmpA = (refSys.mol(molecule).pos(Acceptors[j])-refSys.mol(molecule).pos(Hatoms[i])); 
       Vec tmpB = (refSys.mol(molecule).pos(hto[i])-refSys.mol(molecule).pos(Hatoms[i]));
       angle = acos((tmpA.dot(tmpB))/(sqrt(tmpA.dot(tmpA))*(sqrt(tmpB.dot(tmpB)))))*180/3.1416;
       if (angle >= minangle){
	 refHatoms.push_back(Hatoms[i]);
         refAcceptors.push_back(Acceptors[j]);
       }
      }
     }
     }
    }
    


  
 
  //build donor, acceptor list, spit it out
   if (Hatoms.size() == 0) {throw gromos::Exception("Hbond", "No Donors found!\n");}
   if (Acceptors.size() == 0) {throw gromos::Exception("Hbond", "No Acceptors found!\n");}


     int num =1;
     for(int i=0;i< int (refHatoms.size());++i){
	  num +=1;
     }
    
                  
cout << "Doing an INTRAmolecular NATIVE H-bond analysis using molecule" << molecule+1 
     << " from the topology..." << endl;

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
      

      
      numHbpframe=0;
      num =0; 
    for(int i=0;i< int (refHatoms.size());++i){
       num += 1;
      Vec tmp = (sys.mol(molecule).pos(refHatoms[i])-sys.mol(molecule).pos(refAcceptors[i])); 
       dist = sqrt(tmp.dot(tmp));
      if (dist <= maxdist){
	Neighbours nei(sys,molecule,refHatoms[i]);
       Vec tmpA = (sys.mol(molecule).pos(refAcceptors[i])-sys.mol(molecule).pos(refHatoms[i])); 
       Vec tmpB = (sys.mol(molecule).pos(nei[0])-sys.mol(molecule).pos(refHatoms[i]));
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
     for(int i=0;i< int (refHatoms.size());++i){
	  k +=1;
         if (occur[k-1]!=0){
        Neighbours ne(sys,molecule,refHatoms[i]);
        cout.setf(ios::floatfield, ios_base::fixed);
        cout.setf(ios::right, ios::adjustfield);
        cout.precision(3); 
	cout << setw(5) << num << setw(7) << k 
             << setw(4) << sys.mol(molecule).topology().resNum(refHatoms[i])+1
             << sys.mol(molecule).topology().resName(sys.mol(molecule).topology().resNum(refHatoms[i])) 
             << setw(2) << "-"
             << setw(4) << sys.mol(molecule).topology().resNum(refAcceptors[i])+1
	     << sys.mol(molecule).topology().resName(sys.mol(molecule).topology().resNum(refAcceptors[i])) 
             << setw(6) << ne[0]+1 << setw(4) << sys.mol(molecule).topology().atom(ne[0]).name() << "-"
             << setw(6) << refHatoms[i]+1  						    
             << setw(4) << sys.mol(molecule).topology().atom(refHatoms[i]).name() << "-"
	     << setw(6) << refAcceptors[i]+1 
             << setw(4) << sys.mol(molecule).topology().atom(refAcceptors[i]).name() 
	     << setw(7) << (disttot[k-1]/occur[k-1]);

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
	  per90 << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
        else if (per > 70){
          per70 << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
        else if (per > 50){
          per50 << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
        else if (per > 30){
          per30 << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
        else if (per > 10){
          per10 << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
        else if (per > 0){
          per0  << sys.mol(molecule).topology().resNum(refHatoms[i])+1 << " " <<  sys.mol(molecule).topology().resNum(refAcceptors[i])+1 << endl;}
         
          cnum.push_back(num);
          ck.push_back(k);
          num +=1;
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
