// dg_ener.cc This program calculates a free energy difference by reading
//            in the output files of ener, for state A and B

#include "../src/args/Arguments.h"
#include "../src/gmath/physics.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"temp", "stateA", "stateB"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@temp <temperature for perturbation>\n";
  usage += "\t@stateA <energy file for state A>\n";
  usage += "\t@stateB <energy file for state B>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set some values
    args.check("temp", 1);
    double temp=0;
    
    {
      Arguments::const_iterator iter=args.lower_bound("temp");
      if(iter!=args.upper_bound("temp"))
        temp=atof(iter->second.c_str());
    }
    ifstream stateA, stateB;

    //open files
    {
      Arguments::const_iterator iter=args.lower_bound("stateA");
      if(iter!=args.upper_bound("stateA"))
        stateA.open((iter->second).c_str());
      else
        throw gromos::Exception("dg_ener", "energy file for state A missing\n");
      iter=args.lower_bound("stateB");
      if(iter!=args.upper_bound("stateB"))
	stateB.open((iter->second).c_str());
      else
	throw gromos::Exception("dg_ener", "energy file for state B missing\n");
    }
    cout.precision(12);
    
    //print title
    cout << "# Time"
	 << setw(12) << "DE_tot"
	 << setw(12) << "probability"
	 << setw(12) << "DG_AB" 
	 << endl;
    
    string sdum;
    
    double time=0;
    int cont=1, num =0;
    // first, get rid of possible titles in the files
    while(cont==1){
      if((stateA >> sdum)==0) cont=0;
      if (sdum!="#") {
	time = atof(sdum.c_str());
	cont=0;
      } else 
	stateA >> sdum >> sdum >> sdum >> sdum;
	
      if((stateB >> sdum)==0) cont=0;
      if(sdum!="#") {
	time = atof(sdum.c_str());
	cont=0;
      }
      else
	stateB >> sdum >> sdum >> sdum >> sdum;
      
    }
 
    // now, read in the total energies
    cont =1;
    double covA, nbA, totA, covB, nbB, totB;
    double dh, sum=0, p, dg, ave;
    
    
    num =1;
    //first store them in a vector and calculate the average
    vector<double> delta_v;
    vector<double> t;
    t.push_back(time);
    
    double ave_v=100000.0;
    
    while(cont==1){
      stateA >> covA >> nbA >> totA;
      stateB >> covB >> nbB >> totB;
      dh=totB-totA;
      delta_v.push_back(dh);
      if(dh<ave_v) ave_v=dh;
      

      num++;
      
      // now, we try to read in the next time
      if ((stateA >> sdum)==0) cont=0;
      if (sdum=="#") cont=0;
      if ((stateB >> sdum)==0) cont=0;
      if (sdum=="#") cont =0;
      time = atof(sdum.c_str());
      t.push_back(time);
      
    }
    //ave_v/=num;
    
    
    // now loop over all values that were read in 
    for(unsigned int i=0; i<delta_v.size(); i++){
      dh=delta_v[i]-ave_v;
      p = exp(-dh/(BOLTZ*temp));
      sum += p;
      ave = sum/i;
      dg = ave_v - BOLTZ*temp*log(ave);
      
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << t[i]
           << setw(12) << delta_v[i]
           << setw(12) << p
           << setw(12) << dg
           << endl;

    }
    cout << "# substracted value " << ave_v << endl;
    cout << "# final result: " << dg << endl;
    
    stateA.close();
    stateB.close();
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




