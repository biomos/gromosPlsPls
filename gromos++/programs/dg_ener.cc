// dg_ener.cc This program calculates a free energy difference by reading
//            in several output files of ener, for state A and B

#include "../src/args/Arguments.h"
#include <fstream>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"nfiles", "stateA", "stateB"};
  int nknowns = 3;

  string usage = argv[0];
  usage += "\n\t@nfiles <total number of files (A+B)\n";
  usage += "\t@stateA <energy files for state A>\n";
  usage += "\t@atateB <energy files for state B>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set some values
    int nfiles=0, na=0, nb=0;
    double kT = 0.3*8.31441;
    
    
    {
      Arguments::const_iterator iter=args.lower_bound("nfiles");
      if(iter!=args.upper_bound("nfiles"))
        nfiles=atoi(iter->second.c_str());
    }
    ifstream files[nfiles];

    //open files
    for(Arguments::const_iterator 
        iter=args.lower_bound("stateA"),
        to=args.upper_bound("stateA");
        iter!=to; ++iter){
      files[na].open((iter->second).c_str());
      na++;
    }
    for(Arguments::const_iterator 
        iter=args.lower_bound("stateB"),
        to=args.upper_bound("stateB");
        iter!=to; ++iter){
      files[na+nb].open((iter->second).c_str());
      nb++;
    }

    double time;
    double ea[6], eb[6], fdum, dh, dg, p, sum=0, ave;
    int cont=1, num=1;
    

    //first, read in the first value of time from all files
    for(int i=0;i<nfiles;i++)
      if((files[i] >> time)==0) cont=0;
    while(cont){
      for(int j=0;j<6;j++){ ea[j]=0.0; eb[j]=0.0; }
      
      for(int i=0;i<na;i++)
        for(int j=0;j<6;j++){
          files[i] >> fdum;
          ea[j]+=fdum;
	}
      for(int i=na;i<nfiles;i++)
        for(int j=0;j<6;j++){
          files[i] >> fdum;
	  eb[j]+=fdum;
	}
      dh = eb[4]+eb[5]-ea[4]-ea[5];
      p = exp(-dh/kT);
      sum +=p;
      ave = sum/num;
      num++;
      dg = -kT*log(ave);
      
      cout.precision(5);
      cout.setf(ios::right, ios::adjustfield);
      cout << setw(6) << time
           << setw(12) << dh
           << setw(12) << p
           << setw(12) << dg
           << endl;

      for(int i=0;i<nfiles;i++)
        if((files[i] >> time)==0) cont=0;
    }
    

    for(int i=0;i<nfiles;i++)
      files[i].close();
      
    
    
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}




