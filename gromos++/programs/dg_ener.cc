// dg_ener.cc This program calculates a free energy difference by reading
//            in several output files of ener, for state A and B
//            and correcting for the double counting of pairenergies
//            from the output of pairener

#include "../src/args/Arguments.h"
#include <fstream>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace args;


int main(int argc, char **argv){

  char *knowns[] = {"nfiles", "stateA", "pairA", "stateB", "pairB"};
  int nknowns = 5;

  string usage = argv[0];
  usage += "\n\t@nfiles <total number of files (A+B)\n";
  usage += "\t@stateA <energy files for state A>\n";
  usage += "\t@pairA  <pair energy for atoms included in A>\n";
  usage += "\t@stateB <energy files for state B>\n";
  usage += "\t@pairB  <pair energy for atoms included in B>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);


    // set some values
    int nfiles=0, na=0, nca=0, nb =0, ncb=0;
    double kT = 0.3*8.31441;
    
    args.check("nfiles", 1);
    
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
        iter=args.lower_bound("pairA"),
        to=args.upper_bound("pairA");
        iter!=to; ++iter){      
      files[na+nca].open((iter->second).c_str());
      nca++;
    }
    
    for(Arguments::const_iterator 
        iter=args.lower_bound("stateB"),
        to=args.upper_bound("stateB");
        iter!=to; ++iter){
      files[na+nca+nb].open((iter->second).c_str());
      nb++;
    }
    for(Arguments::const_iterator
        iter=args.lower_bound("pairB"),
        to=args.upper_bound("pairB");
        iter!=to; ++iter){      
      files[na+nca+nb].open((iter->second).c_str());
      ncb++;
    }

    double time;
    double ea[6], eb[6], eca[2], ecb[2], fdum, dh, dg, p, sum=0, ave;
    int cont=1, num=1;
    

    //first, read in the first value of time from all files
    for(int i=0;i<nfiles;i++)
      if((files[i] >> time)==0) cont=0;
    while(cont){
      for(int j=0;j<6;j++){ ea[j]=0.0; eb[j]=0.0; }
      // read energies for state A
      for(int i=0;i<na;i++)
        for(int j=0;j<6;j++){
          files[i] >> fdum;
          ea[j]+=fdum;
	}
      // read pair corrections for state A
      for(int i=0;i<nca;i++)
        files[na+i] >> eca[0] >> eca[1];
      // read energies for state B
      for(int i=0;i<nb;i++)
        for(int j=0;j<6;j++){
          files[na+nca+i] >> fdum;
	  eb[j]+=fdum;
	}
      //read pair corrections for state B
      for(int i=0;i<ncb;i++)
        files[na+nca+nb+i] >> ecb[0] >> ecb[1];

      dh = eb[4]-ecb[0]+eb[5]-ecb[1]-ea[4]+eca[0]-ea[5]+eca[1];
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




