// mcti.cc  program to analyse and prepare MCTI runs

#include "../src/args/Arguments.h"
#include <fstream>
#include <strstream>
#include <iostream>
#include <iomanip>
#include <stdlib.h>
#include <vector>
#include <cmath>

using namespace args;

class error
{
  vector<int> d_blocksize;
  vector<double> d_currav, d_blocksum, d_blockssum;
  int d_num, d_counter;

public:
  error(int num);
  ~error(){}
  void addval(double val);
  double ee();
};

void printscript(int run, double lambda, int ok, ifstream &inscript);

int main(int argc, char **argv){

  char *knowns[] = {"files", "lambda", "run", "num", "script", "shake"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@lambda <current lambda value>\n";
  usage += "\t@files <dH/dl-files>\n";
  usage += "\t@run   <number of the run>\n";
  usage += "\t@num   <number of data points>\n";
  usage += "\t@script <generalized job-script>\n";
  usage += "\t@shake <ok/fail>\n";
  

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // set some values
    args.check("lambda", 1);
    double lam=0;
    int ilam=-1;
    
    {
      Arguments::const_iterator iter=args.lower_bound("lambda");
      if(iter!=args.upper_bound("lambda"))
        lam=atof(iter->second.c_str());
    }

    args.check("run",1);
    int run=0;
    {
      Arguments::const_iterator iter=args.lower_bound("run");
      if(iter!=args.upper_bound("run"))
	run=atoi(iter->second.c_str());
    }

    args.check("num",1);
    int num=0;
    {
      Arguments::const_iterator iter=args.lower_bound("num");
      if(iter!=args.upper_bound("num"))
	num=atoi(iter->second.c_str());
    }

    args.check("script",1);
    ifstream genscript(args["script"].c_str());
    
    int ok=1;
    {
      Arguments::const_iterator iter=args.lower_bound("shake");
      if(iter!=args.upper_bound("shake"))
        if (iter->second=="fail")
          ok=0;
    }
    
    // Now start the analysis

    
    cout << "MCTI: Analysing run " << run << endl;

    // read data so far
    vector<double> l, dhdl, err;
    double ave=0, rmsd=0, ee=0;
    
    {
      ostrstream os;
      os << "data_" << run-1 << ".dat" << ends;
      ifstream data(os.str());
      double fdum[3];
      data >> fdum[0] >> fdum[1] >> fdum[2];
      while(!data.eof()){
	if(fdum[0]==lam) ilam=l.size();
	l.push_back(fdum[0]);
	dhdl.push_back(fdum[1]);
	err.push_back(fdum[2]);
	data >> fdum[0] >> fdum[1] >> fdum[2];
      }
      data.close();
      if(ilam < 0 ) throw gromos::Exception("mcti", 
				    " lambda value not in dataset.\n");
      cout << "read in " << l.size() << " lambda values from " 
           << os.str() << endl;
    }

    // if there was no shake-failure in the run, start analysis
    if(ok){
      
      // create the error estimate-class
      error::error e(num);
    
      // read in dhdl-files
      double sum=0.0, sum2=0.0, fdum;
      int nr_val=0;
      for(Arguments::const_iterator iter=args.lower_bound("files"), 
  	  to=args.upper_bound("files"); iter!=to; ++iter){
        ifstream file;
       
        file.open((iter->second).c_str());
        while((file >> fdum)!=0){
  	  sum += fdum;
	  sum2 += fdum*fdum;
	  nr_val++;
          e.addval(fdum);
	}
        file.close();
      }

      if(nr_val!=num) 
        cout << "\nWARNING: read number of dh/dl values is different from"
	     << " specified number!\n\n";
      
      //cout << sum << endl << sum2 << endl << nr_val << endl;
  
      ave=sum/nr_val;
      rmsd=sqrt((sum2 - sum*sum/nr_val)/(nr_val-1));
      ee=e.ee();

      cout << "\nFor current lambda = " << lam << endl;
      cout << "read " << nr_val << " values of dH/dl\n";
      cout << "average " << ave << "\trmsd " << rmsd 
           << "\terror estimate " << ee << endl;
      
      // update the data set
      dhdl[ilam]=ave;
      err[ilam]=ee;
    }
    else{
      cout << "\nShake failure reported in run " << run << endl;
      cout << "no new data taken into account\n";
      ee = err[ilam];
      ave = dhdl[ilam];
      // give rmsd an arbitrary value, because it is used in srand
      rmsd = err[ilam];
      
    }
    
    
    // get current estimate of the free energy difference
    // we also calculate the curvature, although this is currently not used
    double dg=0.0;
    vector<double> curv;
    {
      // create vectors for the used values of l and dhdl
      vector<double> lu, dhdlu, curvu;
      double unlikely=1e6;
      curvu.push_back(unlikely);
      for(unsigned int i=0; i< l.size();i++){
	if(err[i]!=-1){
	  lu.push_back(l[i]);
	  dhdlu.push_back(dhdl[i]);
	}
      }
      
      for(unsigned int i=0; i<lu.size()-1; i++)
        dg+=0.5*(dhdlu[i]+dhdlu[i+1])*(lu[i+1]-lu[i]);

      // calculate curvature for used points
      for(unsigned int i=1; i<lu.size()-1; i++){
	double gr1=(dhdlu[i]-dhdlu[i-1])/(lu[i]-lu[i-1]);
        double gr2=(dhdlu[i+1]-dhdlu[i])/(lu[i+1]-lu[i]);
        double crv=2*(gr2-gr1)/(lu[i+1]-lu[i-1]);
	
	curvu.push_back(crv);
      }
      curvu.push_back(curvu[lu.size()-2]);

      // extrapolate curvature
      int iu=1;
      if(curvu.size()>2){
	
        for(unsigned int i=0; i<l.size(); i++){
	  if(l[i]<=lu[0])
	    curv.push_back(curvu[1]);
	  else if(l[i]>=lu[lu.size()-1])
            curv.push_back(curvu[lu.size()-2]);
          else if(l[i]==lu[iu]){
	    curv.push_back(curvu[iu]);
	    iu++;
  	  }
	  else {
	    //extrapolate between i-1 and iu
            double crv=(curvu[iu]-curv[i-1])*(l[i]-l[i-1])/(lu[iu]-l[i-1]);
	    curv.push_back(crv+curv[i-1]);
	  }
	}
      }
      else
        for(unsigned int i=0;i<l.size();i++)
	  curv.push_back(1);

      // transform curvature into cumulative absolute values
      for(unsigned int i=1; i<l.size();i++)
	curv[i]=curv[i-1]+abs(curv[i]);
    }
    
    cout << "\nCurrent estimate of DG: " << dg << endl;

    // write the new data set to file
    ostrstream os;
    os << "data_" << run << ".dat" << ends;
    
    ofstream data(os.str());
    for(unsigned int i=0; i<l.size();i++)
      data << l[i] << "\t" << dhdl[i] << "\t" << err[i] << endl;
    data.close();
    cout << "\nWrote data-set to " << os.str() << endl << endl;
    
    // get a random lambda value

    // get a random value between 0 and l.size()
    srand(int(1000*run*abs(rmsd)));
    int r=rand();
    int itry=(r*l.size())/RAND_MAX;
    cout << "random lambda value: " << l[itry] << endl << endl;
    
    /*
    // curvature biased random lambda
    // get a random double between 0 and curv[l.size()-1]
    srand(int(1000*run*abs(rmsd)));
    int r=rand();
    double rnd=double(r*curv[l.size()-1])/double(RAND_MAX);
    //cout << "random number: " << rnd << endl;
    // find out to which lambda this belongs
    int itry=0;
    while (curv[itry]<rnd) itry++;
    cout << "curvature biased random lambda: " << l[itry] << endl<< endl;
    */
    cout << "lambda " << l[itry] << "\t" << lam << endl;
    cout << "error  " << err[itry] << "\t" << err[ilam] << endl;
    
    // determine whether we accept
    double lnew;
    
    if(err[itry]>err[ilam]||err[itry]==-1)
      lnew=l[itry];
    else{
      
      double met =exp((err[itry]-err[ilam])/2.4777);
      double test =double(rand())/double(RAND_MAX);
      cout << "\nMetropolis exp: " << met << "\trandom: " << test << endl;
      if(test<met)
        lnew=l[itry];
      else
        lnew=l[ilam];
    }
    
    cout << "\nnew lambda " << lnew << endl;
    
    // print out new input script
    printscript(run, lnew, ok, genscript);
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}


void printscript(int run, double lambda, int ok, ifstream &inscript)
{

  int ilam=int((lambda+0.0001)*100);
  ostrstream os;
  os << "jmd" << run+1 << ".sh" << ends;
  ofstream outscript(os.str());
  cout << "\nWriting jobscript " << os.str() << endl;
  outscript << "#!/bin/sh\n";
  outscript << "\n# VARIABLES SET BY MCTI\n";
  outscript << "DLAM="  << lambda << "\t # lambda value\n";
  outscript << "ILAM="  << ilam   << "\t # lambda value * 100\n";
  outscript << "OORUN=" << run-1  << "\t # two runs ago\n";
  outscript << "ORUN="  << run    << "\t # previous run\n";
  outscript << "RUN="   << run+1  << "\t # current run\n";
  outscript << "OSHK=";
  if(!ok) outscript << "fail";
  else outscript << "ok  ";
  outscript << " \t # shake failure in previous run\n";
  
  outscript << "# END OF MCTI VARIABLES\n\n";

  while(!inscript.eof()){
    char tmp[100];
    //inscript >> tmp;
    
    
    inscript.getline(tmp, 100);
    tmp[99]='\0';
    
    outscript << tmp << endl;
  }

}


error::error(int num)
{
  d_num=num;
  double blksz=2;
  
  int old=2;
  
  while(4*blksz<d_num){
    //cout << int(blksz) << endl;
    
    d_blocksize.push_back(int(blksz));
    d_currav.push_back(0.0);
    d_blocksum.push_back(0.0);
    d_blockssum.push_back(0.0);
    old=int(blksz);
    while(old==int(blksz)) blksz = blksz*1.07177; 
  }
  d_counter=0;
  
}

void error::addval(double val)
{
  int Nblocks=d_blocksize.size();
  double ave=0;
  
  for(int i=0;i<Nblocks; i++){
    d_currav[i]+=val;
    if((d_counter+1)%d_blocksize[i]==0){
      ave=d_currav[i]/d_blocksize[i];
      d_blocksum[i]+=ave;
      d_blockssum[i]+=ave*ave;
      d_currav[i]=0;
    }
  }
  d_counter++;
}

double error::ee()
{
  int Nblocks=d_blocksize.size();
  double rmsd[Nblocks];
  for(int i=0; i<Nblocks; i++){
    int Nblcki=d_num/d_blocksize[i];
    
    rmsd[i] = sqrt((d_blockssum[i] - d_blocksum[i] * d_blocksum[i]/Nblcki)/
      (Nblcki-1)/Nblcki);
    //cout << i+1 << "\t" << d_blocksize[i] << "\t" << d_blocksum[i]/Nblcki << "\t" << rmsd[i] << endl;
  }
  
  //now minimize (y[i] - rmsd[i])^2
  double t,e1, e2, v1, v2, y, dydt1=0, dydt2=0, dyda=0, t1=d_num/100, t2=d_num/1e6, a=0.5, ls=1000, lsold=0;
  while(abs((lsold-ls)/lsold)>=0.0001){
    ls=lsold;
    lsold=0;
    if(dydt1!=0) t1-=dydt1*0.001*t1/abs(dydt1);
    if(dydt2!=0) t2-=dydt2*0.001*t2/abs(dydt2);
    if(dyda !=0) a-=dyda*0.001*a/abs(dyda);
    //cout << "t1 " << t1 << " t2 " << t2 << " a " << a << endl;

    dydt1=0;
    dydt2=0;
    dyda=0;
    //cout << "ls " << ls << endl;
    
    for(int i=0; i< Nblocks; i++){
      t=d_blocksize[i];
      e1=exp(-t/t1)-1;
      e2=exp(-t/t2)-1;
      v1=2*t1*(e1*t1/t + 1);
      v2=2*t2*(e2*t2/t + 1);
      y=sqrt(a*v1+(1-a)*v2);
      dydt1 += 2*(y-rmsd[i])*(v1/t1+e1)/y;
      dydt2 += 2*(y-rmsd[i])*(v2/t2+e2)/y;
      dyda  +=   (y-rmsd[i])*(v1-v2)/y;
      lsold += (y-rmsd[i])*(y-rmsd[i]);
      //cout << "y " << y << " lsold " << lsold << endl;
      // cout <<  t << "\t" << y <<endl;
           
    }
    //cout << "dydt1 " << dydt1 << " dydt2 " << dydt2 << " dyda " << dyda << endl;
    //cout << "y lsold " << lsold << endl;
     
      
    
    
  }
  
  // make the plot
  /*  for(int i=0; i< Nblocks; i++){
      t=d_blocksize[i];
      e1=exp(-t/t1)-1;
      e2=exp(-t/t2)-1;
      v1=2*t1*(e1*t1/t + 1);
      v2=2*t2*(e2*t2/t + 1);
      y=sqrt(a*v1+(1-a)*v2);
      cout <<  t << "\t" << y <<endl;
      
    }
  */
  //error estimate
  double err=sqrt(2*(a*t1 + (1-a)*t2));
  
  return err;
}


    
