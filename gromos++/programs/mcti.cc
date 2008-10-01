// mcti.cc  program to analyse and prepare MCTI runs

#include <cassert>

#include "../src/args/Arguments.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace std;
using namespace args;

class stat
{
  vector<int> d_blocksize;
  vector<double> d_ener, d_vals;
  int d_num, d_counter;
  double d_ave,d_rmsd, d_kT;
  
public:
  stat(int num);
  ~stat(){}
  void addval(double val, double ener);
  void setKT(double kT);
  double rmsd();
  double ave();
  double subave(int b, int e);
  double Z(double &min, int b, int e);
  double ee();
  int n();
  
};

void printscript(int run, double lambda, int ok, ifstream &inscript);

int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "files" << "lambda" << "run" << "num" << "script" << "shake";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@lambda <current lambda value>\n";
  usage += "\t@files <dH/dl-files>\n";
  usage += "\t@run   <number of the run>\n";
  usage += "\t@num   <number of data points>\n";
  usage += "\t@script <generalized job-script>\n";
  usage += "\t@shake <ok/fail>\n";
  

  try{
    Arguments args(argc, argv, knowns, usage);

    // set some values
    double kT = 0.29815*8.31441;
    
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
      ostringstream os;
      os << "data_" << run-1 << ".dat";
      ifstream data(os.str().c_str());
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
      stat s(num);
      s.setKT(kT);
      
      // read in dhdl-files
      double fdum1, fdum2;
      
      int nr_val=0;
      for(Arguments::const_iterator iter=args.lower_bound("files"), 
  	  to=args.upper_bound("files"); iter!=to; ++iter){
        ifstream file;
       
        file.open((iter->second).c_str());
        while((file >> fdum1)!=0){
          file >> fdum2;
          s.addval(fdum2,fdum1);
	  nr_val++;
	}
        file.close();
      }
      if(nr_val!=num) 
        cout << "\nWARNING: read number of dh/dl values is different from"
	     << " specified number!\n\n";
      
      ave=s.ave();
      rmsd=s.rmsd();
      ee=s.ee();

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
    ostringstream os;
    os << "data_" << run << ".dat";
    
    ofstream data(os.str().c_str());
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
  ostringstream os;
  os << "jmd" << run+1 << ".sh";
  ofstream outscript(os.str().c_str());
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

stat::stat(int num)
{
  d_num=num;
  double blksz=50;
  int old=2;
  while(4*blksz<d_num){
    d_blocksize.push_back(int(blksz));
    old=int(blksz);
    while(old==int(blksz)) blksz = blksz*1.07177; 
  }
  for(int i=0;i<num;i++){
    d_ener.push_back(0.0);
    d_vals.push_back(0.0);
  }
  d_counter=0;
  d_kT=0.29815*8.31441;
  d_ave=1e6;
  d_rmsd=1e6;
}

void stat::addval(double val, double ener)
{
  d_vals[d_counter]=val;
  d_ener[d_counter]=ener;
  d_counter++;
}
void stat::setKT(double kT)
{
  d_kT=kT;
}


double stat::ave()
{
  if(d_ave==1e6)
    d_ave = this->subave(0,d_num);
  return d_ave;
}

double stat::subave(int b, int e)
{
  double ave=0;
  double min=1e6;
  //calculate min and Z
  double Z=this->Z(min, b, e);
  
  //calculate the average
  for(int i=b;i<e;i++){
    ave+=d_vals[i]*exp(-(d_ener[i]-min)/d_kT)/Z;
  }
  return ave;
}

double stat::Z(double &min, int b, int e)
{
  min=1e6;
  for(int i=b; i<e; i++)
    if(d_ener[i]<min)min=d_ener[i];
  double Z=0;
  for(int i=b; i<e; i++)
    Z+=exp(-(d_ener[i]-min)/d_kT);
  return Z;
}

double stat::rmsd()
{
  if(d_rmsd==1e6){
  
    double rmsd=0;
    double min=1e6;
    //get min and Z
    double Z=this->Z(min,0,d_counter);
    double ave=this->ave();
  
    //calculate the average of the square of the difference
    for(int i=0; i<d_counter; i++){
      rmsd+= (d_vals[i]-ave)*(d_vals[i]-ave)*exp(-(d_ener[i]-min)/d_kT)/Z;
    }
    d_rmsd=sqrt(rmsd);
  }
  
  return d_rmsd;
    
}

int stat::n()
{
  return d_counter;
}

double stat::ee()
{
  int Nblocks=d_blocksize.size();
  double rmsd2, ave=0, min, en, Z, Zblocks;
  double runave=this->ave();
  double runrmsd=this->rmsd();
  std::vector<double> fit(Nblocks), x(Nblocks);
  
  for(int j=0; j<Nblocks; j++){
    int Nblcki=d_num/d_blocksize[j];
    std::vector<double> Eave(Nblcki);

    // for every block, calculate the average energy
    for(int i=0; i<Nblcki; i++){
      min=1e6;
      Eave[i]=0;
      Z=this->Z(min, i*d_blocksize[j],(i+1)*d_blocksize[j]);
      
      for(int k=0; k<d_blocksize[j]; k++){
	en=d_ener[i*d_blocksize[j]+k];
	Eave[i]+=en*exp(-(en-min)/d_kT);
      }
      Eave[i]/=Z;
    }
    
    // now calculate the 'partition function' over the blocks
    min=1e6;
    for(int i=0; i<Nblcki; i++) if(Eave[i]<min) min=Eave[i];
    Zblocks=0.0;
    for(int i=0; i<Nblcki; i++)
      Zblocks += exp(-(Eave[i]-min)/d_kT);
    
    // The rmsd of the property we are interested in, weighted by the
    // average energy of the blocks
    rmsd2=0;
    for(int i=0; i<Nblcki; i++){
      ave = this->subave(i*d_blocksize[j],(i+1)*d_blocksize[j]);
      rmsd2+=(ave-runave)*(ave-runave)*exp(-(Eave[i]-min)/d_kT);
    }
    rmsd2/=Zblocks;
    fit[j]=d_blocksize[j]*rmsd2/runrmsd/runrmsd;
    x[j]=1.0/d_blocksize[j];
    
  }
  double sx=0, sf=0,sfx=0,sxx=0;
  for(int i=0; i<Nblocks;i++){
    sx+=x[i];
    sf+=fit[i];
    sfx+=x[i]*fit[i];
    sxx+=x[i]*x[i];
  }
  
  double a, b;
  a=(sf*sx/Nblocks-sfx)/(sx*sx/Nblocks-sxx);
  b = (sf - a*sx)/Nblocks;
  //  for(int i=0; i<Nblocks; i++)
  //    cout << x[i] << "\t" << fit[i] << "\t" << a*x[i]+b << endl;
  
  double error=sqrt(b/d_counter)*runrmsd;
  return error;
}
