#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"
#include "../src/gio/InTopology.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gmath/Stat.h"
#include "../src/gmath/Physics.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <vector>
#include <cmath>

using namespace args;
using namespace gio;
using namespace gcore;

using namespace std;


int main(int argc, char **argv){

  Argument_List knowns; 
  knowns << "topo" << "time" << "files" << "temp" << "volume" << "inttime" 
         << "range";

  string usage = "# " + string(argv[0]);
  usage += "\n\t@topo  <topology>\n";
  usage += "\t@time    <t and dt>\n";
  usage += "\t@files   <pressure files>\n";
  usage += "\t@temp    <temperature>\n";
  usage += "\t@volume  <volume>\n";
  usage += "\t@inttime <integration timestep>\n";
  usage += "\t@range   <determining range>\n";
  

  try{
    Arguments args(argc, argv, knowns, usage);

    // get simulation time
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

    // get the range
    double tbegin=0, tend=-1;
    {
      Arguments::const_iterator iter=args.lower_bound("range");
      if(iter!=args.upper_bound("range")){
        tbegin=atof(iter->second.c_str())-time;
        ++iter;
      }
      if(iter!=args.upper_bound("range"))
        tend=atof(iter->second.c_str())-time;
    }
    
    // and the temperature, volume and integration time
    double temp=atof(args["temp"].c_str());
    double vol=atof(args["volume"].c_str());
    double intt=atof(args["inttime"].c_str());
    
    // setup the virial, kinetic energy and pressure arrays
    double vir[6], kin[6], press[6];

    // and the statistical information
    gmath::Stat<double> p[3];

    // set some constants
    double fac=vol/(2.0*gmath::boltz*temp);
    double conv=0.0016605655;
    string label[3]={"xy","xz","yz"};

    // write a header
    cout << "#" << setw(9) << "time";
    for(int k=0; k<3; k++)
      cout << setw(14) << "<DG("+label[k]+")>"
	   << setw(14) << "<DG("+label[k]+")^2>";
    cout << endl;
    
    // loop over the files
    Ginstream gin;
    string line;
    vector<string> buffer;
    
    Arguments::const_iterator iter=args.lower_bound("files"),
      to=args.upper_bound("files");
    for(;iter!=to; ++iter){
      gin.open((iter->second).c_str());
      //cout << "opened file" << endl;
      gin.getline(line);
      
      while(!gin.stream().eof()){
	while(line!="VISCOSITY"){
	  gin.getline(line);
	  
	  if(gin.stream().eof())
	    throw gromos::Exception("visco_ana", "couldn't find VISCOSITY block in file " + iter->second);
	  
	}
	
	//cout << "found the block " << line << endl;
	//gin.getblock(buffer);
	
	//cout << "block size " << buffer.size() << endl;
	//for(unsigned int j=0; j< buffer.size(); j++){
	//  cout << buffer[j] << endl;
	//}
	
	//while(dum!="ENERGY") gin >> dum;
	//string s=gio::concatenate(buffer.begin(), buffer.end()-1, s);
	//cout << "string :\n" << s << endl;
	
	//istringstream is(s);
	
	for(int i=0; i<6; i++){
	  gin.getline(line); vir[i]=atof(line.c_str());
	}
	
	for(int i=0; i<6; i++){
	  gin.getline(line); kin[i]=atof(line.c_str());
	}
	
	for(int i=0; i<6; i++){
	  gin.getline(line); press[i]=atof(line.c_str());
        
	
	  //cout << "read the numbers " << vir[i] << " " << kin[i] << " " << press[i]<<  endl;
	}
	
        p[0].addval((press[0]+press[3])/2);
	p[1].addval((press[1]+press[4])/2);
	p[2].addval((press[2]+press[5])/2);
	//cout << "added the numbers" << endl;
	gin.getline(line);
	gin.getline(line);
      }
      gin.close();
      
    }

    //cout << "finished reading the file correctly" << endl;
    
    if(tend==-1)tend=p[0].n()*dt;
    
    //now, do the integration
    //calculate how many steps we can do
    int num_per_int=int(intt/dt);
    //cout << intt << " " << dt << " " << num_per_int << endl;
    int max_time_num=int(tend/dt);
    int begin_time_num=int(tbegin/dt);
    //cout << max_time_num << " " << begin_time_num<< endl;
    
    
    int max_int=max_time_num/num_per_int;
    //cout << p[0].n() << " " << max_int << endl;
    //set some values for the regression
    double sx[3]={0,0,0},sy[3]={0,0,0},sxx[3]={0,0,0},sxy[3]={0,0,0};
    int nreg=0;
    
    //cout << "cAlculated things" << endl;
    
    for(int i=1; i<max_int; i++){
      //cout << "loop over integration points " << i <<  endl;
      
      double sum[3]={0.0,0.0,0.0};
      double ssum[3]={0.0,0.0,0.0};
      double integrant=0;
      
      //we start every i*num_per_int/2
      int start=i*num_per_int/2;
      int count=0;
      
      for(int j=0; j+i*num_per_int<p[0].n(); j+=start, count++){
	
	//cout << "loop over j " << j << endl;
      
	for(int k=0; k<3; k++){
	  integrant=p[k].subave(j,j+i*num_per_int)*i*num_per_int*dt;
	  sum[k]+=integrant;
	  ssum[k]+=integrant*integrant;
	}
      }
      
      if(i>=begin_time_num/num_per_int){
	for(int k=0; k<3; k++){
	  
	  sx[k]+=i*intt;
	  sy[k]+=ssum[k]/count;
	  sxx[k]+=i*i*intt*intt;
	  sxy[k]+=i*intt*ssum[k]/count;
	  
	}
	nreg++;
      }
      
      cout << setw(10) << i*intt;
      for(int k=0; k<3; k++)
	cout << setw(14) << sum[k]/count
	     << setw(14) << ssum[k]/count;
      cout << endl;
    }
    double a[3], b[3];

    cout << endl;
    cout << "# Least square fit over time = " << time+tbegin 
	 << " to " << time+tend << " ps" << endl;
    for(int k=0; k<3; k++){
      
      a[k] = (sxy[k] - sx[k]*sy[k]/nreg)/(sxx[k]-sx[k]*sx[k]/nreg);
      b[k] = -(a[k]*sx[k] -sy[k])/nreg;
      cout <<endl;
      
      cout << "# <DG(" << label[k] << ")^2> = " << b[k] << " + " << a[k] 
	   << " * t" << endl;
      cout << "# viscosity (" << label[k] << ") = " << a[k]*fac 
	   << " u/nm/ps = " << a[k]*fac*conv << " g/m/ps (= cP)" << endl;
    }
    cout << endl;
    double ave=(a[0]+a[1]+a[2])*fac/3;
    
    cout << "# average viscosity: " << ave << " u/nm/ps" << endl;
    cout << "# average viscosity: " << ave * conv 
	 << " g/m/s (= cP)" << endl;
    
    
    
  }
  
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}





