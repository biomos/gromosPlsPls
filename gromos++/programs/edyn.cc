// essential dynamics edyn.cc 
//--mika

#include <cassert>

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/AtomTopology.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"
#include "../src/gmath/Matrix.h"
#include "../src/gmath/Vec.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <functional>
#include <iostream>

using namespace std;
using namespace gmath;
using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;

void writePdb(const char *name, int dim, vector<int> atoml, System &sys);

int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "class", "pbc", "ref", "mol"};
  int nknowns = 6;
  int mol = 0;
  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc <boundary type>\n";
  usage += "\t@mol <molecules to be considered.>\n";
  usage += "\t@class <classes of atoms to consider>\n";
  usage += "\t@ref <reference coordinates>\n";
  usage += "\t@traj <trajectory files>\n";


  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // read topology
    InTopology it(args["topo"]);
    System refSys(it.system());
    
    // read reference coordinates...
    InG96 ic;

    try{
      args.check("ref",1);
      ic.open(args["ref"]);
    }
    catch(const Arguments::Exception &){
      args.check("traj",1);
      ic.open(args["traj"]);
    }
    ic >> refSys;
    ic.close();


    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(refSys, args);
    //parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);
    // gather reference system
    (*pbc.*gathmethod)();
    delete pbc;
    
    Reference ref(&refSys);

    // Adding references
    int added=0;
    // which molecules considered?
    vector<int> mols;
    if(args.lower_bound("mol")==args.upper_bound("mol"))
      for(int i=0;i<refSys.numMolecules();++i)
	mols.push_back(i);
    else
      for(Arguments::const_iterator it=args.lower_bound("mol");
	  it!=args.upper_bound("mol");++it){
	if(atoi(it->second.c_str())>refSys.numMolecules())
	  throw Arguments::Exception(usage);
	mols.push_back(atoi(it->second.c_str())-1);
      }
    // add classes for fit and essential dynamics
    vector<string> name;
    for(Arguments::const_iterator it=args.lower_bound("class");
	it != args.upper_bound("class"); ++it){
     name.push_back(it->second);
      for(vector<int>::const_iterator mol=mols.begin();
	  mol!=mols.end();++mol)
	ref.addClass(*mol,it->second);
      added=1;
    }

    // System for calculation
    System sys(refSys);

    // get atoms from class
    vector<int> atomlist;
    int tmp = name.size();
    for (int i=0;i<tmp;++i){
      ref.makePosList(sys, mols[0], name[i], atomlist);
    } 
   
     // sort atom numbers
    std::sort(atomlist.begin(), atomlist.end());
  
    // did we add anything at all?
    if(!added)
      throw Arguments::Exception(usage);

    // Parse boundary conditions for sys
    pbc = BoundaryParser::boundary(sys, args);

    RotationalFit rf(&ref);


 int numFrames = 0;
 int size = atomlist.size();
 std::vector<double> avpos(size*3);
 for (int i = 0;i<size*3;++i){avpos[i]=0;}
 std::vector<double> pvec(size*3);
 Matrix cov(size*3,size*3,0.0);

 cout << "reading trajectory..."<< endl;
 for(Arguments::const_iterator iter=args.lower_bound("traj");
  iter!=args.upper_bound("traj"); ++iter){
  ic.open(iter->second);
    // loop over all frames
  while(!ic.eof()){
   numFrames++;
   ic >> sys;	
   (*pbc.*gathmethod)();
   rf.fit(&sys);

    // calculate average positions, put coords into one array
   for (int i=0, x=0; i<size;++i,x+=2) {
     int j = atomlist[i];
    pvec[i+x]=sys.mol(mol).pos(j)[0];
    pvec[i+1+x]=sys.mol(mol).pos(j)[1];
    pvec[i+2+x]=sys.mol(mol).pos(j)[2];
    avpos[i+x]+=pvec[i+x];
    avpos[i+1+x]+= pvec[i+1+x];
    avpos[i+2+x]+= pvec[i+2+x];
   }  

   //build one triangle of the covariance matrix
   for (int ii=0; ii<size*3;++ii){
    for (int jj=0; jj<=ii;++jj){
     cov(ii,jj)=cov(ii,jj)+(pvec[ii]*pvec[jj]);
    }
   }
	   
  }
  ic.close();
 }
   
  //average positions, write them out
    for (int i=0, x=0, y=0,z=0; i<size*3;++i) {
     y=atomlist[z];
     avpos[i] = avpos[i]/numFrames;
     sys.mol(0).pos(y)[x]=avpos[i];
    x+=1; if (x == 3){x=0;z+=1;}
    }

     char outFI[]="AVE";
     ostringstream o;
     o <<outFI<<".pdb";
     writePdb(o.str().c_str(),size,atomlist, sys);
  
  //check matrix dimension vs. total number of frames   
    cout << "building matrix..."<<endl;  
    int NDIM = size*3;
    if (numFrames < NDIM) {
     cout << "The number of dimensions (" << NDIM 
          << ") is larger than the total number of frames (" 
          << numFrames << ")!\n"
          << "This may lead to poor results.\n" << endl;}
     
    for (int ii=0;ii<NDIM;++ii){
     for (int jj=0;jj<=ii;++jj){
      cov(ii,jj)=(cov(ii,jj)/numFrames)-(avpos[ii]*avpos[jj]);
     }
    }  
   
  //calculate trace of the symmetric matrix
    double tcov=0;
    for (int z=0;z<NDIM;++z){
     tcov+=cov(z,z);
    }
    if (tcov < 0){ throw gromos::Exception(
    "edyn", " trace of the covariance matrix is negative. "
    "Cannot go on. Something might be wrong with your trajectory.\n");
    }
           
  //build up the other triangle of the cov matrix
    for (int ii=0;ii<NDIM;++ii){
     for (int jj=0;jj <= ii;++jj){
      cov(jj,ii)=cov(ii,jj);
     }
    }
    cout << "diagonalizing matrix..." << endl;
      
  //diagonalize the matrix
    std::vector<double> eigen(NDIM);
    for (int g=0;g<NDIM;++g){eigen[g]=0.0;}
     cov.diagonaliseSymmetric(&eigen[0]);

  //calculate trace of the diagonalized matrix
    double  tdcov=0;
    for (int z=0;z<NDIM;++z){	 
     tdcov+=eigen[z];
    }
    if (tdcov < 0){throw gromos::Exception("edyn", " trace of the diagonalized matrix "
                                           "is negative. Cannot go on. Something might "
                                           "be wrong with your trajectory.\n");}


  //compare traces
    if (abs(tdcov-tcov) > (0.01*(tdcov+tcov))) {throw gromos::Exception("edyn", " trace of "
                                           "the covariance and the diagonalized matrix "
                                           "deviate too much. Cannot go on. Something went "
                                           "wrong during the diagonalization. Check your "
                                           "trajectory.\n");}
  //spit out eigenvalues         
    ofstream oeig;
    oeig.open("EIVAL");
    oeig << "Eigenvalues" << endl;
    for (int i=0;i<NDIM;++i){         
     oeig <<  i+1 << " " << eigen[i] << endl;
    }
    oeig.close();
       
  //spit out relative fluctuations
    ofstream orel;
    orel.open("EIFLUC");       
    orel << "Relative Fluctuations of the Eigenvalues" << endl;
    double refl=0;
    for (int i=0;i<NDIM;++i){
     refl+=(eigen[i]/tdcov);
     orel << i+1 << " " << refl << endl;
    }
    orel.close();

  //spit out eigenvectors
    ofstream oeiv;
    oeiv.open("EIVEC");
    oeiv << "Eigenvectors" << endl;  
    for (int ii=0, x=0;ii<=99;++ii){
     for (int jj=0;jj<NDIM;++jj){
      double eivec0 = cov(jj,ii);
      oeiv.setf(ios::right, ios::adjustfield);     
      oeiv << setw(13) << eivec0 << " ";
      x+=1;
      if (x == 6){oeiv << endl;x=0;}
     }
    }
    oeiv.close();
   
  //eigenvector components of the first ten EVs
    for (int i=0;i<10;++i){
     char outFile[]="EVCOMP";
     ostringstream out;
     ofstream outf;
     out <<outFile<<"_"<<i+1;
     outf.open(out.str().c_str());
     for (int j=1;j<=NDIM/3;++j){
      double tmp = 0;
      for (int k=3;k>=1;--k){
       tmp = tmp + cov((3*j-k),i)*cov((3*j-k),i);
      }
      tmp = sqrt(tmp);
      outf <<  j << "  " << tmp << endl;
     }
     outf.close();
    }
    
   //prepare for next loop
 
 int maxsel=20;
 std::vector<double> minprj(maxsel),maxprj(maxsel),avs(maxsel),avsq(maxsel),sig(maxsel); 
 double norm=0.0;
 for (int i=0;i<maxsel;++i){
  minprj[i]=100000.0;
  maxprj[i]=-100000.0;
  avs[i]=0;avsq[i]=0;sig[i]=0; 
}
 int sel[13]={0,1,2,3,4,5,6,7,8,9,19,49,99};

//put selected EV's in separate Matrix
 cout << "putting EV's in EIG" << endl;
 Matrix eig(NDIM,13,0.0);
 for (int i=0;i<13;++i){
  for (int j=0;j<NDIM;++j){
    int vec = sel[i];
    eig(j,i)=cov(j,vec);
  }
 }
 cout << "done" << endl;
 //start loop

 numFrames = 0; 
for(Arguments::const_iterator iter=args.lower_bound("traj");
  iter!=args.upper_bound("traj"); ++iter){
  ic.open(iter->second);
      
// loop over all frames
 while(!ic.eof()){
  numFrames++;
  ic >> sys;	
  (*pbc.*gathmethod)();
  rf.fit(&sys);

//substract average from frame
  for (int i=0, x=0; i<size;++i,x+=2) {
    int j = atomlist[i];
   pvec[i+x]=sys.mol(mol).pos(j)[0]-avpos[i+x];
   pvec[i+1+x]=sys.mol(mol).pos(j)[1]-avpos[i+1+x];
   pvec[i+2+x]=sys.mol(mol).pos(j)[2]-avpos[i+2+x];
  }  

 
  //  Matrix mproj (numFrames,13,0.0);
  for (int i=0;i<13;++i){  
   double proj=0.0;
   for (int j=0;j<NDIM;++j){
    proj+=pvec[j]*eig(j,i);
   }
   proj=-proj;
   int f=sel[i];
   ofstream outfile;
   ostringstream ou;
   char outFile[]="EVPRJ";
   ou <<outFile<<"_"<<f+1;
   outfile.open(ou.str().c_str(), ofstream::app);
   outfile << numFrames << ' ' << proj << endl; 
   outfile.close();
        
   maxprj[i]=((proj)>(maxprj[i]) ? (proj) : (maxprj[i]));
   minprj[i]=((proj)<(minprj[i]) ? (proj) : (minprj[i]));

   //distribution properties
   avs[i]+=proj;
   avsq[i]+=proj*proj;      
  }
 }
}
 double dx=0.0;
  for (int i=0;i<13;++i){
   int f=sel[i];
   ofstream disout;
   ostringstream di;
   char disOut[]="DXPRJ";
   di <<disOut <<"_"<<f+1;
   disout.open(di.str().c_str());
   for (int z=0;z<60;z++){
     dx=(z-4)*((maxprj[i]-minprj[i])/50.0)+minprj[i];
     disout << dx << endl;
   }
   disout.close();
  }
   
//determine averages
    ofstream output;
    output.open("ESSDYN.out"); 
 for (int i=0;i<13;++i){
   int f = sel[i];
  avs[i]=avs[i]/numFrames;
  avsq[i]=avsq[i]/numFrames;
  sig[i]=sqrt((avsq[i]-avs[i]));//*(avsq[i]-avs[i]));
  norm=1/(sqrt(2*3.1415)*(sig[i]));

  output << "Average Projection per frame for EV " << f+1 << ": " << avs[i] << endl;
  output << "Corresponding SD: " << sig[i] << endl; 
 }  


   //write out the extreme structures
   for (int i=0;i<6;++i){
    output << "Maximum and Minimum projections:" << endl;
    output << "Max. Proj. for EV " << i+1 << ": " << maxprj[i] << endl;
    output << "Min. Proj. for EV " << i+1 << ": " << minprj[i] << endl;  
    
    for (int j=0, x=0,y=0,z=0;j<NDIM;++j){
      y=atomlist[z];
     sys.mol(0).pos(y)[x]=avpos[j]+(eig(j,i)*maxprj[i]);
     x+=1; if (x == 3){x=0;z+=1;}
     }
    char outFile[]="PRJMAX";
    ostringstream out;
    out <<outFile<<"_"<<i+1<<".pdb";
    writePdb(out.str().c_str(), size, atomlist, sys);

    for (int j=0, x=0,y=0,z=0;j<NDIM;++j){
      y=atomlist[z];  
     sys.mol(0).pos(y)[x]=avpos[j]+(eig(j,i)*minprj[i]);
            x+=1; if (x == 3){x=0;z+=1;}
     }

    char outF[]="PRJMIN";
    ostringstream ou;
    ou <<outF<<"_"<<i+1<<".pdb";
    writePdb(ou.str().c_str(),size,atomlist, sys);

   }
            
  }

  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void writePdb(const char *name, int dim, vector<int> atoml, System &sys){ 
  ofstream out;
  out.open(name);
  out.setf(ios::fixed, ios::floatfield);
  out.setf(ios::unitbuf);
  out.precision(3);
  int d_count=0;
  int d_resoff=1;
  for (int h=0; h<dim; ++h){
    ++d_count;
    int j=atoml[h];
    int res=sys.mol(0).topology().resNum(j);
    out << "ATOM";
    out.setf(ios::right, ios::adjustfield);
    out << setw(7) << d_count;
    out.setf(ios::left, ios::adjustfield);
    out << "  " <<setw(4) << sys.mol(0).topology().atom(j).name().c_str();
    out << setw(4) << sys.mol(0).topology().resName(res).c_str();
    out.setf(ios::right, ios::adjustfield);
    out << setw(5) << res+d_resoff << "    "
         << setw(8) << sys.mol(0).pos(j)[0]*10
         << setw(8) << sys.mol(0).pos(j)[1]*10
         << setw(8) << sys.mol(0).pos(j)[2]*10
             << "  1.00  0.00" << endl;
  }
  d_resoff+=sys.mol(0).topology().numRes();
  out.close();
}
