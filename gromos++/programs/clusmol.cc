//clustering of molecules 
// this programs reads in the output specified by
// "out" of the programs neighbours


#include "../src/args/Arguments.h"
#include "../src/utils/AtomSpecifier.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gio/InTopology.h"
#include "../src/gmath/Matrix.h"
#include <string>
#include <iomanip>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace utils;
using namespace args;
using namespace gio;
using namespace gcore;
using namespace gmath;

int main(int argc, char **argv){

  char *knowns[] = {"topo","file", "print", "nsm"};
  int nknowns = 4;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@file <output of neighbours>\n";
  usage += "\t@print <all/series>\n";
  usage += "\t@nsm   <number of solute molecules>\n";
  
  
 
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    string sdum;
    int nr_atoms;
    int maxclus=0;
    int frame=0;
    int warning=0;
    
    //read topology
    InTopology it(args["topo"]); 
    System sys(it.system()); 

    //set number of solute atoms
    int nsm=1;
    
    {
      Arguments::const_iterator iter=args.lower_bound("nsm");
      if(iter!=args.upper_bound("nsm")){
        nsm=atoi(iter->second.c_str());
      }
    }
    for(int i=1;i<nsm;i++)
      for(int j=0;j<it.system().numMolecules();j++)
        sys.addMolecule(it.system().mol(j).topology());
        
    int nr_mol=sys.numMolecules();
    
    //open file
    string s=args["file"];
    ifstream fin(s.c_str());
   
    //set printoption
    int print=0;
    {
      Arguments::const_iterator iter=args.lower_bound("print");
      if(iter!=args.upper_bound("print"))
        if(iter->second=="all") print=1;
    }

    //read the first line
    fin >> sdum >> nr_atoms;
    double time=0;

    if(!print) cout << "time\t#clusters\tclustersizes\n";
    
    //loop over the file
    while((fin >> sdum >> time )!=0){
      fin >> sdum >> sdum >> sdum;
      Matrix m(nr_mol, nr_mol, 0);
      if(print) cout << "\nTime: " << time << endl;
      
      for(int i=0;i<nr_atoms;i++){
	int nr_neigh;
        string atom, neigh;
	
	fin >> atom >> nr_neigh;
        AtomSpecifier as(sys,atom);
	for(int j=0; j<nr_neigh;j++){
	  fin >> neigh;
          if(atoi(neigh.c_str())>0){
	    AtomSpecifier bs(sys,neigh);
            m(as.mol(0), bs.mol(0))=1;
	  }
	}
      }
      for(int i=0;i<nr_mol;i++){
        for(int j=0;j<nr_mol;j++){
          if (m(i,j)!=m(j,i)){
            m(i,j)=1;
	    m(j,i)=1;
	    if(print)
              cout << "matrix m fixed for molecules " 
                   << i << " and " << j << endl;
	    else 
	      warning=1;
	  }
	}
      }
      // now work with the matrix to get the clustr
      int cl[nr_mol];
      int curr=-1, max=-1;
      
      for(int i=0;i<nr_mol;i++) cl[i]=-1;
      // for every molecule
      for(int i=0;i<nr_mol;i++){
	//check if it is not yet in a cluster
	if(cl[i]==-1){
	  int found=0;
	  for(int p=0;p<i-1&&found==0;p++){
	    if(m(i,p)==1&&cl[p]!=-1){
	      found=1;
	      curr=cl[p];
	    }
	  }
	  if(found==0){
	    max++;
	    curr=max;
	  }
	  cl[i]=curr;
	}
	else{
	  curr=cl[i];
	}
	for(int j=0;j<nr_mol;j++){
	  if (m(i,j)==1){
	    if(cl[j]==-1) cl[j]=curr;
	    if(cl[j]!=curr){
	      int tmp=cl[j];
	      for(int p=0;p<nr_mol;p++){
		if(cl[p]==tmp)cl[p]=curr;
	      }
	    }
	  }
        }
      }
      //and print out all the clusters
      int count=0;
      int clusno=1;
      vector <int> cluster;
      
      for(int i=0;i<=max;i++){
	int clustersize=0;
	for(int j=0;j<nr_mol;j++){
	  if(cl[j]==i){
	    if(clustersize==0){
	      if(print) cout << "cluster " << clusno++ << ":\n";
	    }
	    clustersize++;
	    if(print) cout << " " << j+1;
	    count++;
	  }
	}
	if(clustersize!=0) {

	  vector<int>::iterator itr=cluster.begin();
	  vector<int>::iterator to =cluster.end();
	  if(itr!=to)
            for(;(clustersize<*itr)&&itr!=to;itr++);
	  cluster.insert(itr,clustersize);
	  
	  if(print) cout << "\nnumber of molecules: "<< clustersize << endl<<endl;
	}
      }
      clusno--;
      if(clusno>maxclus) maxclus=clusno;
      if(print){
	
        cout << "\nsum over all clusters    " << count << " molecules\n";
        cout << "total number of clusters " << clusno << endl;
        cout << "average clustersize      " << double(count)/double(clusno) 
             <<  " molecules\n";
      }
      
      if(!print){
        cout << time << "\t" << cluster.size() << "\t";
        for(unsigned int k=0;k<cluster.size();k++)
          cout << "\t" << cluster[k];
        cout << endl;
      }
      
      frame++;
      
    }
  
    fin.close();
    if(!print&&warning) 
      cout << "\nWarning: Not all molecular neighbour-matrices were symmetric\n"
           << "consider rerunning with @print all\n";
    
    
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
