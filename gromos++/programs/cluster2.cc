// cluster.cc

#include <cassert>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <iostream>

#include "../src/args/Arguments.h"
#include "../src/gio/Ginstream.h"

using namespace std;
using namespace args;

// change this typedef if you have more than 65536 structures
typedef unsigned short the_type;
//typedef unsigned int the_type;

class cluster_parameter
{
public:
  int num;
  int maxstruct;
  int skip;
  int stride;
  double t0;
  double dt;
  double cutoff;
  bool free;
  int precision;
  int number_cluster;
  cluster_parameter() : num(0), maxstruct(-1), skip(0), stride(1), 
			t0(0.0), dt(1.0), cutoff(0.0),
			free(true), precision(10000), number_cluster(0) {};
};

void read_matrix(string const filename, vector< vector < the_type > > &matrix,
		 bool const human, cluster_parameter & cp);

int main(int argc, char **argv){

  char *knowns[] = {"rmsdmat","cutoff", "human", "ref", "maxstruct", "time"};
  int nknowns = 6;

  string usage = argv[0];
  usage += "\n\t@rmsdmat    <rmsd matrix file name>\n";
  usage += "\t@cutoff     <cutoff>\n";
  usage += "\t@time       <t0> <dt>\n";
  usage += "\t[@maxstruct <maximum number of structures to consider>]\n";
  usage += "\t[@human]    (use a human readable matrix)\n";
  usage += "\t[@ref]      (force clustering on the reference structure)\n";
  
  
  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // create the cluster parameters
    cluster_parameter cp;

    // get the cutoff
    cp.cutoff = atof(args["cutoff"].c_str());
    if (cp.cutoff <= 0.0){
      throw gromos::Exception("cluster", "cutoff should be > 0.0");
    }

    // get the time
    {
      Arguments::const_iterator iter=args.lower_bound("time"), 
	to=args.upper_bound("time");
      if(iter!=to){
	cp.t0 = atof(iter->second.c_str());
	++iter;
      }
      if(iter!=to){
	cp.dt = atof(iter->second.c_str());
      }
    }
    double time = cp.t0;
    
    // try to read the maxstruct
    if(args.count("maxstruct") > 0 ) 
      cp.maxstruct = atoi(args["maxstruct"].c_str());
    
    // is the matrix human readable
    bool human=false;
    if(args.count("human") >= 0) human=true;
    
    // is the clustering free
    if(args.count("ref") >=0) cp.free=false;
    
    // create the data structure
    vector< vector <the_type> > pairs;
    
    // read matrix
    read_matrix(args["rmsdmat"], pairs, human, cp);
    
    // now we are almost done
    size_t num=pairs.size();
    vector<int> taken(num, -1);
    vector<the_type> central_member;
    vector<vector <the_type> > cluster;
    int remaining=num;
    int clustercount=0;

    // set the first cluster
    central_member.push_back(0);
    cluster.push_back(pairs[0]);
    if(cp.free)
      pairs[0].clear();
    taken[0]=0;

    if(!cp.free){
      // mark them as taken
      for(size_t j=0, jto=pairs[0].size(); j != jto; ++j){
	taken[pairs[0][j]]=0;
	pairs[pairs[0][j]].clear();
      }
      pairs[0].clear();
     
      // take them out 
      remaining-=pairs[0].size();
      for(size_t i=1; i < num; ++i){
	vector< the_type > temp;
	for(size_t j = 0, jto=pairs[i].size(); j!=jto; ++j){
	  if(taken[pairs[i][j]]==-1){
	    temp.push_back(pairs[i][j]);
	  }
	}
	pairs[i]=temp;
      }
    }
    
    
    //while(remaining){
    while(true){
      
      // search for the one with the largest number of neighbours
      size_t maxsize=0;
      size_t maxindex=0;
      for(size_t i=0; i < num; ++i){
	if(pairs[i].size() > maxsize) {
	  maxsize = pairs[i].size();
	  maxindex = i;
	}
      }
      if(!maxsize) break;
      
      // put them in
      clustercount++;
      central_member.push_back(maxindex);
      cluster.push_back(pairs[maxindex]);
      
      // and take them out
      
      remaining-=pairs[maxindex].size();

      for(size_t j=0, jto=pairs[maxindex].size(); j != jto; ++j){
	taken[pairs[maxindex][j]]=clustercount;
	pairs[pairs[maxindex][j]].clear();
      }
      pairs[maxindex].clear();
      for(size_t i=0; i < num; ++i){
	vector< the_type > temp;
	for(size_t j = 0, jto=pairs[i].size(); j!=jto; ++j){
	  if(taken[pairs[i][j]] == -1){
	    temp.push_back(pairs[i][j]);
	  }
	}
	pairs[i]=temp;
      }
    } // while remaining
    
    // NOW PRODUCE OUTPUT
    // time series
    {
      ofstream fout("cluster_ts.dat");
      for(size_t i =1; i < num; ++i){
	fout << setw(10) << time
	     << setw(10) << i
	     << setw(10) << taken[i] << endl;
	time += cp.dt;
      }
    }
    // cluster contents
    cp.number_cluster = cluster.size();
    if(cp.free) cp.number_cluster--;
    
    {
      ofstream fout("cluster_structures.dat");
      fout << "TITLE\n"
	   << "\tClustering from " << args["rmsdmat"] << "\n"
	   << "\tCutoff : " << cp.cutoff << "\n"
	   << "\tTotal number of structures : " << num << "\n";
      if(cp.free) fout << "\tFree clustering performed\n";
      else fout << "\tFirst cluster forced to contain reference structure\n";
      fout << "\tTotal number of clusters found : " << cp.number_cluster
	   << "\nEND\n";
      fout << "CLUSTERPARAMETERS\n"
	   << "#  structs   maxstruct        skip      stride\n"
	   << setw(10) << cp.num
	   << setw(12) << cp.maxstruct
	   << setw(12) << cp.skip
	   << setw(12) << cp.stride << "\n"
	   << "#       t0          dt\n"
	   << setw(10) << cp.t0
	   << setw(12) << cp.dt << "\n"
	   << "#   cutoff        free   precision\n"
	   << setw(10) << cp.cutoff
	   << setw(12) << ((cp.free)?1:0)
	   << setw(12) << rint(log(double(cp.precision)) / log(10.0)) << "\n"
	   << "# clusters\n"
	   << setw(10) << cp.number_cluster << "\n"
      	   << "END\n";
      
      fout << "CLUSTER\n";
      fout << "#    clu  center    time    size\n";
      for(size_t i=0, ito=cluster.size(); i!=ito; ++i){
	fout << setw(8) << i
	     << setw(8) << central_member[i];
	if(central_member[i]==0)
	  fout << setw(8) << "ref";
	else
	  fout << setw(8) << cp.t0 + (central_member[i] -1)* cp.dt;
	fout << setw(8) << cluster[i].size() << endl;
      }
      fout << "END\n";
      fout << "MEMBERS\n";
      for(size_t i=0, ito=cluster.size(); i!=ito; ++i){
	fout << setw(6) << i;
	for(size_t j=0, jto=cluster[i].size(); j!=jto; ++j){
	  if(j%10==0 && j!=0) fout << "\n      "; 
	  fout << setw(6) << cluster[i][j];
	}
	if((cluster[i].size())%10!=0) fout << "\n";
      }
      if((cluster[cluster.size()-1].size())%10==0) fout << "\n";
      fout << "END\n";
    }
    
    // standard output
    {
      cout << "# Clustering from " << args["rmsdmat"] << "\n"
	   << "# Cutoff: " << cp.cutoff << "\n"
	   << "# Total number of structures: " << num << "\n";
      if(cp.free) cout << "# Free clustering performed\n";
      else cout << "# First cluster forced to contain reference structure\n";
      cout << "# Total number of clusters found: " << cp.number_cluster << "\n";
      
      cout << "#\n";
      cout << "#    clu    size\n";
      if(cp.free) cout << setw(8) << "ref";
      else cout << setw(8) << 0;
      cout << setw(8) << cluster[0].size() << endl;
      
      for(size_t i=1, ito=cluster.size(); i!=ito; ++i){
	cout << setw(8) << i
	     << setw(8) << cluster[i].size() << endl;
      }
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}



void read_matrix(string const filename, vector< vector < the_type > > &matrix,
		 bool const human, cluster_parameter & cp)
{
  if(human){
    gio::Ginstream gin(filename);
    if(!gin.stream()){
      throw gromos::Exception("cluster", "Error opening rmsdmat file\n");
    }
    
    string sdum;
    
    gin.getline(sdum);
    if(sdum!="RMSDMAT")
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "expected RMSDMAT block, got "+sdum);
    gin.getline(sdum);
    istringstream is(sdum);
    if(!(is >> cp.num >> cp.skip >> cp.stride))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read number of structures\n");
    if(cp.num<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative number structures\n"); 
    if(cp.num > pow(2.0, 8.0 * sizeof(the_type)))
      throw gromos::Exception("cluster", "GROMOS96 ERROR: number of "
			      "structures is too large for data type. "
			      "Change typedef in program and recompile");
      
    if(cp.skip<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative number of skipped structures\n");
    if(cp.stride<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative number of structure stride\n");
    gin.getline(sdum);
    is.clear();
    is.str(sdum);
    if(!(is >> cp.precision))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read precision\n");
    if(cp.precision<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative precision\n");

    double cuttest=cp.cutoff * cp.precision / 10;
    if(fabs(double(int(cuttest))-cuttest) > cp.cutoff / 100.0)  
      throw gromos::Exception("cluster", "A cutoff with this precision "
			      "requires a higher precision in the rmsd "
			      "matrix. \nYes that means that you have to "
			      "redo your matrix $%^$#$@$%!!");
    
    if(cp.maxstruct < 0) cp.maxstruct = cp.num;
    else cp.maxstruct++;
    
    int icutoff=int(cp.cutoff * cp.precision);
    
    matrix.resize(cp.maxstruct);
    
    int ii,jj, rmsd;
    
    for(int i=0; i < cp.num; ++i){
      if(i < cp.maxstruct){
	//cout << "diagonal " << i << endl;
	
	matrix[i].push_back(i);
      }
      
      for(int j=i+1; j < cp.num; ++j){
	gin.getline(sdum);
	is.clear();
	is.str(sdum);
	if(!(is >> ii >> jj >> rmsd))
	  throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
				  "could not read line:" + is.str());
	if(ii!=i || jj!=j){
	  ostringstream os;
	  os << "Now we are really confused!\n"
	     << i << " = " << ii << " or " << j << " = " << jj << "\n";
	  throw gromos::Exception("cluster", os.str());
	}
	
	if(ii < cp.maxstruct && jj < cp.maxstruct && rmsd < icutoff){
	  //cout << "ij " << ii << " " << jj << endl;
	  
	  matrix[ii].push_back(jj);
	  if(ii>0)
	    matrix[jj].push_back(ii);
	}
      }
    }
    if(!(gin.getline(sdum)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "where is the END?");
    if(sdum.substr(0,3)!="END")
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "Trailing data on file. Expected 'END', got"
			      + sdum);
  }
  else{
    ifstream fin(filename.c_str(), ios::in | ios::binary);
    if(!fin){
      throw gromos::Exception("cluster", "Error opening rmsdmat file\n");
    }
    
    if(!fin.read((char*)&cp.num, sizeof(int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read number of structures\n");
    if(cp.num<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative number of structures\n");
    if(cp.num > pow(2.0, 8.0 * sizeof(the_type)))
      throw gromos::Exception("cluster", "GROMOS96 ERROR: number of "
			      "structures is too large for data type. "
			      "Change typedef in program and recompile");
    if(!fin.read((char*)&cp.skip, sizeof(int)))
	throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read skip\n");
    if(cp.skip<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read number of skipped structures\n");
    if(!fin.read((char*)&cp.stride, sizeof(int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read stride\n");
    if(cp.stride<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative number of stride structures\n");

    if(!fin.read((char*)&cp.precision, sizeof(int)))
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "could not read precision\n");
    if(cp.precision<0)
      throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
			      "read negative precision\n");

    double cuttest=cp.cutoff * cp.precision / 10;
    //cout << "cuttest " << cuttest << " " << cp.cutoff << " " << cp.precision << endl;
    if(fabs(double(int(cuttest))-cuttest) > cp.cutoff / 100.0)
      throw gromos::Exception("cluster", "A cutoff with this precision "
			      "requires a higher precision in the rmsd "
			      "matrix. \nYes that means that you have to "
			      "redo your matrix $%^$#$@$%!!");
    if(cp.maxstruct < 0) cp.maxstruct = cp.num;
    else cp.maxstruct++;
    
    unsigned int icutoff=unsigned(cp.cutoff * cp.precision);
    
    matrix.resize(cp.maxstruct);
    

    if(cp.precision < 1e5){
      typedef unsigned short ushort;
      
      ushort rmsd;
      
      for(int i=0; i < cp.num; ++i){
	if(i<cp.maxstruct)
	  matrix[i].push_back(i);
	
	for(int j=i+1; j < cp.num; ++j){
	  if(!fin.read((char*)&rmsd,sizeof(ushort)))
	    throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
				    "file corrupt");
	  if(i < cp.maxstruct && j < cp.maxstruct && rmsd < short(icutoff)){
	    matrix[i].push_back(j);
	    if(i>0)
	      matrix[j].push_back(i);
	  }
	}
      }
    
    }
    else{
      unsigned rmsd;
      
      for(int i=0; i < cp.num; ++i){
	if(i<cp.maxstruct)
	  matrix[i].push_back(i);
	
	for(int j=i+1; j < cp.num; ++j){
	  if(!fin.read((char*)&rmsd,sizeof(unsigned)))
	    throw gromos::Exception("cluster", "Error while reading rmsdmat file\n"
				    "file corrupt");
	  if(i < cp.maxstruct && j < cp.maxstruct && rmsd < icutoff){
	    matrix[i].push_back(j);
	    if(i>0)
	      matrix[j].push_back(i);
	  }
	}
      }
    }
  }
}
