// cluster_analyze.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/utils/Rmsd.h"
#include "../src/fit/Reference.h"
#include "../src/fit/RotationalFit.h"
#include "../src/gio/InG96.h"
#include "../src/gio/OutG96.h"
#include "../src/gio/OutPdb.h"
#include "../src/gcore/Box.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/Boundary.h"
#include "../src/gio/OutPdb.h"
#include "../src/gmath/Vec.h"

#include <vector>
#include <set>
#include <fstream>
#include <strstream>
#include <iomanip>
#include <algorithm>
#include <iostream>

using namespace gcore;
using namespace gio;
using namespace utils;
using namespace bound;
using namespace args;
using namespace fit;


class PDB_plot_data {
  public:
    string pdb_file_name;
    int  index_in_ctraj;
};
bool comes_before(const PDB_plot_data &x1,const PDB_plot_data &x2) {
   return ((x1.index_in_ctraj) < (x2.index_in_ctraj));
}

int main(int argc, char **argv){

  char *knowns[] = {"topo","pbc","mat","ctraj","cut","nclus","plot"};
  int nknowns = 7;

  string usage = argv[0];
  usage += "\n\t-topo <topology>\n";
  usage += "\t-pbc <boundary type>\n";
  usage += "\t-ctraj <clustering trajectory file>\n";
  usage += "\t-mat <clustering matrix file>\n";
  usage += "\t-cut <cutoff>\n";
  usage += "\t-nclus <number of clusters wanted>\n";
  usage += "\t-plot <number of pdb structure files to write per cluster>\n";

  // prepare output
  cout.setf(ios::right, ios::adjustfield);
  cout.setf(ios::fixed, ios::floatfield);

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // get cutoff
    float cutoff=0.0;
    args.check("cut",1);
    {
      Arguments::const_iterator iter=args.lower_bound("cut");
      if(iter!=args.upper_bound("cut")){
	cutoff=atof(iter->second.c_str());
      }
    }
    if (cutoff <= 0.0) {
      throw gromos::Exception("cluster_analyze", "Invalid cutoff" );
    }

    // get number of clusters wanted
    int want_clusters=0;
    args.check("nclus",1);
    {
      Arguments::const_iterator iter=args.lower_bound("nclus");
      if(iter!=args.upper_bound("nclus")){
	want_clusters=atoi(iter->second.c_str());
      }
    }
    if (want_clusters < 1) {
      throw gromos::Exception("cluster_analyze", "Invalid number of clusters" );
    }

    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());
    
    // read coordinates...
    bool  usetraj=false;
    InG96 ic;

    try{
      args.check("ctraj",1);
      ic.open(args["ctraj"]);
      usetraj=true;
    }
    catch(const Arguments::Exception &){
      usetraj=false;
    }
    int frames_in_traj=0;
    if (usetraj) {
      while (!ic.eof()) {
        ic >> sys; 
        frames_in_traj++;
      }
      ic.close();
    }

    // do we want plotting of structures in clusters ?
    bool makeplot=false;
    int want_plots=0;
    try{
      args.check("plot",1);
      makeplot=true;
      if (!usetraj) {
        throw gromos::Exception("cluster_analyze", "Plotting needs trajectory" );
      }
    }
    catch(const Arguments::Exception &){
      makeplot=false;
    }
    if (makeplot) {
      Arguments::const_iterator iter=args.lower_bound("plot");
      if(iter!=args.upper_bound("plot")){
	want_plots=atoi(iter->second.c_str());
      }
    }
    if (makeplot && (want_plots < 1)) {
      throw gromos::Exception("cluster_analyze", "Invalid number of plots" );
    }

    // what names should plots have ?
    string pdb_name_prefix="./clusterpdb";

    // open matrix file
    ifstream mf;
    string mfname_str="";
    char *mfname;
    try{
      args.check("mat",1);
      mfname_str=""+args["mat"];
      int cstrlen=mfname_str.length()+1;
      mfname=new char[cstrlen];
      copy(mfname_str.c_str(),mfname_str.c_str()+cstrlen,mfname);
      mf.open(mfname);
    }
    catch(const Arguments::Exception &){
      mfname_str="./cluster_matrix_file.dat";
      int cstrlen=mfname_str.length()+1;
      mfname=new char[cstrlen];
      copy(mfname_str.c_str(),mfname_str.c_str()+cstrlen,mfname);
      mf.open(mfname);
    }
    if (!mf) {
      char message_buffer[120] = "";
      ostrstream ErMs(message_buffer,100);
      ErMs << "  Error opening matrix file " << mfname_str << endl;
      throw gromos::Exception("cluster_analyze", message_buffer );
    }
    
    // count through matrix file
    int matsiz=0;
    int nframes=0;
    {
      int i,j;
      float f;
      while (mf >> i >> j >> f) {
        matsiz++;
        if ((i+1) > nframes) nframes=i+1;
        if ((j+1) > nframes) nframes=j+1;
      }
    }
    int propersize = (nframes*nframes - nframes)/2;
    if (matsiz != propersize) {
       throw gromos::Exception("cluster_analyze", "Wrong size matrix file" );
    }
    if (usetraj) {
      if (frames_in_traj != nframes) {
        throw gromos::Exception("cluster_analyze", "Wrong size trajectory file" );
      }
    }
    mf.close();
    mf.clear();

    // read in the rmsd matrix
    float *rmsdmat;
    rmsdmat = new float[matsiz];
    mf.open(mfname);
    if (!mf) {
      char message_buffer[120] = "";
      ostrstream ErMs(message_buffer,100);
      ErMs << "  Error opening matrix file " << mfname_str << endl;
      throw gromos::Exception("cluster_analyze", message_buffer );
    }
    {
      int i,j,k,l,ix;
      float f;
      while (mf >> i >> j >> f) {
        if (i < j) {
          k=i; l=j;
        }
        else if (j < i) {
          k=j; l=i;
        }
        else {
          throw gromos::Exception("cluster_analyze", "Inconsistent matrix file" );
        }
        ix=k*nframes-((k*k+k)/2)+l-k-1;
        rmsdmat[ix]=f;
      }
    }
    mf.close();
    delete [] mfname;


 
    // clustering algorithm
    bool *taken;
    taken  = new bool[nframes];
    for (int ii=0; ii<nframes; ii++) taken[ii]=false;

    int  *neighb;
    neighb = new int[nframes];

    vector<int> *the_clusters;
    the_clusters = new vector<int>[want_clusters];

    int cluster_no = 0;
    while (cluster_no < want_clusters) {
      //reset number of neighbors for each clustering turn
      for (int ii=0; ii<nframes; ii++) neighb[ii]=0;

      int ix;

      // count neighbors closer than cutoff among untaken
      for (int ii=0; ii < nframes-1; ii++) {
        for (int jj=ii+1; jj < nframes; jj++) {
          if ((!taken[ii]) && (!taken[jj])) {
            ix=ii*nframes-((ii*ii+ii)/2)+jj-ii-1;
            if (rmsdmat[ix] < cutoff) {
              neighb[ii]++;
              neighb[jj]++;
            }
          }
        }
      }

      // find the one with most neighbors
	      int hasmost=-1;
	      int maxneigh=-1;
      for (int ii=0; ii < nframes; ii++) {
        if (!taken[ii]) {
          if (neighb[ii] > maxneigh) {
            hasmost=ii;
            maxneigh=neighb[ii];
          }
        }
      }
      if (hasmost < 0) {
        throw gromos::Exception("cluster_analyze",
                                "Given cutoff yields too few clusters" );
      }
      the_clusters[cluster_no].push_back(hasmost);
      taken[hasmost]=true;

      // find its neighbors, make a cluster      
      for (int ii=0; ii < nframes; ii++) {
        if ((ii != hasmost) && !taken[ii]) {
          int k,l;
          if (ii < hasmost) {
            k=ii; l=hasmost;
          }
          else {
            k=hasmost; l=ii;
          }
          ix=k*nframes-((k*k+k)/2)+l-k-1;
          if (rmsdmat[ix] < cutoff) {
            the_clusters[cluster_no].push_back(ii);
            taken[ii]=true;
          }
        }
      }

      cluster_no++;
    }



    // write results
    cout << endl;
    cout << "Clustering results" << endl;
    cout << endl;
    cout << endl;
    for (int ii=0; ii < want_clusters; ii++) {
      cout << "cluster no " << setw(7) << ii
           << " has " << setw(7) << the_clusters[ii].size()
           << " members" << endl << endl;
      cout << "center is  ";
      const int maxontheline = 10;
      int ontheline = maxontheline - 1;
      for (vector<int>::size_type jj=0; jj != the_clusters[ii].size() ; ++jj) {
        cout << setw(7) << the_clusters[ii][jj];
        ontheline++;
        if (ontheline >= maxontheline) {
          ontheline = 0;
          cout << endl;
        } 
      }
      cout << endl << endl;
    }

    // write PDB files of structures, if requested
    if (makeplot) {
      char filename_buffer[400] = "";
      PDB_plot_data outelement;
      vector<PDB_plot_data> outlist;
      for (int ii=0; ii < want_clusters; ii++) {
         int numplots=want_plots;
         set<int> taken,notyet;
         for (vector<int>::const_iterator xx=the_clusters[ii].begin();
              xx != the_clusters[ii].end(); xx++) {
           notyet.insert(*xx);
         }
         if (static_cast<int>(the_clusters[ii].size()) < numplots) {
           numplots = the_clusters[ii].size();
         }
         {
           // always take the center for plotting first
           ostrstream PDBFileName(filename_buffer,350);
           PDBFileName << pdb_name_prefix;
           PDBFileName << "_c" << ii;
           PDBFileName << "_n0";
           PDBFileName << "_i" << the_clusters[ii][0];
           PDBFileName << ".pdb";
           PDBFileName << '\0';
           outelement.pdb_file_name=filename_buffer;
           outelement.index_in_ctraj=the_clusters[ii][0];
           outlist.push_back(outelement);
           taken.insert(the_clusters[ii][0]);
           notyet.erase(the_clusters[ii][0]);
         }
         for (int jj=1;jj < numplots; jj++) {
           // find most distant from the ones taken
           double maxdist=-1.0;
           int maxdist_nr=-1;
           for (set<int>::const_iterator xx=notyet.begin(); xx!=notyet.end(); xx++) {
             double dist=0;
             for (set<int>::const_iterator yy=taken.begin(); yy!=taken.end(); yy++) {
               int k,l,ix;
               if (*xx < *yy) {
                 k=*xx; l=*yy;
               }
               else {
                 k=*yy; l=*xx;
               }
               ix=k*nframes-((k*k+k)/2)+l-k-1;
               dist += rmsdmat[ix]*rmsdmat[ix];
             }
             if (dist > maxdist) {
               maxdist=dist;
               maxdist_nr= *xx;
             }
           }
           ostrstream PDBFileName(filename_buffer,350);
           PDBFileName << pdb_name_prefix;
           PDBFileName << "_c" << ii;
           PDBFileName << "_n" << jj;
           PDBFileName << "_i" << maxdist_nr;
           PDBFileName << ".pdb";
           PDBFileName << '\0';
           outelement.pdb_file_name=filename_buffer;
           outelement.index_in_ctraj=maxdist_nr;
           outlist.push_back(outelement);
           taken.insert(maxdist_nr);
           notyet.erase(maxdist_nr);
         }
         cout << ii << endl;
         for (set<int>::const_iterator xx=taken.begin(); xx !=taken.end(); xx++) {
           cout << *xx << " ";
         }
         cout << endl;
         for (set<int>::const_iterator xx=notyet.begin(); xx !=notyet.end(); xx++) {
           cout << *xx << " ";
         }
         cout << endl;
      }
      sort(outlist.begin(),outlist.end(),comes_before);

      ofstream pf;
      ic.open(args["ctraj"]);
      {
        int current_frame=0;
        vector<PDB_plot_data>::const_iterator iter=outlist.begin();
        while ((!ic.eof()) && (iter!=outlist.end()) ) {
          ic >> sys;
          if ( (iter->index_in_ctraj) == current_frame) {
            pf.open(iter->pdb_file_name.c_str());
            if (!pf) {
              throw gromos::Exception("cluster_analyze",
                                  "Error writing pdb file");
            }
            OutPdb PDBOutFile(pf);
            PDBOutFile.writeTitle("Output from cluster_analyze program");
            PDBOutFile << sys;
            pf.close();
            iter++;
          }
          current_frame++;
        }
      }
      ic.close();
    }


  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

