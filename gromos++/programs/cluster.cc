// cluster.cc

#include "../src/args/Arguments.h"
#include "../src/utils/RmsdMat.h"
#include "../src/utils/Cluster.h"

#include <vector>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace utils;
using namespace args;

int main(int argc, char **argv){

  char *knowns[] = {"mat","cut"};
  int nknowns = 2;

  string usage = argv[0];
  usage += "\n\t@mat <rmsd matrix file name>\n";
  usage += "\t@cut <cutoff>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    Arguments::const_iterator iter;
    string errmsg;

    // get cutoff
    float cutoff = 0.0;
    args.check("cut", 1);
    iter = args.lower_bound("cut");
    if(iter != args.upper_bound("cut"))
	cutoff = atof(iter->second.c_str());
    if (cutoff <= 0.0){
      errmsg = usage;
      errmsg += "Invalid cutoff: ";
      errmsg += "Must be larger than zero.\n";
      throw Arguments::Exception(errmsg);
    }

    // open matrix file and count the number of frames
    args.check("mat",1);
    ifstream mf;
    mf.open(args["mat"].c_str());
    if (!mf) {
      string errmsg = "Error opening matrix file: ";
      errmsg += args["mat"];
      errmsg += "\n";
      throw Arguments::Exception(errmsg);
    }
    unsigned int nframes = 0;
    unsigned int i = 0;
    unsigned int j = 0;
    float f = 0;
    while (!i) {
      mf >> i >> j >> f;
      nframes++;
    }
    mf.close();
    mf.clear();

    // read in the rmsd matrix
    RmsdMat rmsdmat(nframes);
    mf.open(args["mat"].c_str());
    if (!mf) {
      string errmsg = "Error reopening matrix file: ";
      errmsg += args["mat"];
      errmsg += "\n";
      throw gromos::Exception("cluster_analyze", errmsg);
    }
    while (mf >> i >> j >> f)
      rmsdmat.insert(i, j, f);
    mf.close();

    // clustering algorithm
    vector<Cluster> clusters(nframes, Cluster());
    for(i = 0; i < nframes; i++) clusters[i].center = i;

    unsigned int maxNumNeighbors = 0;
    unsigned int hasMostNeighbors = 0;
    bool firstLoop = true;
    while(maxNumNeighbors || firstLoop){
      firstLoop = false;
      maxNumNeighbors = 0;
      hasMostNeighbors = 0;
      // reset number of neighbors
      for(i = 0; i < clusters.size(); i++) 
        clusters[i].neighbors.clear();

      // count neighbors closer than cutoff among untaken
      for(i = 0; i < clusters.size() - 1; i++){
        for(j = i + 1; j < clusters.size(); j++){
          if(!clusters[i].is_taken && 
            !clusters[j].is_taken && 
            rmsdmat.retrieve(i, j) < cutoff){
            clusters[i].neighbors.push_back(j);
            clusters[j].neighbors.push_back(i);
          }
        }
      }

      // find the one with most neighbors
      for(i = 0; i < nframes; i++){
        if(!clusters[i].is_taken && 
          clusters[i].neighbors.size() > maxNumNeighbors){
          hasMostNeighbors = i;
          maxNumNeighbors = clusters[i].neighbors.size();
        }
      }
      if(!maxNumNeighbors) break;

      #define CLUSTER clusters[hasMostNeighbors]

      CLUSTER.is_taken = true;
      for(i = 0; i != CLUSTER.neighbors.size(); ++i)
        clusters[CLUSTER.neighbors[i]].is_taken = true;

      /*
       * When e.g. we want to look at the clusters, we want
       * cluster member i to be the most distant from all the
       * previous members (such that the visualization represents
       * the cluster well). So we have to sort the neighbors 
       * according to this criterion.
      */
      vector<int> sortedNeighbors;
      double rmsdSum;
      double maxRmsdSum;
      int mostDistNeighbor;
      while(CLUSTER.neighbors.size()){
        rmsdSum = 0;
        maxRmsdSum = 0;
        mostDistNeighbor = 0;
        for(i = 0; i < CLUSTER.neighbors.size(); i++){
          rmsdSum = rmsdmat.retrieve(CLUSTER.center, CLUSTER.neighbors[i]);
          for(j = 0; j < sortedNeighbors.size(); j++)
            rmsdSum += rmsdmat.retrieve(sortedNeighbors[j], CLUSTER.neighbors[i]);
          if(rmsdSum > maxRmsdSum){
            maxRmsdSum = rmsdSum;
            mostDistNeighbor = i;
          }
        }
        sortedNeighbors.push_back(CLUSTER.neighbors[mostDistNeighbor]);
        CLUSTER.neighbors.erase(CLUSTER.neighbors.begin() + mostDistNeighbor);
      }
      CLUSTER.neighbors = sortedNeighbors;
      



      cout << "MEMBERS" << endl;
      cout << setw(7) << CLUSTER.neighbors.size() + 1 << endl;
      cout << "END" << endl;
      cout << "CENTER" << endl;
      cout << setw(7) << CLUSTER.center + 1 << endl;
      cout << "END" << endl;
      cout << "NEIGHBORS" << endl;
      for(j = 0; j != CLUSTER.neighbors.size(); ++j)
        cout << setw(7) << CLUSTER.neighbors[j] + 1 << endl;
      cout << "END" << endl;
      #undef CLUSTER

    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}
