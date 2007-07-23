// gio_InTopology.cc

#include <cassert>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "Ginstream.h"
#include "../gcore/System.h"
#include "../gcore/PtTopology.h"
#include "InPtTopology.h"

using namespace std;
using namespace gcore;
using gio::InPtTopology_i;
using gio::InPtTopology;

// Implementation class
class gio::InPtTopology_i: public gio::Ginstream
{
  friend class InPtTopology;
  gcore::PtTopology d_pttopo;
  std::map<std::string, std::vector<std::string> > d_blocks;
  /**
   * the init function reads in the whole file into the map of blocks and 
   * reads in the topology version
   */
  void init();
  /**
   * parseForceField takes all blocks that end up in the forcefield
   * and stores the information in d_pttopo
   */
  void parsePtTopology(int start);
  /**
   * _initBlock is a function that initialized the reading of
   * a block. It checks whether the block is read in, and returns 
   * the number of elements that are to be read in from it.
   */
  int _initBlock(std::vector<std::string> &buffer,
		 std::vector<std::string>::const_iterator &it,
		 const std::string blockname);

  InPtTopology_i(std::string &s):d_blocks()
  {
    this->open(s);
  }
  
};

gio::InPtTopology::InPtTopology(std::string name, int start){
  d_this = new InPtTopology_i(name);
  d_this->init();
  d_this->parsePtTopology(start);
}

gio::InPtTopology::~InPtTopology()
{
  delete d_this;
}

const string InPtTopology::title()const
{
  return d_this->title();
}

const gcore::PtTopology &InPtTopology::ptTopo()const{
  return d_this->d_pttopo;
}

int gio::InPtTopology_i::_initBlock(std::vector<std::string> &buffer,
	       std::vector<std::string>::const_iterator &it,
	       const string blockname)
{
  int n;
  
  buffer.clear();
  buffer=d_blocks[blockname];
  if(buffer.size() < 3)
    throw InPtTopology::Exception("Topology file"+name()+
				" is corrupterd. No (or empty) "+blockname+
				" block!");
  if(buffer[buffer.size()-1].find("END")!=0)
    throw InPtTopology::Exception("Topology file " + name() +
				     " is corrupted. No END in "+blockname+
				     " block. Got\n"
				     + buffer[buffer.size()-1]);

  it=buffer.begin() +1;
  _lineStream.clear();
  _lineStream.str(*it);
  _lineStream >> n;
  ++it;
  return n;
}

void gio::InPtTopology_i::init(){
  
  if(!stream())
    throw InPtTopology::Exception("Could not open topology file."+name());
  
  // First read the whole file into the map
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  while(!stream().eof()){
    getblock(buffer);
    if(buffer.size()){
      d_blocks[buffer[0]]=buffer;
      buffer.clear();
    }
  }
}

void gio::InPtTopology_i::parsePtTopology(int start)
{
  std::vector<std::string> buffer;
  std::vector<std::string>::const_iterator it;
  
  std::string nm;
  int k,l, n, iiac, numat, numpt;
  double mass, dq;
  
  if(d_blocks.count("PERTATOM")){
    int numat = _initBlock(buffer, it, "PERTATOM");
    
    d_pttopo.setSize(numat,1);
    for (n=0; it < buffer.end()-1; ++it, ++n){
      _lineStream.clear();
      _lineStream.str(*it);
      _lineStream >> k;
      if(n==0&&start>=0) start-=k;
      d_pttopo.setAtomNum(n,k+start);
      _lineStream >> l >> nm >> iiac >> mass >> dq >> l >> l;
      d_pttopo.setAtomName(n, nm);
      d_pttopo.setIac(n,0,iiac-1);
      d_pttopo.setCharge(n,0,dq);
    }
    if(n!=numat)
      throw InPtTopology::Exception("Perturbation topology file "+name()+
				    " is corrupted. Failed to read all atoms");
    
  }
  else if(d_blocks.count("MPERTATOM")){
    buffer.clear();
    buffer=d_blocks["MPERTATOM"];
    
    it=buffer.begin()+1;
    
    if(buffer.size() < 3 ) 
      throw InPtTopology::Exception("Topology file "+name()+
                          " is corrupted. Empty MPERTATOM block!");
    if(buffer[buffer.size()-1].find("END")!=0)
      throw InPtTopology::Exception("Topology file " + name() +
				  " is corrupted. No END in MPERTATOM"
				  " block. Got\n"
				  + buffer[buffer.size()-1]);
     _lineStream.clear();
     _lineStream.str(*it);
     _lineStream >> numat >> numpt;
     d_pttopo.setSize(numat, numpt);
     ++it;
     _lineStream.clear();
     _lineStream.str(*it);
     for(int j=0; j<d_pttopo.numPt(); ++j){
       _lineStream >> nm;
       d_pttopo.setPertName(j,nm);
     }
     ++it;
     
     for (n=0; it < buffer.end()-1; ++it, ++n){
       _lineStream.clear();
       _lineStream.str(*it);
       _lineStream >> k;
       if(n==0&&start>=0) start-=k;
       d_pttopo.setAtomNum(n,k+start);
       _lineStream >> nm;
       d_pttopo.setAtomName(n, nm);
       for(int j=0; j<d_pttopo.numPt(); ++j){
	 _lineStream >> iiac >> dq;
	 d_pttopo.setIac(n,j,iiac-1);
	 d_pttopo.setCharge(n,j,dq);
       }
     }
     if(n!=numat)
       throw InPtTopology::Exception("Perturbation topology file "+name()+
				     " is corrupted. Failed to read all atoms");
  }
  else
    throw InPtTopology::Exception("No PERTATOM or MPERTATOM block in "
				  "perturbation topology file\n");
  
}





