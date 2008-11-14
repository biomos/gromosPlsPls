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
  int k,l, n, iiac[2], numat, numpt;
  double mass[2], dq[2], alphaLJ, alphaCRF;
  
  if(d_blocks.count("PERTATOMPARAM")) {
    int numat = _initBlock(buffer, it, "PERTATOMPARAM");

    d_pttopo.setSize(numat, 1);
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> k;
      if (n == 0 && start >= 0) start -= k;
      d_pttopo.setAtomNum(n, k + start);

      _lineStream >> l >> nm >> iiac[0] >> mass[0] >> dq[0]
                  >> iiac[0] >> mass[0] >> dq[0] >> alphaLJ >> alphaCRF;

      if (_lineStream.fail()) {
        ostringstream os;
        os << "Bad line in PERTATOMPARAM block (line " << n + 1 << ").";
        throw InPtTopology::Exception(os.str());
      }

      d_pttopo.setAtomName(n, nm);
      d_pttopo.setIac(n, 0, iiac[0] - 1);
      d_pttopo.setMass(n, 0, mass[0]);
      d_pttopo.setCharge(n, 0, dq[0]);
      d_pttopo.setIac(n, 1, iiac[1] - 1);
      d_pttopo.setMass(n, 1, mass[1]);
      d_pttopo.setCharge(n, 1, dq[1]);
      d_pttopo.setAlphaLJ(n, alphaLJ);
      d_pttopo.setAlphaCRF(n, alphaCRF);
    }
    if (n != numat)
      throw InPtTopology::Exception("Perturbation topology file " + name() +
            " is corrupted. Failed to read all atoms");

    if (d_blocks.count("PERTPOLPARAM")) {
          int numat = _initBlock(buffer, it, "PERTPOLPARAM");

    if (numat != d_pttopo.numAtoms()) {
      ostringstream os;
      os << "Error in PERTPOLPARAM block: The block contains " << numat
         << " atoms but the PERTATOMPARAM block " << d_pttopo.numAtoms()
         << " atoms.";
      throw InPtTopology::Exception(os.str());
    }
    for (n = 0; it < buffer.end() - 1; ++it, ++n) {
      _lineStream.clear();
      _lineStream.str(*it);

      _lineStream >> k;
      if (n == 0 && start >= 0) start -= k;
      const int atom_num = k + start;
      if (d_pttopo.atomNum(n) != atom_num) {
        ostringstream os;
        os << "Error in PERTPOLPARAM block: The atom " << atom_num + 1 << " was "
           << "not found in perturbation topology. Found atom " << d_pttopo.atomNum(n)
           << " instead.";
        throw InPtTopology::Exception(os.str());
      }

      double alpha[2], dampingLevel[2];
      _lineStream >> l >> nm >> alpha[0] >> dampingLevel[0] >> alpha[1]
                  >> dampingLevel[1];

      if (_lineStream.fail()) {
        ostringstream os;
        os << "Bad line in PERTPOLPARAM block (line " << n + 1 << ").";
        throw InPtTopology::Exception(os.str());
      }

      for(unsigned int i = 0; i < 2; ++i) {
        d_pttopo.setPolarisability(n, i, alpha[i]);
        d_pttopo.setDampingLevel(n, i, dampingLevel[i]);
      }
    }
    if (n != numat)
      throw InPtTopology::Exception("Perturbation topology file " + name() +
            " is corrupted. Failed to read all atoms");
      
    } // PERTPOLPARAM
  } else if(d_blocks.count("MPERTATOM")){
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
	 _lineStream >> iiac[0] >> dq[0];
	 d_pttopo.setIac(n,j,iiac[0]-1);
	 d_pttopo.setCharge(n,j,dq[0]);
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





