// gcore_BuildingBlock.h

#ifndef INCLUDED_GCORE_BUILDINGBLOCK
#define INCLUDED_GCORE_BUILDINGBLOCK

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif

namespace gcore{

  class BbSolute;
  class BbEnd;
  class SolventTopology;

  class BuildingBlock{
    std::vector<BbSolute*> d_bb;
    std::vector<BbEnd*> d_be;
    std::vector<SolventTopology*> d_bs;
    double d_fpepsi;
    double d_hbar;
    int d_linkExclusions;
    

  public:
    //Constructors
    BuildingBlock();
    BuildingBlock(const BuildingBlock &bld);
    ~BuildingBlock();

    // Methods
    BuildingBlock &operator=(const BuildingBlock &bld);
    void addBbSolute(const BbSolute &mol);
    void addBbSolvent(const SolventTopology &sol);
    void addBbEnd(const BbEnd &mol);
    void setFpepsi(const double a){d_fpepsi=a;};
    void setHbar(const double a){d_hbar=a;};
    void setLinkExclusions(const int i){d_linkExclusions=i;};
    
    
    // Accessors
    const BbSolute &bb(int i)const;
    BbSolute &bb(int i);
    const BbEnd &be(int i)const;
    BbEnd &be(int i);
    const SolventTopology &bs(int i)const;
    SolventTopology &bs(int i);

    int numBbSolutes()const;
    int numBbSolvents()const;
    int numBbEnds()const;
    double Fpepsi()const;
    double Hbar()const;
    int LinkExclusions()const;
    int findBb(std::string s);
    int findBs(std::string s);
    
    
};

  inline const BbSolute &BuildingBlock::bb(int i)const{
    assert (i < this->numBbSolutes());
    return *d_bb[i];
  }

  inline BbSolute &BuildingBlock::bb(int i){
    assert (i < this->numBbSolutes());
    return *d_bb[i];
  }

  inline const BbEnd &BuildingBlock::be(int i)const{
      assert (i < this->numBbEnds());
      return *d_be[i];
  }

  inline BbEnd &BuildingBlock::be(int i){
      assert (i < this->numBbEnds());
      return *d_be[i];
  }

  inline const SolventTopology &BuildingBlock::bs(int i)const{
      assert (i< this->numBbSolvents());
      return *d_bs[i];
  }

  inline SolventTopology &BuildingBlock::bs(int i){
      assert (i < this->numBbSolvents());
      return *d_bs[i];
  }

  inline int BuildingBlock::numBbSolutes()const{
    return d_bb.size();
  }

  inline int BuildingBlock::numBbEnds()const{
      return d_be.size();
  }

  inline int BuildingBlock::numBbSolvents()const{
    return d_bs.size();
  }

  inline double BuildingBlock::Fpepsi()const{
    return d_fpepsi;
  }
  
  inline double BuildingBlock::Hbar()const{
    return d_hbar;
  }
  
  inline int BuildingBlock::LinkExclusions()const{
    return d_linkExclusions;
  }
  
  
}
#endif

