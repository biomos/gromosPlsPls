// utils_FfExpert.h

// Class that contains statistically information on the occurrence
// of types and charges in a force field building block

#ifndef INCLUDED_UTILS_FFEXPERT
#define INCLUDED_UTILS_FFEXPERT

namespace gcore
{
  class BuildingBlock;
  class Bond;
  class Angle;
  class Improper;
  class Dihedral;
}

namespace utils
{
  /**
   * Class FfExpert
   * contains statistical data on the occurence of types in a BuildingBlock
   * usefull for checking of consistency and suggesting types.
   *
   * Description:
   * This is a low-level expert system that knows which types of bonds, 
   * angles etc. are most commonly seen with certain atomic Iac values. The 
   * IAC are connected to the first letters of atom names
   *
   * @class FfExpert
   * @author C. Oostenbrink
   * @ingroup utils
   * @sa gcore::BuildingBlock
   */
  class FfExpert{
  public:
    struct counter
    {
      int type, occurence;
      counter(int i, int j)
      {
	type=i;
	occurence = j;
      }
      
    };
  private:
    std::multimap<std::string,counter> d_name2iac;
    std::multimap<int,counter> d_iac2mass;
    std::multimap<int,counter> d_iac2charge;
    std::multimap<gcore::Bond,counter> d_iac2bond;
    std::multimap<gcore::Angle,counter> d_iac2angle;
    std::multimap<gcore::Improper,counter> d_iac2improper;
    std::multimap<gcore::Dihedral,counter> d_iac2dihedral;
    std::vector<double> d_chargeType;
  public:
    /**
     * Standard constructor
     */
    FfExpert(){};
    /**
     * Constructor with learning of a file
     */
    FfExpert(gcore::BuildingBlock const & mtb);
    /**
     * Function to learn about a BuildingBlock
     */
    void learn(gcore::BuildingBlock const & mtb);
    /**
     * Accessor to the names
     */
    void name2iac(std::string s, std::vector<counter> &v);
    /**
     * Accessor to the masses
     */
    void iac2mass(int i, std::vector<counter> &v);
    /**
     * Accessor to charge types
     */
    void iac2charge(int i, std::vector<counter> &v);
    /**
     * Accessor to charge as function of type
     */
    double charge(int i);
    /**
     * Accessor to bonds
     */
    void iac2bond(gcore::Bond const & b, std::vector<counter> &v);
 /**
     * Accessor to angles
     */
    void iac2angle(gcore::Angle const & b, std::vector<counter> &v);
    /**
     * Accessor to impropers
     */
    void iac2improper(gcore::Improper const & b, std::vector<counter> &v);
 /**
     * Accessor to dihedrals
     */
    void iac2dihedral(gcore::Dihedral const & b, std::vector<counter> &v);
 

  };

  static int sort(std::vector<FfExpert::counter> &v, bool tt=true);


  // Inline functions and methods
  inline FfExpert::FfExpert(gcore::BuildingBlock const & mtb){
      learn(mtb);
  }
  inline void FfExpert::name2iac(std::string s, std::vector<counter> &v)
    {
      v.clear();
      if(d_name2iac.count(s)) {
	for(std::multimap<std::string, FfExpert::counter>::const_iterator iter=d_name2iac.lower_bound(s), to=d_name2iac.upper_bound(s); iter!=to; ++iter){
	  v.push_back(iter->second);
	}
      }
    }
  
	
  inline void FfExpert::iac2mass(int i, std::vector<counter> &v)
    {
      v.clear();
      if(d_iac2mass.count(i)){
	for(std::multimap<int, FfExpert::counter>::const_iterator iter=d_iac2mass.lower_bound(i), to=d_iac2mass.upper_bound(i); iter!=to; ++iter){
	  v.push_back(iter->second);
	}
      }
    }
  inline void FfExpert::iac2charge(int i, std::vector<counter> &v)
    {
      v.clear();
      if(d_iac2charge.count(i)){
	for(std::multimap<int, FfExpert::counter>::const_iterator iter=d_iac2charge.lower_bound(i), to=d_iac2charge.upper_bound(i); iter!=to; ++iter){
	  v.push_back(iter->second);
	}
      }
    }
  inline double FfExpert::charge(int i)
  {
    return d_chargeType[i];
  }
  
  inline void FfExpert::iac2bond(gcore::Bond const & b, 
				 std::vector<counter> &v)
  {
    v.clear();
    if(d_iac2bond.count(b)){
      for(std::multimap<gcore::Bond, FfExpert::counter>::const_iterator 
	    iter=d_iac2bond.lower_bound(b), to=d_iac2bond.upper_bound(b); 
	  iter!=to; ++iter){
	v.push_back(iter->second);
      }
    }
  }
  
  inline void FfExpert::iac2angle(gcore::Angle const & b, 
				 std::vector<counter> &v)
  {
    v.clear();
    if(d_iac2angle.count(b)){
      for(std::multimap<gcore::Angle, FfExpert::counter>::const_iterator 
	    iter=d_iac2angle.lower_bound(b), to=d_iac2angle.upper_bound(b); 
	  iter!=to; ++iter){
	v.push_back(iter->second);
      }
    }
  }  

  inline void FfExpert::iac2improper(gcore::Improper const & b, 
				 std::vector<counter> &v)
  {
    v.clear();
    if(d_iac2improper.count(b)){
      for(std::multimap<gcore::Improper, FfExpert::counter>::const_iterator 
	    iter=d_iac2improper.lower_bound(b), to=d_iac2improper.upper_bound(b); 
	  iter!=to; ++iter){
	v.push_back(iter->second);
      }
    }
  }  

  inline void FfExpert::iac2dihedral(gcore::Dihedral const & b, 
				 std::vector<counter> &v)
  {
    v.clear();
    if(d_iac2dihedral.count(b)){
      for(std::multimap<gcore::Dihedral, FfExpert::counter>::const_iterator 
	    iter=d_iac2dihedral.lower_bound(b), to=d_iac2dihedral.upper_bound(b); 
	  iter!=to; ++iter){
	v.push_back(iter->second);
      }
    }
  }  

 int utils::sort(std::vector<FfExpert::counter> &v, bool tt)
  {
    int max_occur=0, max_index=0;
    if(tt){
      for(unsigned int i=1; i<v.size(); i++){
	
	FfExpert::counter t=v[i];
	int j=i-1;
	while ((j>=0) && t.type < v[j].type){
	  v[j+1]=v[j];
	  j--;
	}
	v[j+1]=t;
      }
      for(unsigned int i=0; i< v.size(); i++){
	if(v[i].occurence > max_occur){
	  max_occur=v[i].occurence;
	  max_index=i;
	}
      }
      return max_index;
    }
    else{
      for(unsigned int i=1; i<v.size(); i++){
	FfExpert::counter t=v[i];
	int j=i-1;
	while ((j>=0) && t.occurence > v[j].occurence){
	  v[j+1]=v[j];
	  j--;
	}
	v[j+1]=t;
      }
      return 0;
    }
  }

} // namespace utils

#endif