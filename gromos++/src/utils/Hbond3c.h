// utils_Hbond.h

// Class that contains a sequential list of specific atoms

#ifndef INCLUDED_HBOND3C
#define INCLUDED_HBOND3C

#ifndef INCLUDED_VECTOR
#include <vector>
#define INCLUDED_VECTOR
#endif
#ifndef INCLUDED_STRING
#include <string>
#define INCLUDED_STRING
#endif
#ifndef INCLUDED_GROMOS_EXCEPTION
#include "../gromos/Exception.h"
#endif







namespace utils
{

  class Hbond3c{
    double d_dist_don_A1, d_dist_don_A2, d_angle_b_don_A1, d_angle_b_don_A2,
      d_angle_sum, d_dihedral;
    double d_mean_dist_don_A1, d_mean_dist_don_A2, d_mean_angle_b_don_A1, d_mean_angle_b_don_A2,
      d_mean_angle_sum, d_mean_dihedral;
    int d_don_ind, d_b_ind, d_ac1_ind, d_ac2_ind;
    int d_num;
   
  public: 
    // Constructor
    Hbond3c() { 
      d_dist_don_A1 = 0;
      d_dist_don_A2 = 0;
      d_angle_b_don_A1 = 0;
      d_angle_b_don_A2 = 0;
      d_angle_sum = 0;
      d_dihedral = 0;
      d_num = 0;
      };

    Hbond3c(int ind_don, int ind_b, int ind_ac1, int ind_ac2) {
      setIndices(ind_don, ind_b, ind_ac1, ind_ac2);
    };
     

    /**
     * Hbondcalc Deconstructor
     */
    ~Hbond3c(){}
    
    
    void adddistances(double i, double j);
    
    void addangles(double i, double j);
    
    void addanglesum(double i);
    
    void adddihedral(double i);

    void setIndices(int donor, int boundto, int acceptor1, int acceptor2);
    
    void addnum();
    
    int don_ind();
    
    int b_ind();
    
    int ac1_ind();
    
    int ac2_ind();

    int num();
    
    double meandist_don_a1();
    double meandist_don_a2();
    
    double meanangle_b_don_a1();
    double meanangle_b_don_a2();
    
    double meanangle_sum();
    double meandihedral();
    
    void calcmean();

    void clear();

    /**
     * @struct Exception
     * Throws an exception if something is wrong
     */
    struct Exception: public gromos::Exception{
      /**
       * @exception If called says Hbondcalc, followed by the argument
       * @param what The string that is thrown
       */
      Exception(const std::string &what): gromos::Exception("Hbond3c", what){}
    };
  protected:
    //Internal function
    
    
    

  }; //end class Hbond

  inline int Hbond3c::don_ind()
    {
      return d_don_ind;
    }
  
  inline int Hbond3c::b_ind()
    {
      return d_b_ind;
    }
  
  inline int Hbond3c::ac1_ind()
  {
    return d_ac1_ind;
  }
  inline int Hbond3c::ac2_ind()
  {
    return d_ac2_ind;
  }
  

  inline int Hbond3c::num()
  {
    return d_num;
  }
  
  inline double Hbond3c::meandist_don_a1()
  {
    return d_mean_dist_don_A1;
  }
  inline double Hbond3c::meandist_don_a2()
  {
    return d_mean_dist_don_A2;
  }
  
  inline double Hbond3c::meanangle_b_don_a1()
  {
    return d_mean_angle_b_don_A1;
  }

  inline double Hbond3c::meanangle_b_don_a2()
  {
    return d_mean_angle_b_don_A2;
  }
  inline double Hbond3c::meanangle_sum()
  {
    return d_mean_angle_sum;
  }
  inline double Hbond3c::meandihedral()
  {
    return d_mean_dihedral;
  }
  
}
#endif
