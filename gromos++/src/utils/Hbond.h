// utils_Hbond.h

// Class that contains a sequential list of specific atoms

#ifndef INCLUDED_HBOND
#define INCLUDED_HBOND

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

  /**
   * Class Hbond
   * purpose: contains specific atoms of the system, keeping track of 
   * molecule and atom numbers
   *
   * Description:
   * The Hbondcalc can be used to look over a specific set of atoms,
   * possibly spanning different molecules. A 'specifier' is a string with 
   * following format <mol>[:<atom>[-<atom>]]. For example "1:3-5" means atoms
   * 3,4 and 5 of molecule 1; "2:5" means atom 5 of molecule 2; "3" means all
   * atoms of molecule 3.
   *
   * @class Hbondcalc
   * @author M. Kastenholz
   * @ingroup utils
   */
  class Hbond{
    
    double d_mean_dist, d_mean_angle, d_dist, d_angle;
    int d_don_ind, d_b_ind, d_ac_ind;
    int d_num;
   
  public: 
    // Constructor
    Hbond() { 
      d_dist = 0;
      d_angle = 0;
      d_num = 0;
      };

    Hbond(int indx, int indy, int indz) {
      setIndices(indx, indy, indz);
    };
     

    /**
     * Hbondcalc Deconstructor
     */
    ~Hbond(){}
   

    //    int addType(std::string s);
    /**
     * Method to sort the atoms ascending order. Some applications might
     * need the atoms to be ordered. This is just a simple bubble sort
     */
    //    void sort();


    /**
     * Accessor, returns the molecule number of the i-th atom in the
     * Hbondcalc
     */    
     void adddistance(double i);

     void addangle(double i);

     void setIndices(int donor, int boundto, int acceptor);

     void addnum();
     
     int don_ind();

     int b_ind();

     int ac_ind();
   
     int num();

     double meandist();

     double meanangle();

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
      Exception(const std::string &what): gromos::Exception("Hbond", what){}
    };
  protected:
    //Internal function
   
    


  }; //end class Hbond

  inline int Hbond::don_ind()
    {
      return d_don_ind;
    }
  
  inline int Hbond::b_ind()
    {
      return d_b_ind;
    }
  
  inline int Hbond::ac_ind()
    {
      return d_ac_ind;
    }

  inline int Hbond::num()
    {
      return d_num;
    }
  
  inline double Hbond::meandist()
    {
      return d_mean_dist;
    }
  
  inline double Hbond::meanangle()
    {
      return d_mean_angle;
    }
  
}
#endif
