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
    int d_Hmol, d_Htomol, d_Acmol;
    int d_Hatom, d_Htoatom, d_Acatom;
    int d_don_ind, d_b_ind, d_ac_ind;
    int d_num;
   
  public: 
    // Constructor
    /**
     * Hbondcalc Constructor
     * @param sys The Hbondcalc needs to know about the system. It 
     *            does not know about any atoms yet.
     * @param args all arguments are passed into Hbondcalc. 
     */
    Hbond() { 
      d_dist = 0;
      d_angle = 0;
      d_num = 0;
      };

    Hbond(int i, int j, int k, int l, int m, int n, int indx, int indy, int indz) {
      setHBond(i,j,k,l,m,n);
      setindices(indx, indy, indz);
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

     void setHBond(int Hmol, int Hatom, int Htomol, int Htoatom, int Acmol, int Acatom);

     void setindices(int donor, int boundto, int acceptor);

     void addnum();
     
     int don_ind();

     int b_ind();

     int ac_ind();
   
     int num();

     int Hmol();
   
     int Htomol();
  
     int Acmol();

     int Hatom(); 
  
     int Htoatom();

     int Acatom();

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

  inline int Hbond::Hmol()
  {
  return d_Hmol;
  }

  inline int Hbond::Htomol()
  {
  return d_Htomol;
  }

  inline int Hbond::Acmol()
  {
  return d_Acmol;
  }

  inline int Hbond::Hatom()
  {
  return d_Hatom;
  }

  inline int Hbond::Htoatom()
  {
  return d_Htoatom;
  }

  inline int Hbond::Acatom()
  {
  return d_Acatom;
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
