// gcore_Improper.h

#ifndef INCLUDED_GCORE_IMPROPER
#define INCLUDED_GCORE_IMPROPER

namespace gcore{
  /**
   * Class Improper
   * Purpose: contains a gromos96 Improper dihedral angle
   *
   * Description:
   * Contains the atoms and type making up an improper dihedral angle. The 
   * aotms are sorted in such a way that the two middle atoms are in 
   * ascending order (similar to normal Dihedral)
   *
   * @class Improper
   * @author R. Buergi
   * @ingroup gcore
   * @sa gcore::ImproperType
   * @sa gcore::MoleculeTopology
   */
class Improper{
  int d_a[4];
  int d_type;
  // not implemented
  Improper();
 public:
  /**
   * Improper constructor
   * @param a, b, c, d atom numbers defining the improper. Will be stored 
   *                   so that b<c (i.e. abcd of dcba)
   */
  Improper(int a, int b, int c, int d, bool warn = true);
  /**
   * Improper copy constructor
   * @param & Improper to be copied
   */
  Improper(const Improper &);
  /**
   * Improper deconstructor
   */
  ~Improper(){}
  /**
   * Method to set the gromos96 type of the Improper
   */
  void setType(int i){d_type=i;}
  /**
   * Member operator = copies one Improper into the other
   */
  Improper &operator=(const Improper &);
  /**
   * Accessor, returns the atom number of the i-th atom in the Improper
   * as a const
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Accessor, returns the atom number of the i-th atom in the Improper
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the type of the Improper
   */
  int type()const{return d_type;}
};
/**
 * @relates Improper
 * Operator < compares two improper angles a and b
 *
 * a is defined to be smaller than b if <ol>
 * <li> The atom number of the second atom in Improper a is smaller than
 *      the atom number of the second atom in Improper b (a[1]<b[1]), or
 * <li> The second atoms in a and b have the same atom number, but the 
 *      third atom number in a is lower than the third atom number in b
 *      (a[1]==b[1] && a[2]<b[2]), or
 * <li> The middle atom numbers of a and b are all the same, but the first
 *      atom number in a is lower than the first atom number in b
 *      (a[1]==b[1] && a[2]==b[2] && a[0]<b[0]), or
 * <li> The first three atom numbers in a and b are all the same, but the
 *      fourth atom number in a is lower than the fourth atom number in b
 *      (a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]<b[3]), or  
 * <li> All atom numbers in a and b are the same, but a has a lower bond
 *      type than b (a[0]==b[0] && a[1]==b[1] && a[2]==b[2] && a[3]==b[3] 
 *      && a.type()<b.type())
 *</ol>
 */
int operator<(const Improper &a, const Improper &b);

}
#endif
