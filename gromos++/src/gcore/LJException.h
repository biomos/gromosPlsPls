// gcore_LJException.h

#ifndef INCLUDED_GCORE_LJEXCEPTION
#define INCLUDED_GCORE_LJEXCEPTION

namespace gcore{

  /**
   * Class LJException
   * Purpose: contains a LJ exception
   *
   * Description:
   * Contains the atoms and type making up a LJException. The atoms are sorted
   * to have the lowest atom number first.
   *
   * @class LJException
   * @author A. Eichenberger
   * @ingroup gcore
   */
class LJException{
  int d_a[2];
  int d_type;
  // not implemented
  LJException();
 public:
  /**
   * LJException constructor
   * constructs a new LJException defined by atoms a and b. These atoms are stored
   * such that a<b.
   * @param a,b atom numbers of the atoms making the LJException
   */
  LJException(int a, int b);
  /**
   * LJException copy constructor
   * Constructs a new LJException and copies the specied LJException into it
   * @param & LJException to be copied
   */
  LJException(const LJException &);
  /**
   * LJException deconstructor
   */
  ~LJException(){}
  /**
   * Member operator =, copies one LJException into the other
   */
  LJException &operator=(const LJException &b);
  /**
   * Method setType sets the LJException type of the LJException
   * @param i LJException type to be set
   */
  void setType(int i){d_type=i;}
  /**
   * Member operator [], returns the atom number of the i-th atom in the LJException
   * @param i atom index in the LJException (0 or 1)
   * @return Atom number of the i-th atom in the LJException
   */
  int operator[](int i)const{return d_a[i];}
  /**
   * Member operator [], returns the atom number of the i-th atom in the LJException
   * @param i atom index in the LJException (0 or 1)
   * @return Atom number of the i-th atom in the LJException
   */
  int &operator[](int i){return d_a[i];}
  /**
   * Accessor, returns the LJException type of the LJException
   */
  int type()const{return d_type;}
};
/**
 * @relates LJException
 * Operator < compares two LJExceptions a and b
 *
 * a is defined to be smaller than b if<br><ol>
 * <li> The atom number of the first atom in LJException a is lower than the
 *      atom number of the first atom in LJException b (a[0]<b[0]), or
 * <li> The atom numbers of the first atoms in both LJExceptions are the same
 *      but the second atom has a lower atom number in LJException a (a[0]==b[0]
 *      && a[1]<b[1]), or
 * <li> All atoms are the same, but the LJException type of LJException a has a lower
 *      number than the LJException type of LJException b (a[0]==b[0] && a[1]==b[1] &&
 *      a.type() < b.type()).
 * </ol>
 * &param a, b LJExceptions to be compared
 * &return 1 if a<b; 0 otherwise
 */
int operator<(const LJException &a, const LJException &b);
}
#endif
