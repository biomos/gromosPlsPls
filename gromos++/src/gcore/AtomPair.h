#ifndef INCLUDED_ATOMPAIR
#define INCLUDED_ATOMPAIR

namespace gcore{

  /**
   * Class AtomPair
   * Purpose, define a pair of atoms
   *
   * Description:
   * stores two atoms, that make up a pair. The atoms are sorted in such 
   * a way that the atom with the lowest atom number comes first.
   * 
   * @class AtomPair
   * @ingroup gcore
   * @author R. Buergi
   */
class AtomPair{
  int d_a[2];
  // not implemented
  AtomPair();
 public:
  /**
   * AtomPair constructor. The atoms are sorted in such a way that the 
   * atom with the lowest atom number comes first.
   * @param a,b atom numbers of atoms that make up a pair
   */
  AtomPair(int a, int b);
  /**
   * AtomPair copyconstructor
   * @param & AtomPair to be copied
   */
  AtomPair(const AtomPair &);
  /**
   * Member operator[]
   * @param i element of the atom pair to be returned (0 or 1)
   * @return The atom number of atom i in the pair
   */
  int operator[](int i)const{return d_a[i];}
};
/**
 * @relates AtomPair
 * Operator < compares to AtomPairs a and b<br>
 * a is defined to be smaller than b if<br><ol>
 * <li> the first atom of pair a has a lower atomnumber than the first 
 *      atom of pair b (a[0]<b[0]), or 
 * <li> the first atoms of both pairs have the same atomnumbers and the
 *      second atom of pair a has the lowest atomnumber
 * </ol>
 * @param a, b atomPairs to be compared
 * @return 1 if a < b and 0 otherwise
 */
int operator<(const AtomPair &a, const AtomPair &b);
}
#endif
