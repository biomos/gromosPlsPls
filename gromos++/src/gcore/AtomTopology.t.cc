// gcore_AtomTopology.t.cc

#include "AtomTopology.h"
#include "Exclusion.h"
#include <iostream>

using namespace gcore;
using namespace std;

int main() {
  Exclusion e;
  e.insert(1);
  e.insert(2);
  e.insert(4);

  AtomTopology at;
  at.setName("A");
  at.setIac(1);
  at.setChargeGroup(3);
  at.setMass(.1);
  at.setCharge(.2);
  at.setExclusion(e);

  cout << "Name: " << at.name() << endl
          << "IAC: " << at.iac() << endl
          << "ChargeGroup: " << at.chargeGroup() << endl
          << "mass: " << at.mass() << endl
          << "charge: " << at.charge() << endl
          << at.exclusion().size() << " Exclusions: ";


  for (int i = 0; i < at.exclusion().size(); ++i)
    cout << at.exclusion().atom(i) << ' ';
  cout << endl;
  AtomTopology bt = at;

  cout << "Name: " << bt.name() << endl
          << "IAC: " << bt.iac() << endl
          << "ChargeGroup: " << at.chargeGroup() << endl
          << "mass: " << bt.mass() << endl
          << "charge: " << bt.charge() << endl
          << at.exclusion().size() << " Exclusions: ";


  for (int i = 0; i < at.exclusion().size(); ++i)
    cout << at.exclusion().atom(i) << ' ';
  cout << endl;

  bt.setName("B");
  bt.setCharge(.3);
  at = bt;

  cout << "Size: " << sizeof (at) << ' ' << sizeof (bt) << endl;

  cout << "Name: " << at.name() << endl
          << "IAC: " << at.iac() << endl
          << "ChargeGroup: " << at.chargeGroup() << endl
          << "mass: " << at.mass() << endl
          << "charge: " << at.charge() << endl
          << at.exclusion().size() << " Exclusions: ";


  for (int i = 0; i < at.exclusion().size(); ++i)
    cout << at.exclusion().atom(i) << ' ';
  cout << endl;

  return 0;
}

/* OUTPUT:
Name: A
IAC: 1
ChargeGroup: 3
mass: 0.1
charge: 0.2
3 Exclusions: 1 2 4 
Name: A
IAC: 1
ChargeGroup: 3
mass: 0.1
charge: 0.2
3 Exclusions: 1 2 4 
Size: 4 4
Name: B
IAC: 1
ChargeGroup: 3
mass: 0.1
charge: 0.3
3 Exclusions: 1 2 4 
 */


