#include "MoleculeTopology.h"
#include "BbEnd.h"
#include "AtomTopology.h"
#include "Exclusion.h"
#include "Bond.h"
#include "Angle.h"
#include <iostream>

using namespace gcore;
using namespace std;

int main(){
  BbEnd mt;
  AtomTopology at;
  at.setIac(3);
  mt.addAtom(at);
  at.setIac(4);
  mt.addAtom(at);
  mt.setRep(2);
  
  cout << "AtomTopology IAC: ";
  for(int i=0; i<mt.numAtoms();++i)
    cout << mt.atom(i).iac() << ' ';
  cout << endl;
  cout << mt.rep();
  cout << endl;
  
  
  Bond b(1,2);
  mt.addBond(b);
  b=Bond(2,3);
  mt.addBond(b);
  
  cout << "Bonds: ";
  BeBondIt bi(mt);
  for(;bi;++bi)
    cout << '(' << bi()[0] << ',' << bi()[1] << ") ";
  cout << endl;
  
  Angle ang(1,2,3);
  mt.addAngle(ang);
  ang=Angle(4,5,6);
  mt.addAngle(ang);

  cout << "Angles: ";
  BeAngleIt ai(mt);
  for(;ai;++ai)
    cout << '(' << ai()[0] << ' ' << ai()[1] << ' ' << ai()[2] << ") ";
  cout << endl;

  mt.setResName("ASP");
  cout << "Resnames: ";
  for(int i=0; i<5; ++i)
    cout << i << ' '<< mt.resName() << ' ';
  cout << endl;
  
  return 0;
}
