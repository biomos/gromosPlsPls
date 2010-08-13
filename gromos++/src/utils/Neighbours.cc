#include <cassert>
#include <set>

#include "../gcore/System.h"
#include "../gcore/Molecule.h"
#include "../gcore/LJException.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/SolventTopology.h"
#include "../gcore/Solvent.h"
#include "../gcore/Bond.h"
#include "../gcore/Constraint.h"

#include "Neighbours.h"

using namespace gcore;
using namespace std;
using utils::Neighbours;

static void init(Neighbours *n, const MoleculeTopology &mt, int i){
  BondIterator iter(mt);
  for(;iter;++iter){
    if(iter()[0]==i)n->push_back(iter()[1]);
    if(iter()[1]==i)n->push_back(iter()[0]);
  }
}

static void inits(Neighbours *n, const SolventTopology &mt, int i){
  ConstraintIterator iter(mt);
  for(;iter;++iter){
    if(iter()[0]==i)n->push_back(iter()[1]);
    if(iter()[1]==i)n->push_back(iter()[0]);
  }
}

Neighbours::Neighbours(const System &sys, int mol, int i):
  vector<int>()
{
  //  int j=0, atNum=0;
  //  while(i > (atNum+=sys.mol(j).numAtoms())) ++j;
  //  int k=i-atNum+sys.mol(j).numAtoms();

  init(this,sys.mol(mol).topology(),i);
}

Neighbours::Neighbours(const System &sys, int mol, int i, int j):
  vector<int>()
{
  //  int j=0, atNum=0;
  //  while(i > (atNum+=sys.mol(j).numAtoms())) ++j;
  //  int k=i-atNum+sys.mol(j).numAtoms();

  inits(this,sys.sol(mol).topology(),i);
}

Neighbours::Neighbours(const MoleculeTopology &mt, int i): vector<int>()
{
  init(this,mt,i);
}

Neighbours::Neighbours(const Molecule &mol, int k):
  vector<int>()
{
  init(this,mol.topology(),k);
}





