// neigbour.cc

#include "../src/args/Arguments.h"
#include "../src/args/BoundaryParser.h"
#include "../src/args/GatherParser.h"
#include "../src/gio/InG96.h"
#include "../src/gcore/System.h"
#include "../src/gcore/Molecule.h"
#include "../src/gcore/MoleculeTopology.h"
#include "../src/gcore/Solvent.h"
#include "../src/gio/InTopology.h"
#include "../src/bound/TruncOct.h"
#include "../src/bound/RectBox.h"
#include "../src/bound/Vacuum.h"
#include "../src/gmath/Matrix.h"
#include "../src/utils/AtomSpecifier.h"
#include <vector>
#include <set>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <strstream>

using namespace gcore;
using namespace gio;
using namespace gmath;
using namespace bound;
using namespace args;
using namespace std;
using namespace utils;

class potential
{
  int d_mol;
  int d_atom;
  Vec d_pos;
  
public:
  potential(int mol, int atom);
  ~potential(){}
  potential &operator=(const potential &p);
  Vec setPos(Vec a, Vec b);
  //accessors  
  int mol()const {return d_mol;}
  int atom()const {return d_atom;}
  Vec pos()const {return d_pos;}
};

class vertex
{
  Vec d_pos;
  int d_atoms[3];

public:
  vertex(int a, int b, int c);
  ~vertex(){}
  
  int calcVert(Vec a, Vec b, Vec c);
  vertex &operator=(const vertex &v);
  //accessors
  Vec pos() { return d_pos; }
  int operator[](int i) {return d_atoms[i];}
  int atoms(int i){return d_atoms[i];}
};

void scan_atoms(System &sys, int m, int a, double cut, vector <potential> &p);
void calc_dummies(Vec v, vector <potential> &p);
void init_vertices(vector <potential> p, vector <vertex> &v);
void calc_vertices(vector <potential> p, vector <vertex> &v);
set <int> get_neighbours(vector <vertex> v);



int main(int argc, char **argv){

  char *knowns[] = {"topo", "traj", "nsm", "pbc", "time", "cut", "atoms", "out"};
  int nknowns = 8;

  string usage = argv[0];
  usage += "\n\t@topo <topology>\n";
  usage += "\t@pbc  <boundary type>\n";
  usage += "\t@nsm  <number of molecules in traj>\n";
  usage += "\t@time <time> <dt>\n";
  usage += "\t@cut  <initial cut-off>\n";
  usage += "\t@atoms <atomspecifier>\n";
  usage += "\t@traj <trajectory files>\n";
  usage += "\t@out  <filename for output>\n";

  try{
    Arguments args(argc, argv, nknowns, knowns, usage);

    // get the initial cut-off
    double cut=0.5, max=0.0;
    {
      Arguments::const_iterator iter=args.lower_bound("cut");
      if(iter!=args.upper_bound("cut")){
	cut=atof(iter->second.c_str());
      }
    }
    
    // get simulation time
    double time=0, dt=1;
    int frame=0;
    
    {
      Arguments::const_iterator iter=args.lower_bound("time");
      if(iter!=args.upper_bound("time")){
	time=atof(iter->second.c_str());
	++iter;
      }
      if(iter!=args.upper_bound("time"))
	dt=atof(iter->second.c_str());
    }
    // read topology
    InTopology it(args["topo"]);
    System sys(it.system());

    // get number of solute molecules
    int nsm=1;
    
    {
      Arguments::const_iterator iter=args.lower_bound("nsm");
      if(iter!=args.upper_bound("nsm")){
	nsm=atoi(iter->second.c_str());
      }
    }
    for(int i=1;i<nsm;i++){
      MoleculeTopology mt=sys.mol(0).topology();
      Molecule mol(mt);
      
      sys.addMolecule(mol);
    }
    int res[nsm*sys.numMolecules()][100];
    for(int i=0;i<100;i++)
      for(int j=0;j<nsm*sys.numMolecules();j++)
	res[j][i]=0;
    
    // get atomspecifier
    utils::AtomSpecifier which(sys);
    {
      Arguments::const_iterator iter=args.lower_bound("atoms");
      Arguments::const_iterator to  =args.upper_bound("atoms");
      if(iter==to){
	for(int i=0;i<sys.numMolecules();i++)
            for(int j=0;j<sys.mol(i).numAtoms();j++)
              which.addAtom(i,j);
      }
      else{
	for(;iter!=to;iter++) {
          string s=iter->second.c_str();
          if (s=="all") {
            for(int i=0;i<sys.numMolecules();i++)
              for(int j=0;j<sys.mol(i).numAtoms();j++)
                which.addAtom(i,j);
            break;
  	  }
          else
            which.addSpecifier(s);
	}
      }
    }
    
    // get output filename
    int file=0;
    ofstream fout;
    {
      Arguments::const_iterator iter=args.lower_bound("out");
      if(iter!=args.upper_bound("out")){
	file=1;
	fout.open(iter->second.c_str());
      }
    }
    if(file) fout << "nr_atoms " << which.size() << endl;
    
    
    // Parse boundary conditions
    Boundary *pbc = BoundaryParser::boundary(sys, args);
    // parse gather method
    Boundary::MemPtr gathmethod = args::GatherParser::parse(args);

    
    // loop over the trajectories
    

    // write out some info
    cout << "atoms are specified as mol:atom\n";
    cout << "if mol = 0 this implies a solvent atom\n\n";
    
    for(Arguments::const_iterator iter=args.lower_bound("traj");
    	  iter!=args.upper_bound("traj"); ++iter){
      InG96 ic;    
      ic.open(iter->second.c_str());
      ic.select("ALL");
      
      // loop over all frames
      while(!ic.eof()){ 
        cout << "Neighbour list at time : " << time << endl;
	if(file) {
          fout << "time " << time << endl;
	  fout << "atom\tnr_neighbours\tneighbours\n";
	}
	
        ic >> sys;

        // loop over all specified atoms
        for(int spec=0;spec<which.size();spec++){
	  int i=which.mol(spec);
	  int j=which.atom(spec);

          vector<potential> p;
          Vec curr=sys.mol(i).pos(j);
	  // gather around this atom
          pbc->setReference(0,curr);
	  (*pbc.*gathmethod)();
	   
	  // generate dummy atoms to build a first polyhedron
          calc_dummies(curr, p);
	  
	  // get the possible neighbours
          scan_atoms(sys, i, j, cut, p);
	  
          // now, to get the vertices.
          vector<vertex> vert;
          // The dummy atoms create four vertices
          init_vertices(p, vert);
	    
          // now do the real atoms
          calc_vertices(p, vert);

          // loop over remaining vertices and print out the neighbours
	  set<int> neigh=get_neighbours(vert);

	  // print out neighbours
	  cout << "\t" << i+1 << ":" << j+1 << "|";
          if(file) fout << i+1 << ":" << j+1 << "\t" << neigh.size() << "\t";
          int cnt=1;
	    
          for(set<int>::const_iterator itt=neigh.begin(), to=neigh.end();
              itt!=to;itt++){
            cout << "\t" << p[*itt].mol()+1 << ":" << p[*itt].atom()+1;
            if(file) fout << "\t" << p[*itt].mol()+1 << ":" << p[*itt].atom()+1;
            if (p[*itt].pos().abs()>max) max=p[*itt].pos().abs();
            if(cnt++%7==0) {
              cout << "\n\t";
              if(file) fout << "\n\t\t";
	    }
	    
	  }
          cout << endl << endl;
          if(file) fout << endl << endl;
	  
        
        }
        time+=dt;
        frame++;
	
      }
      
      ic.close();
      
    }
    if(cut - max <= 0.01){
      cout << "Warning! You specified a cut-off of " << cut << " nm, ";
      cout << "but the largest neighbour-neighbour distance is " << max;
      cout << "\nConsider rerunning with larger cut-off\n";
    }
  }
  catch (const gromos::Exception &e){
    cerr << e.what() << endl;
    exit(1);
  }
  return 0;
}

void scan_atoms(System &sys, int m, int a, double cut, vector<potential> &p)
{
  Vec curr=sys.mol(m).pos(a);

  // first, the molecules
  int numm=sys.numMolecules();
  for(int k=0;k<numm;k++){
    int numa=sys.mol(k).numAtoms();
    for(int l=0;l<numa;l++){
      if(k!=m || l!=a){
        potential pp(k,l);
        pp.setPos(curr, sys.mol(k).pos(l));
		  
        // do we consider this neighbour
        if(pp.pos().abs()<cut)
	  p.push_back(pp);
      }
    }
  }
  // and now do the same for the solvent
  for(int k=0;k<sys.sol(0).numCoords();k++){
    potential pp(-1,k);
    pp.setPos(curr, sys.sol(0).pos(k));
	      
    // do we consider this neighbour
    if(pp.pos().abs()<cut)
      p.push_back(pp);
  }
}


void calc_dummies(Vec v, vector <potential> &p)
{
  potential dum1(-2,0);
  potential dum2(-2,1);
  potential dum3(-2,2);
  potential dum4(-2,3);
  Vec v1 = v + Vec( 1, 1, 1);
  Vec v2 = v + Vec( 1,-1,-1);
  Vec v3 = v + Vec(-1, 1,-1);
  Vec v4=  v + Vec(-1,-1, 1);
  dum1.setPos(v,v1);
  dum2.setPos(v,v2);
  dum3.setPos(v,v3);
  dum4.setPos(v,v4);
  p.push_back(dum1);
  p.push_back(dum2);
  p.push_back(dum3);
  p.push_back(dum4);
  
}

void init_vertices(vector <potential> p, vector <vertex> &v)
{
  vertex v1(0,1,2);
  v1.calcVert(p[0].pos(), p[1].pos(), p[2].pos());
  vertex v2(1,2,3);
  v2.calcVert(p[1].pos(), p[2].pos(), p[3].pos());
  vertex v3(2,3,0);
  v3.calcVert(p[2].pos(), p[3].pos(), p[0].pos());
  vertex v4(3,0,1);
  v4.calcVert(p[3].pos(), p[0].pos(), p[1].pos());
  v.push_back(v1);
  v.push_back(v2);
  v.push_back(v3);
  v.push_back(v4);
}
	      
void calc_vertices(vector <potential> p, vector <vertex> &vert)
{
  double eps=1e-5;
  
  // loop over all the (real) atoms
  int num=p.size();
    
  for(int a=4; a<num; a++){
    vector<vertex> vadd;
	      
    Vec v=p[a].pos();
    double vabs2=v.abs2();
	      
    //loop over the existing vertices
    for(vector<vertex>::iterator itr=vert.begin();
      itr!=vert.end();itr++){
      double temp=v.dot(itr->pos());
		
      if(temp>0.5*vabs2){
        //add three vertices and remove this one
        vertex v1(a,itr->atoms(0), itr->atoms(1));
	vadd.push_back(v1);
	vertex v2(a,itr->atoms(0), itr->atoms(2));
	vadd.push_back(v2);
        vertex v3(a,itr->atoms(1), itr->atoms(2));
	vadd.push_back(v3);
        vert.erase(itr--);
      }
    }

    // add and calculate added vertices
    int vaddsize=vadd.size();
	      
    for(int q=0;q<vaddsize;q++){
      // check if exists
      int exist=0;
      for(unsigned int r=0;r<vert.size()&&!exist;r++){
	if((vert[r].atoms(0)==vadd[0].atoms(0))&&
           (vert[r].atoms(1)==vadd[0].atoms(1))&&
           (vert[r].atoms(2)==vadd[0].atoms(2)))
	  exist=1;
      }
      if(!exist){
	if(vadd[0].calcVert(p[vadd[0].atoms(0)].pos(),
                            p[vadd[0].atoms(1)].pos(),
	                   p[vadd[0].atoms(2)].pos())==1){
	  int sure=1;

          for(int s=0;s<a;s++){
	    Vec vs=p[s].pos();
	    double tmp=vs.dot(vadd[0].pos());
            double tmp2=vs.abs2();
	    if(tmp>0.5*tmp2+eps) sure=0;
	  }
	  if(sure) vert.insert(vert.end(), vadd[0]);
        }
      }
      vadd.erase(vadd.begin());
    }
  }
}

set <int> get_neighbours(vector<vertex> vert)
{ 
  set<int> neigh;
  
  for(vector<vertex>::iterator itr=vert.begin(), to=vert.end();
     itr!=to;itr++){
    for(int k=0;k<3;k++){
      neigh.insert(itr->atoms(k));
    }
  }
  return neigh;
  
}


potential::potential(int mol, int atom)
{
  d_mol=mol;
  d_atom=atom;
}

potential &potential::operator=(const potential &p)
{
  if(this!=&p){
    d_mol=p.d_mol;
    d_atom=p.d_atom;
    d_pos=p.d_pos;
  }
  return *this;;
}


Vec potential::setPos(Vec a, Vec b)
{
  d_pos=b-a;
  return d_pos;
}

vertex::vertex(int a, int b, int c)
{
  int d=a, e=b, f=c;
  if(a<b&&a<c){
    d=a;
    if(b<c){ e=b; f=c; }
    else   { e=c; f=b; }
  }
  if(b<a&&b<c){
    d=b;
    if(a<c){ e=a; f=c; }
    else   { e=c; f=a; }
  }
  if(c<a&&c<b){
    d=c;
    if(a<b){ e=a; f=b; }
    else   { e=b; f=a; }
  }
  d_atoms[0]=d; d_atoms[1]=e; d_atoms[2]=f;
}

int vertex::calcVert(Vec a, Vec b, Vec c)
{
  double B;
  Vec val(-0.5*a.abs2(),-0.5*b.abs2(), -0.5*c.abs2());
  Vec tmp0(a[0], b[0], c[0]);
  Vec tmp1(a[1], b[1], c[1]);
  Vec tmp2(a[2], b[2], c[2]);
  Matrix m(tmp0, tmp1, tmp2);
  B=m.det();
  if(abs(B)>1e-12){
    m=Matrix(val, tmp1, tmp2);
    d_pos[0]=-m.det()/B;
    m=Matrix(tmp0, val, tmp2);
    d_pos[1]=-m.det()/B;
    m=Matrix(tmp0, tmp1, val);
    d_pos[2]=-m.det()/B;
    return 1;
  }
  else{
    return 0;
  }
}

vertex &vertex::operator=(const vertex &v)
{
  if(this!= &v){
    for(int i=0; i<3;i++) d_atoms[i]=v.d_atoms[i];
    d_pos=v.d_pos;
  }
  return *this;
}
