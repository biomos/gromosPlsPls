/// gio_InParameter.cc

#include "InParameter.h"
#include "Ginstream.h"
#include "../gcore/BondType.h"
#include "../gcore/Bond.h"
#include "../gcore/AngleType.h"
#include "../gcore/Angle.h"
#include "../gcore/DihedralType.h"
#include "../gcore/Dihedral.h"
#include "../gcore/ImproperType.h"
#include "../gcore/Improper.h"
#include "../gcore/LJType.h"
#include "../gcore/AtomPair.h"
#include "../gcore/Exclusion.h"
#include "../gcore/AtomTopology.h"
#include "../gcore/MoleculeTopology.h"
#include "../gcore/Molecule.h"
#include "../gcore/System.h"
#include "../gcore/GromosForceField.h"

#include <map>
#include <deque>
#include <set>
#include <sstream>

using namespace std;
using namespace gcore;
using gio::InParameter_i;
using gio::InParameter;

// Implementation class
class InParameter_i{
  friend class gio::InParameter;
  Ginstream d_gin;
  GromosForceField d_gff;
  string d_version;
  string d_name;
  int rdline(char s[7][15], int num);
  void init();
  InParameter_i (const char *name):
    d_gin(name), d_version(), d_name(name){
    init();
  }
  ~InParameter_i(){
    d_gin.close();
  }
};

// Constructors

InParameter::InParameter(string name):
  d_this(new InParameter_i(name.c_str())){
}

InParameter::~InParameter(){
  delete d_this;
}

const string &InParameter::title()const{
  return d_this->d_gin.title();
}

const GromosForceField &InParameter::forceField()const{
  return d_this->d_gff;
}
int InParameter_i::rdline(char s[7][15],int num){
  // just a little roiutine to read a line, write num strings in s
  // and return 1 if the line reads "END"
  // should probably move to Ginstream, but the s[7][15] definition
  // should then be generalised. (strings can't be used with atof etc.
  string str;
  
  d_gin.getline(str,100);
  if(str=="END") return 1;
  else{
    for(int l=0; l<num;l++){
      unsigned int k=0;
      while((str[0]!=' '&&str[0]!='\t')&&str!="") {
	s[l][k++]=str[0];
	str.erase(str.begin());
      }
      while((str[0]==' '||str[0]=='\t')&&str!="") str.erase(str.begin());
      s[l][k]='\0';
    }
    return 0;
  }
}

void InParameter_i::init(){

  if(!d_gin)
    throw InParameter::Exception("Could not open parameter file "+d_gin.name()+".");

  // generic variables
  double d[4];
  int num,i[5];
  string s;
  //char ss[7][15];

  // MASSATOMTYPECODE block
  if(! d_gin.check("MASSATOMTYPECODE"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo MASSATOMTYPECODE block!");
  d_gin >> s;
  while(s!="END"){
    
   // i[0]=atoi((char *)s.begin());
    i[0]=atoi(s.c_str());
    // istringstream i(&s.begin());
    d_gin >> d[0];
    d_gin >> s;
    
    //cout << i[0] << "\t"<<d[0]<<"\t"<<s<<endl;
    // maybe do something with this knowledge?
    d_gin >> s;
  }
   
  // BONDTYPECODE block
  if(! d_gin.check("BONDTYPECODE"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo BONDTYPECODE block!");
  d_gin >> s;
  while(s!="END"){
    d_gin >> d[0];
    d_gin >> d[1];
    d_gff.addBondType(BondType(d[0],d[1]));
    
    d_gin >> s;
    
  }
  

  // BONDANGLETYPECOD block
  if(! d_gin.check("BONDANGLETYPECOD"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo BONDANGLETYPECOD block!");
  d_gin >> s;
  
  while(s!="END"){
    d_gin >> d[0];
    d_gin >> d[1];
    d_gff.addAngleType(AngleType(d[0],d[1]));

    d_gin >> s;
  }
  
  // IMPDIHEDRALTYPEC block
  if(! d_gin.check("IMPDIHEDRALTYPEC"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo IMPDIHEDRALTYPEC block!");
  d_gin >> s;
  while(s!="END"){
    d_gin >> d[0];
    d_gin >> d[1];
    d_gff.addImproperType(ImproperType(d[0],d[1]));
    
    d_gin >> s;
  }
  
  // DIHEDRALTYPECODE block
  if(! d_gin.check("DIHEDRALTYPECODE"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo DIHEDRALTYPECODE block!");
  d_gin >> s;
  while(s!="END"){
    d_gin >> d[0];
    d_gin >> d[1];
    d_gin >> i[0];
    d_gff.addDihedralType(DihedralType(d[0],d[1],i[0]));
    
    d_gin >> s;
  }

  // SINGLEATOMLJPAIR block
  if(! d_gin.check("SINGLEATOMLJPAIR"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo SINGLEATOMLJPAIR block!");
  d_gin >> num;
  double sc6[num], sc12[3][num], scs6[num], scs12[num];
  int    pl[num][num];
  for(int j=0;j<num; j++){
    d_gin >> i[0] >> s >> sc6[j] >> sc12[0][j] >> sc12[1][j] >> sc12[2][j];
    d_gin >> scs6[j] >> scs12[j];
    
    for(int k=0; k<num; k++)
      d_gin >> pl[j][k];
    
    d_gff.addAtomTypeName(s);
    for(int k=0; k<=j;k++){
      d[1]=sc6[j]*sc6[k];
      d[0]=sc12[pl[j][k]-1][j]*sc12[pl[k][j]-1][k];
      d[3]=scs6[j]*scs6[k];
      d[2]=scs12[j]*scs12[k];
      
      d_gff.setLJType(AtomPair(j,k),LJType(d[0],d[1],d[2],d[3]));
    }
    
  }
  if(! d_gin.check())
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nSINGLEATOMLJPAIR block is not OK!");
  // MIXEDATOMLJPAIR block
  if(! d_gin.check("MIXEDATOMLJPAIR"))
    throw InParameter::Exception("Parameter file "+d_gin.name()+" is corrupted:\nNo MIXEDATOMLJPAIR block!");
  d_gin >> s;
  while(s!="END"){
  //  i[0]=atoi((char *)s.begin());
   i[0]=atoi(s.c_str());
    d_gin >> i[1];
    d_gin >> d[1];
    d_gin >> d[0];
    d_gin >> d[3];
    d_gin >> d[2];
    d_gff.setLJType(AtomPair(--i[0],--i[1]),LJType(d[0],d[1],d[2],d[3]));
    
    d_gin >> s;
  }
}




