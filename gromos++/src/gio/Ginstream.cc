// gio_Ginstream.cc

#include "Ginstream.h"
#include "../gmath/Vec.h"

using namespace std;
using gio::Ginstream;

// Max. characters / line
const int Ginstream::MAXCHAR=1000;

Ginstream::Ginstream() : ifstream(),d_title(), d_name(){}

Ginstream::Ginstream(const string &name, ios::openmode mode):
  ifstream(name.c_str(), mode), d_title(), d_name(name)
{
  if(!good())
    throw Exception("Could not open file "+ d_name +".");

  readTitle();
}

Ginstream::~Ginstream(){
  close();
}

void Ginstream::open(const char *name, ios::openmode  mode){
  d_name=string(name);
  
  ifstream::open(name, mode);
  if(!good())
    throw Exception("Could not open file " + d_name +".");
  readTitle();
}

void Ginstream::open(const string &name, ios::openmode  mode)
{
  open(name.c_str(),mode);
}


Ginstream &Ginstream::getline(string &str, int max)
{
  char *buffer=new char[MAXCHAR];

  do{
    if (good()) {
      ifstream::getline(buffer,max);
    } 
    else {
      if(eof())throw Exception("Early end of file in " + name() +".\n");
      else throw Exception("Error reading " + name() +".\n");
    }
    str=string(buffer);
    if(str.find("#")<=str.length())str=str.erase(str.find("#"));
    while((str[0]==' '||str[0]=='\t')&&str!="")str.erase(str.begin());
    
  }while(str=="");
  
  return *this;
}


int Ginstream::eof()
{
  int flag=1;
  char a;
  char buffer[Ginstream::MAXCHAR];

  do {
    a = peek();
    switch (a) {
    case '#': ifstream::getline(buffer,Ginstream::MAXCHAR);break;
    case ' ': get();break;
    case '\n': get();break;
    case '\t': get();break;
    case -1: return 1;
    default: flag=0;
    }
  } while (flag);

  return ifstream::eof();
}

Ginstream &Ginstream::operator>>(Vec &v)
{
  for (int i=0;i<3;i++){
    if (!good()) return *this;
    *this>>v[i];
  }
  return *this;
}

Ginstream &Ginstream::operator>>(char s[]) 
{
  char a;
  int flag = 1;
  char buff[Ginstream::MAXCHAR];
  
  if (ifstream::good()) {
    do {
      a = peek();

      if (!good()) return *this;
      
      switch (a) {
      case '#': ifstream::getline(buff,Ginstream::MAXCHAR);break;
      case ' ': get();break;
      case '\n': get();break;
      case '\t': get();break;
      default: flag=0;
      }
      
    } while (flag);
    
    // Read in the actual string
    std::operator>>(*this,s);
    
    return *this;
  }
  else{
    if(eof())throw Ginstream::Exception("Early end of file in >>" + name()+ ".\n");
    else throw Ginstream::Exception("Failed to read " + name()+ ".\n");
  }
  
  return *this;
}

Ginstream &Ginstream::operator>>(string &s) 
{
  char buff[Ginstream::MAXCHAR];
  *this >> buff;
  s=string(buff);
  return *this;
}


Ginstream &Ginstream::operator>>(int &c) 
{
  char buffer[Ginstream::MAXCHAR];

  char a = peek();
  while(a==' '||a=='\n'||a=='\t'){get();a=peek();}
  
  while(peek()=='#') 
    ifstream::getline(buffer,Ginstream::MAXCHAR);
  
  ifstream::operator>>(c);
  return *this;
}

Ginstream &Ginstream::operator>>(double &d) 
{
  char buffer[Ginstream::MAXCHAR];

  char a = peek();
  while(a==' '||a=='\n'||a=='\t'){get();a=peek();}
  
  while(peek()=='#') 
    ifstream::getline(buffer,Ginstream::MAXCHAR);
  
  ifstream::operator>>(d);
  return *this;
}

Ginstream &Ginstream::operator>>(float &d) 
{
  char buffer[Ginstream::MAXCHAR];

  char a = peek();
  while(a==' '||a=='\n'||a=='\t'){get();a=peek();}
  
  while(peek()=='#') 
    ifstream::getline(buffer,Ginstream::MAXCHAR);
  
  ifstream::operator>>(d);
  return *this;
}

bool Ginstream::check(const string &str){
  string line;
    
    *this >> line;
  return (line==str);
}

void Ginstream::readTitle()
{
  // Reading of title block
  if(!check("TITLE"))
    throw Exception("Corrupt Gromos file: "+name()+": No TITLE block.");
  
  d_title="";
  string l;
  getline(l);
  while(l!="END") {
    d_title+=l+"\n";
    getline(l);
  }
  d_title=d_title.substr(0,d_title.size()-1);
}

const string &Ginstream::name()const{
  return d_name;
}

void Ginstream::close(){
  ifstream::close();
  d_title="";
  d_name="";
}

const string &Ginstream::title()const{
  return d_title;
}
