#include "Arguments.h"
#include <algorithm>
#include <sstream>
#include <fstream>
#include <set>
#include <string>

using namespace std;

namespace args{

typedef multimap<string,string>::value_type argType;

// checks if an argument is known
static int isKnown(const string str, set<string> known_args){
  set<string>::const_iterator fr, to=known_args.end();
  for(fr=known_args.begin();fr!=to; fr++)
    if(*fr==str)break;
  if(fr==to)return 0;
  return 1;
}

class Arguments_i{
  friend class Arguments;
  friend istream &operator>>(istream &is, Arguments &args);
  string d_usage;
  string d_prog;
  set<string> d_known;
  Arguments_i():d_usage(""), d_prog(""), d_known(){}
};

Arguments::Arguments(int argc, char **argv, int nknown, 
		     char **known, const string &usage) :
  multimap<string,string>(), d_this(new Arguments_i) 
{

  if(argc)d_this->d_prog=argv[0];
  d_this->d_usage="\n\n"+usage;

  for(int i=0;i<nknown;++i)
    d_this->d_known.insert(string(known[i]));

  string s("");
  
  for(int i=1;i<argc;i++){
    
    if(string(argv[i])=="@f"){
      // input file
      ++i;
      ifstream inp(argv[i]);
      inp>>*this;
      inp.close();
    }
    else
      s+=string(argv[i])+' ';
  }
  
  // istrstream is(s.c_str());
  stringstream is(s.c_str());
  is>>*this;
}

Arguments::~Arguments() {
  delete d_this;
}

istream &args::operator>>(istream &istr, Arguments &args)
{
  // get away the comments
  char buff[1000];
  string s("");
  
  while(istr.good()&&istr.getline(buff,1000)){
    s+=string(buff);
    if(s.find("#")<=s.length())s=s.erase(s.find("#"));
    s+='\n';
  }
  stringstream is(s.c_str());
  

  string str, last;
  
  if(!(is>>last))
    return istr;
  if(last[0]!='@')
    throw Arguments::Exception(args.d_this->d_usage);
  last=last.substr(1);
  if(args.find(last)!=args.end())
    args.erase(args.lower_bound(last), args.upper_bound(last));

  if(!isKnown(last, args.d_this->d_known)) {
    string except = "\n\nArgument @" + last  + " not known! Possible arguments: " + args.d_this->d_usage;
    throw Arguments::Exception(except);
  }
  
  while(is>>str){
    if(str=="@help") throw Arguments::Exception(args.d_this->d_usage);
      
    if(str[0]=='@'){
      if(args.find(last)==args.end())
	args.insert(argType(last,""));
      last=str.substr(1);
      if(args.find(last)!=args.end())
	for(Arguments::iterator l=args.lower_bound(last);
	    l!=args.upper_bound(last);++l)
	  args.erase(l);

      if(!isKnown(last, args.d_this->d_known)) {
        string except = "\n\nArgument @" + last  + " not known! Possible arguments: " + args.d_this->d_usage;	
	throw Arguments::Exception(except);
      
      continue;
      }
    }
    else
      args.insert(argType(last,str));
  }

  // insert the last one without value
  if(args.find(last)==args.end())
    args.insert(argType(last,""));

  
  return istr;
}
  


const string &Arguments::operator[](const string &str)const
{
  const_iterator f=find(str);
  if(f==end()){
    throw Exception(d_this->d_usage);
  }
  
  else return find(str)->second;
}

int Arguments::check(const string &str, int num_args)const
{
  if(find(str)==end())
    throw Exception(d_this->d_usage);
  int num=0;
  for(const_iterator l=lower_bound(str), u=upper_bound(str);
      l!=u;++l)
    if (l->second!="")++num;
  if(num<num_args)
    throw Exception(d_this->d_usage);
  return 0;
}

int Arguments::count(const string &str)const
{
  if(find(str)==end())
    return -1;
  int num=0;
  for(const_iterator l=lower_bound(str), u=upper_bound(str);
      l!=u;++l)
    if (l->second!="")++num;
  return num;
}

}
