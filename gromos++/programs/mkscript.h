// mkscript.h

enum filetype{unknownfile, inputfile, topofile, coordfile, refposfile, 
	      posresspecfile, disresfile, pttopofile, dihresfile, jvaluefile,
	      ledihfile, outputfile, outtrxfile, outtrvfile, outtrefile, outtrgfile, 
	      scriptfile};
int numFiletypes=17;
typedef map<string, filetype>::value_type FT;
const FT filetypes[] ={FT("", unknownfile),
		       FT("input", inputfile),
		       FT("topo", topofile),
		       FT("coord", coordfile),
		       FT("refpos", refposfile),
		       FT("posresspec", posresspecfile),
		       FT("disres", disresfile),
		       FT("pttopo", pttopofile),
		       FT("dihres", dihresfile),
		       FT("jvalue", jvaluefile),
		       FT("ledih", ledihfile),
		       FT("output", outputfile),
		       FT("outtrx", outtrxfile),
		       FT("outtrv", outtrvfile),
		       FT("outtre", outtrefile),
		       FT("outtrg", outtrgfile),
		       FT("script", scriptfile)
		       };
static map<string,filetype> FILETYPE(filetypes,filetypes+numFiletypes);

enum blocktype {unknown, systemblock, startblock,stepblock,boundaryblock,
                submoleculesblock, tcoupleblock, pcoupleblock, 
		centreofmassblock,printblock, 
		writeblock,shakeblock, forceblock,
                plistblock, longrangeblock, posrestblock, perturbblock};
typedef map<string, blocktype>::value_type BT;
const BT blocktypes[] ={BT("",unknown),
			BT("SYSTEM",systemblock),
			BT("START",startblock),
			BT("STEP",stepblock),
			BT("BOUNDARY",boundaryblock),
			BT("SUBMOLECULES",submoleculesblock),
			BT("TCOUPLE",tcoupleblock),
			BT("PCOUPLE",pcoupleblock),
			BT("CENTREOFMASS",centreofmassblock),
			BT("PRINT",printblock),
			BT("WRITE",writeblock),
                        BT("SHAKE",shakeblock),
			BT("FORCE",forceblock),
			BT("PLIST",plistblock),
			BT("LONGRANGE",longrangeblock),
			BT("POSREST",posrestblock),
			BT("PERTURB",perturbblock)
};
static map<string,blocktype> BLOCKTYPE(blocktypes,blocktypes+17);

enum templateelement{unknowntemplate, systemtemplate, numbertemplate,
		     oldnumbertemplate,  
		      start_timetemplate, end_timetemplate, queuetemplate};
typedef map<string, templateelement>::value_type TE;
const TE templateelements[]={TE("", unknowntemplate),
			     TE("system", systemtemplate),
			     TE("number", numbertemplate),
			     TE("oldnumber", oldnumbertemplate),
			     TE("start_time", start_timetemplate),
			     TE("end_time", end_timetemplate),
			     TE("queue", queuetemplate)
};
static map<string,templateelement> TEMPLATE(templateelements, 
					     templateelements+7);

//BLOCKDEFINITIONS
class isystem{
public:
    int found, npm, nsm;
    isystem(){found=0;}
};

class istart{
public:
    int found, ntx,init,ntx0;
    double ig,tempi,heat,boltz;
    istart(){found=0;}
};

class istep{
 public:
    int found, nstlim;
    double t,dt;
    istep(){found=0;}
};

class iboundary{
 public:
    int found, ntb,nrdbox;
    double box[3],beta;
    iboundary(){found=0;}
};

class isubmolecules{
 public:
    int found;
    vector<int> nsp;
    isubmolecules(){found=0;}
};

class itcouple{
 public:
    int found, ntt[3];
    double temp0[3], taut[3];
    itcouple(){found=0;}
};

class ipcouple{
 public:
    int found, ntp;
    double pres0,comp, taup;
    ipcouple(){found=0;}
};

class icentreofmass{
 public:
    int found, ndfmin,ntcm,nscm;
    icentreofmass(){found=0;}
};

class iprint{
 public:
    int found, ntpr, ntpl, ntpp;
    iprint(){found=0;}
};

class iwrite{
 public:
    int found, ntwx,ntwse,ntwv,ntwe,ntwg,ntpw;
    iwrite(){found=0; ntwx=0; ntwse=0; ntwv=0; ntwe=0; ntpw=0;}
};

class ishake{
 public:
    int found, ntc;
    double tol;
    ishake(){found=0;}
};

class iforce{
 public:
    int found, ntf[10];
    vector<int> nre;
    iforce(){found=0;}
};

class iplist{
 public:
    int found, ntnb, nsnb;
    double rcutp, rcutl;
    iplist(){found=0;}
};

class ilongrange{
 public:
    int found;
    double epsrf, appak, rcrf;
    ilongrange(){found=0;}
};

class iposrest{
 public:
    int found, ntr, nrdrx;
    double cho;
    iposrest(){found=0;}
};

class iperturb{
 public:
    int found, ntg,nrdgl,nlam,mmu;
    double rlam, dlamt, rmu, dmut,alphlj,alphc;
    iperturb(){found=0;}
};

class iunknown
{
 public:
  string name;
  string content;
  iunknown(string n):name(n){}
};


class input{
public:
  isystem system;
  istart start;
  istep step;
  iboundary boundary;
  isubmolecules submolecules;
  itcouple tcouple;
  ipcouple pcouple;
  icentreofmass centreofmass;
  iprint print;
  iwrite write;
  ishake shake;
  iforce force;
  iplist plist;
  ilongrange longrange;
  iposrest posrest;
  iperturb perturb;
  vector<iunknown> unknown;
}; 

class fileInfo{
 public:
    double box[3];
    vector<string> blocks;
    vector<int> blockslength;
};

//INSTREAM
istringstream &operator>>(istringstream &is, isystem &s){
    string e;
    s.found=1;
    is >> s.npm >> s.nsm >> e;
    return is;
}
istringstream &operator>>(istringstream &is, istart &s){
    string e;
    s.found=1;
    is >> s.ntx >> s.init >> s.ig >> s.tempi >> s.heat >> s.ntx0 >> s.boltz >> e;
    return is;
}
istringstream &operator>>(istringstream &is, istep &s){
    string e;
    s.found=1;
    is >> s.nstlim >> s.t >> s.dt >> e;
    return is;
}
istringstream &operator>>(istringstream &is, iboundary &s){
    string e;
    s.found=1;
    is >> s.ntb >> s.box[0] >> s.box[1] >> s.box[2] >> s.beta >> s.nrdbox >> e;
    return is;
}
istringstream &operator>>(istringstream &is, isubmolecules &s){
    string e;
    s.found=1;
    int nspm,nsp;
    is >> nspm;
    for(int i=0; i<nspm; i++){is >> nsp; s.nsp.push_back(nsp);}
    is >> e;
    return is;
}
istringstream &operator>>(istringstream &is,itcouple &s){
    string e;
    s.found=1;
    for(int i=0;i<3;i++) is >> s.ntt[i] >> s.temp0[i] >> s.taut[i];
    is >> e;
    return is;
}
istringstream &operator>>(istringstream &is,ipcouple &s){
    string e;
    s.found=1;
    is >> s.ntp >> s.pres0 >> s.comp >> s.taup >> e;
    return is;
}
istringstream &operator>>(istringstream &is,icentreofmass &s){
    string e;
    s.found=1;
    is >> s.ndfmin >> s.ntcm >> s.nscm >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iprint &s){
    string e;
    s.found=1;
    is >> s.ntpr >> s.ntpl >> s.ntpp >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iwrite &s){
    string e;
    s.found=1;
    is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwe >> s.ntwg >> s.ntpw >> e;
    return is;
}
istringstream &operator>>(istringstream &is,ishake &s){
    string e;
    s.found=1;
    is >> s.ntc >> s.tol >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iforce &s){
    string e;
    s.found=1;
    int negr, nre;
    for(int i=0; i<10; i++) is >> s.ntf[i];
    is >> negr;
    for(int i=0; i<negr; i++) {is >> nre; s.nre.push_back(nre);}
    is >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iplist &s){
    string e;
    s.found=1;
    is >> s.ntnb >> s.nsnb >> s.rcutp >> s.rcutl >> e;
    return is;
}
istringstream &operator>>(istringstream &is,ilongrange &s){
    string e;
    s.found=1;
    is >> s.epsrf >> s.appak >> s.rcrf >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iposrest &s){
    string e;
    s.found=1;
    is >> s.ntr >> s.cho >> s.nrdrx >> e;
    return is;
}
istringstream &operator>>(istringstream &is,iperturb &s){
    string e;
    s.found=1;
    is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt >> s.rmu >> s.dmut
       >> s.alphlj >> s.alphc >> s.nlam >> s.mmu >> e;
    return is;
}
istringstream &operator>>(istringstream &is, iunknown &s){
  string e=is.str();
  s.content=e.substr(0,e.find("END"));
  return is;
}


Ginstream &operator>>(Ginstream &is,input &gin){
  vector<string> buffer;
  while (!is.stream().eof()){
    is.getblock(buffer);
    if(buffer.size()){
      
    
      string bufferstring;
      gio::concatenate(buffer.begin()+1, buffer.end()-1, bufferstring);
      istringstream bfstream(bufferstring);
      switch(BLOCKTYPE[buffer[0]]){
	case systemblock:       bfstream >> gin.system;          break;
	case startblock:        bfstream >> gin.start;           break;
	case stepblock:         bfstream >> gin.step;            break;
	case boundaryblock:     bfstream >> gin.boundary;        break;
	case submoleculesblock: bfstream >> gin.submolecules;    break;
	case tcoupleblock:      bfstream >> gin.tcouple;         break;
	case pcoupleblock:      bfstream >> gin.pcouple;         break;
	case centreofmassblock: bfstream >> gin.centreofmass;    break;
	case printblock:        bfstream >> gin.print;           break;
	case writeblock:        bfstream >> gin.write;           break;
	case shakeblock:        bfstream >> gin.shake;           break;
	case forceblock:        bfstream >> gin.force;           break;
	case plistblock:        bfstream >> gin.plist;           break;
	case longrangeblock:    bfstream >> gin.longrange;       break;
	case posrestblock:      bfstream >> gin.posrest;         break;
	case perturbblock:      bfstream >> gin.perturb;         break;
	  
	case unknown:
	  iunknown newblock(buffer[0]);
	  bfstream >> newblock;
	  gin.unknown.push_back(newblock);
	  cout << "Don't know anything about block " << buffer[0]
	       << ". Just storing data.\n";
      }
    }
    
  }
  return is;
}

Ginstream &operator>>(Ginstream &is,fileInfo &s){

    string e;
    string first;
    vector<string> buffer;
    is.getline(first);
    
    while(!is.stream().eof()){
      is.getblock(buffer);
      s.blocks.push_back(first);
      if (first=="BOX") {
	istringstream iss(buffer[0]);
	iss >> s.box[0] >> s.box[1] >> s.box[2] >> e;
      }
      if (first=="TRICLINICBOX"){
	double dummy1, dummy2;
	istringstream iss(buffer[1]);
	iss >> s.box[0];
	iss.str(buffer[2]);
	iss >> dummy1 >> s.box[1];
	iss.str(buffer[3]);
	iss >> dummy1 >> dummy2 >> s.box[2];
      }
      
      s.blockslength.push_back(buffer.size()-1);
      is.getline(first);
    }
    return is;
}

// TEMPLATE handling of (output) filenames

class filename
{
  vector<string> d_parts;
  double d_time, d_dt;
  int d_start;
  string d_system;
  string d_queue;
  string d_template;
  
public:
  filename();
  filename(string s, double t, double dt, int start=1, string q=""){
    d_system=s;
    d_time=t;
    d_dt=dt;
    d_start=start;
    d_queue=q;
  };
  void setInfo(string s, double t, double dt, int start=1, string q=""){
    d_system=s;
    d_time=t;
    d_dt=dt;
    d_start=start;
    d_queue=q;
  };
  
  void setTemplate(string s);
  string temp(){ return d_template;};
  string name(int number);
};

void filename::setTemplate(string s)
{
  d_template=s;
  d_parts.clear();
  string::size_type iter;
  
  string sub;
  iter=s.find('%');
  while(iter !=string::npos){
    sub=s.substr(0,iter);
    s=s.substr(iter+1,s.size()-iter-1);
    iter=s.find('%');
    d_parts.push_back(sub);
    
  }
  d_parts.push_back(s);
}

string filename::name(int number)
{
  ostringstream os;
  for(unsigned int i=0; i<d_parts.size();i++){
    if(i%2){
      switch(TEMPLATE[d_parts[i]]){
	case systemtemplate:     os << d_system;               break;
	case numbertemplate:     os << d_start + number;       break;
	case oldnumbertemplate:  os << d_start + number -1 ;   break;
	case start_timetemplate: os << d_time+number*d_dt; break;
	case end_timetemplate:   os << d_time+(number+1)*d_dt;     break;
        case queuetemplate:      os << d_queue;                break;
	case unknowntemplate:
	  cout << "Do not know how to handle " << d_parts[i] 
	       << " in template. Just printing the words." << endl;
	  os << d_parts[i];
	  break;
      }
    }
   else os << d_parts[i];
  }
  return os.str();
}

// Jobinfo
class jobinfo
{
public:
  map<string, string> param;
  string dir;
  int prev_id;
};


// Writing out of an input file
ostream &operator<<(ostream &os, input &gin)
{
  if(gin.system.found)
    os << "SYSTEM\n"
       << "#      NPM      NSM\n"
       << setw(10) << gin.system.npm
       << setw(9)  << gin.system.nsm
       << "\nEND\n";
  if(gin.start.found)
    os << "START\n"
       << "# values for INIT -- starting procedures in RUNMD.f\n"
       << "#  INIT\n"
       << "#     shake X\n"
       << "#            shake V\n"
       << "# centre of mass motion removal if NTCM=1\n"
       << "#     1   yes    yes   yes\n"
       << "#     2    no    yes   yes\n"
       << "#     3    no     no   yes\n"
       << "#     4    no     no    no\n"
       << "#      NTX     INIT      IG     TEMPI      HEAT  NTXO   BOLTZ\n"
       << setw(10) << gin.start.ntx
       << setw(9)  << gin.start.init
       << setw(8)  << gin.start.ig
       << setw(10)  << gin.start.tempi
       << setw(10) << gin.start.heat
       << setw(6)  << gin.start.ntx0
       << setw(12)  << gin.start.boltz
       << "\nEND\n";
  if(gin.step.found)
    os << "STEP\n"
       << "#   NSTLIM         T        DT\n"
       << setw(10) << gin.step.nstlim
       << setw(10) << gin.step.t
       << setw(10) << gin.step.dt
       << "\nEND\n";
  if(gin.boundary.found)
    os << "BOUNDARY\n"
       << "#      NTB    BOX(1)    BOX(2)    BOX(3)      BETA  NRDBOX\n"
       << setw(10) << gin.boundary.ntb
       << setw(10) << gin.boundary.box[0]
       << setw(10) << gin.boundary.box[1] 
       << setw(10) << gin.boundary.box[2]
       << setw(10) << gin.boundary.beta
       << setw(8) << gin.boundary.nrdbox
       << "\nEND\n";
  if(gin.submolecules.found){
    os << "SUBMOLECULES\n"
       << "#     NSPM  NSP(1.. NSPM)\n"
       << setw(10) << gin.submolecules.nsp.size() << "\n";
    for(unsigned int i=0; i< gin.submolecules.nsp.size(); i++){
      os << setw(6) << gin.submolecules.nsp[i];
      if((i+1)%10==0) os << "\n";
    }
    os << "\nEND\n";
  }
  if(gin.tcouple.found)
    os << "TCOUPLE\n"
       << "#      NTT     TEMP0      TAUT\n"
       << setw(10) << gin.tcouple.ntt[0] 
       << setw(10) << gin.tcouple.temp0[0] 
       << setw(10) << gin.tcouple.taut[0] << "\n"
       << setw(10) << gin.tcouple.ntt[1] 
       << setw(10) << gin.tcouple.temp0[1] 
       << setw(10) << gin.tcouple.taut[1] << "\n"
       << setw(10) << gin.tcouple.ntt[2] 
       << setw(10) << gin.tcouple.temp0[2] 
       << setw(10) << gin.tcouple.taut[2]
       << "\nEND\n";
  if(gin.pcouple.found)
    os << "PCOUPLE\n"
       << "#      NTP     PRES0      COMP      TAUP\n"
       << setw(10) << gin.pcouple.ntp
       << setw(10) << gin.pcouple.pres0
       << setw(10) << gin.pcouple.comp
       << setw(10) << gin.pcouple.taup
       << "\nEND\n";
  if(gin.centreofmass.found)
    os << "CENTREOFMASS\n"
       << "#   NDFMIN      NTCM      NSCM\n"
       << setw(10) << gin.centreofmass.ndfmin
       << setw(10) << gin.centreofmass.ntcm
       << setw(10) << gin.centreofmass.nscm
       << "\nEND\n";
  if(gin.print.found)
    os << "PRINT\n"
       << "#NTPR: print out energies, etc. every NTPR steps\n"
       << "#NTPL: print out C.O.M motion and total energy fitting every NTPL steps\n"
       << "#NTPP: =1 perform dihedral angle transition monitoring\n"
       << "#     NTPR      NTPL      NTPP\n"
       << setw(10) << gin.print.ntpr
       << setw(10) << gin.print.ntpl
       << setw(10) << gin.print.ntpp
       << "\nEND\n";
  if(gin.write.found)
    os << "WRITE\n"
       << "# NTPW = 0 : binary\n"
       << "# NTPW = 1 : formatted\n"
       << "# NTWSE = configuration selection parameter\n"
       << "# =0: write normal trajectory\n"
       << "# >0: chose min energy for writing configurations\n"
       << "#     NTWX     NTWSE      NTWV      NTWE      NTWG      NTPW\n"
       << setw(10) << gin.write.ntwx 
       << setw(10) << gin.write.ntwse
       << setw(10) << gin.write.ntwv 
       << setw(10) << gin.write.ntwe 
       << setw(10) << gin.write.ntwg 
       << setw(10) << gin.write.ntpw
       << "\nEND\n";
  if(gin.shake.found)
    os << "SHAKE\n"
       << "#      NTC       TOL\n"
       << setw(10) << gin.shake.ntc
       << setw(10) << gin.shake.tol
       << "\nEND\n";
  if(gin.force.found){
    os << "FORCE\n"
       << "#      NTF array\n"
       << "# bonds    angles   imp.     dihe     charge nonbonded\n"
       << "# H        H        H        H\n"
       << setw(3) << gin.force.ntf[0] << setw(3) << gin.force.ntf[1] 
       << setw(6) << gin.force.ntf[2] << setw(3) << gin.force.ntf[3]
       << setw(6) << gin.force.ntf[4] << setw(3) << gin.force.ntf[5]
       << setw(6) << gin.force.ntf[6] << setw(3) << gin.force.ntf[7]
       << setw(6) << gin.force.ntf[8] << setw(3) << gin.force.ntf[9]
       << "\n# NEGR    NRE(1)    NRE(2)    ...      NRE(NEGR)\n"
       << setw(6) << gin.force.nre.size();
    for(unsigned int i=0; i< gin.force.nre.size(); i++)
      os << setw(10) << gin.force.nre[i];
    os << "\nEND\n";
  }
  if(gin.plist.found)
    os << "PLIST\n"
       << "#     NTNB      NSNB     RCUTP     RCUTL\n"
       << setw(10) << gin.plist.ntnb
       << setw(10) << gin.plist.nsnb
       << setw(10) << gin.plist.rcutp
       << setw(10) << gin.plist.rcutl
       << "\nEND\n";
  if(gin.longrange.found)
    os << "LONGRANGE\n"
       << "# EPSRF     APPAK      RCRF\n"
       << setw(7) << gin.longrange.epsrf
       << setw(10) << gin.longrange.appak
       << setw(10) << gin.longrange.rcrf
       << "\nEND\n";
  if(gin.posrest.found)
    os << "POSREST\n"
       << "#     values for NTR\n"
       << "#     0: no position re(con)straining\n"
       << "#     1: use CHO\n"
       << "#     2: use CHO/ ATOMIC B-FACTORS\n"
       << "#     3: position constraining\n"
       << "#      NTR       CHO   NRDRX\n"
       << setw(10) << gin.posrest.ntr
       << setw(10) << gin.posrest.cho
       << setw(10) << gin.posrest.nrdrx
       << "\nEND\n";
  if(gin.perturb.found)
    os << "PERTURB\n"
       << "# NTG: 0 no perturbation is applied\n"
       << "#    : 1 calculate dV/dRLAM perturbation\n"
       << "#    : 2 calculate dV/dRMU perturbation\n"
       << "#    : 3 calculate both derivatives\n"
       << "#      NTG     NRDGL     RLAM     DLAMT        RMU   DMUT\n"
       << setw(10) << gin.perturb.ntg
       << setw(10) << gin.perturb.nrdgl
       << setw(9) << gin.perturb.rlam
       << setw(10) << gin.perturb.dlamt
       << setw(11) << gin.perturb.rmu
       << setw(7) << gin.perturb.dmut
       << "\n#   ALPHLJ     ALPHC     NLAM       MMU\n"
       << setw(10) << gin.perturb.alphlj
       << setw(10) << gin.perturb.alphc
       << setw(9) << gin.perturb.nlam
       << setw(10) << gin.perturb.mmu
       << "\nEND\n";

  for(unsigned int i=0; i< gin.unknown.size(); i++){
    os << gin.unknown[i].name << "\n"
       << gin.unknown[i].content
       << "END\n";
  }
  
    
  return os;
  
}
