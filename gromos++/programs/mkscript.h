// mkscript.h

enum filetype{unknownfile, inputfile, topofile, coordfile, refposfile, 
	      posresspecfile, disresfile, pttopofile, dihresfile, jvaluefile,
	      ledihfile};
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
		       FT("ledih", ledihfile)
		       };
static map<string,filetype> FILETYPE(filetypes,filetypes+11);

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
    iwrite(){found=0;}
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
}; 

class fileInfo{
 public:
    double box[3];
    vector<string> blocks;
    vector<int> blockslength;
};

//INSTREAM
Ginstream &operator>>(Ginstream &is, isystem &s){
    string e;
    s.found=1;
    is >> s.npm >> s.nsm >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is, istart &s){
    string e;
    s.found=1;
    is >> s.ntx >> s.init >> s.ig >> s.tempi >> s.heat >> s.ntx0 >> s.boltz >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is, istep &s){
    string e;
    s.found=1;
    is >> s.nstlim >> s.t >> s.dt >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is, iboundary &s){
    string e;
    s.found=1;
    is >> s.ntb >> s.box[0] >> s.box[1] >> s.box[2] >> s.beta >> s.nrdbox >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is, isubmolecules &s){
    string e;
    s.found=1;
    int nspm,nsp;
    is >> nspm;
    for(int i=0; i<nspm; i++){is >> nsp; s.nsp.push_back(nsp);}
    is >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,itcouple &s){
    string e;
    s.found=1;
    for(int i=0;i<3;i++) is >> s.ntt[i] >> s.temp0[i] >> s.taut[i];
    is >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,ipcouple &s){
    string e;
    s.found=1;
    is >> s.ntp >> s.pres0 >> s.comp >> s.taup >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,icentreofmass &s){
    string e;
    s.found=1;
    is >> s.ndfmin >> s.ntcm >> s.nscm >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iprint &s){
    string e;
    s.found=1;
    is >> s.ntpr >> s.ntpl >> s.ntpp >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iwrite &s){
    string e;
    s.found=1;
    is >> s.ntwx >> s.ntwse >> s.ntwv >> s.ntwe >> s.ntwg >> s.ntpw >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,ishake &s){
    string e;
    s.found=1;
    is >> s.ntc >> s.tol >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iforce &s){
    string e;
    s.found=1;
    int negr, nre;
    for(int i=0; i<10; i++) is >> s.ntf[i];
    is >> negr;
    for(int i=0; i<negr; i++) {is >> nre; s.nre.push_back(nre);}
    is >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iplist &s){
    string e;
    s.found=1;
    is >> s.ntnb >> s.nsnb >> s.rcutp >> s.rcutl >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,ilongrange &s){
    string e;
    s.found=1;
    is >> s.epsrf >> s.appak >> s.rcrf >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iposrest &s){
    string e;
    s.found=1;
    is >> s.ntr >> s.cho >> s.nrdrx >> e;
    return is;
}
Ginstream &operator>>(Ginstream &is,iperturb &s){
    string e;
    s.found=1;
    is >> s.ntg >> s.nrdgl >> s.rlam >> s.dlamt >> s.rmu >> s.dmut
       >> s.alphlj >> s.alphc >> s.nlam >> s.mmu >> e;
    return is;
}

Ginstream &operator>>(Ginstream &is,input &gin){
    string s;
    while (!is.eof()){
	is >> s;
        switch(BLOCKTYPE[s]){
	    case systemblock:       is >> gin.system;          break;
	    case startblock:        is >> gin.start;           break;
            case stepblock:         is >> gin.step;            break;
            case boundaryblock:     is >> gin.boundary;        break;
            case submoleculesblock: is >> gin.submolecules;    break;
            case tcoupleblock:      is >> gin.tcouple;         break;
            case pcoupleblock:      is >> gin.pcouple;         break;
            case centreofmassblock: is >> gin.centreofmass;    break;
            case printblock:        is >> gin.print;           break;
            case writeblock:        is >> gin.write;           break;
            case shakeblock:        is >> gin.shake;           break;
            case forceblock:        is >> gin.force;           break;
            case plistblock:        is >> gin.plist;           break;
            case longrangeblock:    is >> gin.longrange;       break;
            case posrestblock:      is >> gin.posrest;         break;
            case perturbblock:      is >> gin.perturb;         break;

            case unknown: 
		cout << "Don't know anything about block " << s
		     << ". Skipping.\n";
		while(s!="END"){is >> s;}
	}
    }
    return is;
}

Ginstream &operator>>(Ginstream &is,fileInfo &s){

    string e;
    int idum=0;
    while(!is.eof()){
	idum=0;
	is >> e;
	s.blocks.push_back(e);
	if (e=="BOX") {
	    is >> s.box[0] >> s.box[1] >> s.box[2] >> e;
	    idum=2;
	}
	else{
	    while(e!="END"){
		is.getline(e,100);
		idum++;
	    }
	}
	s.blockslength.push_back(--idum);
    }
    return is;
}

